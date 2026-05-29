import datetime
import logging
import os
import sys
import traceback
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed

import ee
import pandas as pd
from tqdm import tqdm

from sdrips.utils.utils import (
    load_yaml_config,
    get_irrigation_cmd_area
)
from sdrips.utils.ee_utils import ensure_ee_initialized


def estimate_region_soil_moisture(region, s1_collection, feature_name,
                                  depth_band="0-5cm", depth_mm=100) -> pd.DataFrame:
    """
    Estimate soil moisture (volumetric θ and mm) and percolation (mm) for a command area.

    Workflow (aligned with the validated GEE script):
      1) Sentinel-1 SM_index:
         - Smooth VV using focal_max(30 m)
         - wet = max(Smooth) over time window, dry = min(Smooth) over time window
         - SM_index = (Smooth - dry) / (wet - dry)
         - composite SM_index via median over time window
         - clamp to 0..1
      2) θ_res ≈ wilting point (1500 kPa) from SoilGrids official:
         ISRIC/SoilGrids250m/v2_0/wv1500  (val_<depth>_mean * 1e-3)
      3) θ_fc from SoilGrids official:
         ISRIC/SoilGrids250m/v2_0/wv0033  (val_<depth>_mean * 1e-3)
      4) θ_sat derived from bulk density (bdod) from community asset:
         projects/soilgrids-isric/bdod_mean  (bdod_<depth>_mean)
         θ_sat = 1 - ρb / ρs, ρs = 2650 kg/m3
         bdod assumed in cg/cm3 -> /100 => g/cm3 -> *1000 => kg/m3
      5) Scale SM_index -> θ:
         θ = SM_index * (θ_sat - θ_res) + θ_res
         θ bounded in [θ_res, θ_sat]
      6) Convert to mm with assumed depth (depth_mm):
         SM_mm = θ * depth_mm
         FC_mm = θ_fc * depth_mm
         Percolation_mm = max(0, SM_mm - FC_mm)

    Returns a one-row DataFrame for this region.
    """
    logger = logging.getLogger()
    try:
        region_feature = ee.Feature(region)
        region_id = region_feature.get(feature_name).getInfo()
        region_geometry = region_feature.geometry()

        # -----------------------------
        # 1) Sentinel-1 SM_index (sDRIPS-style)
        # -----------------------------
        def add_smooth(img):
            return img.focal_max(30, "circle", "meters").rename("Smooth") \
                     .copyProperties(img, img.propertyNames())

        s1_smooth = s1_collection.map(add_smooth)  # single-band "Smooth" images

        wet = s1_smooth.max().rename("wet_index")
        dry = s1_smooth.min().rename("dry_index")
        denom = wet.subtract(dry)

        # Avoid division by ~0; also mask invalid pixels
        valid = denom.abs().gt(1e-6)
        denom_safe = denom.where(denom.eq(0), 1e-6)

        sm_index = (
            s1_smooth.median()
            .subtract(dry)
            .divide(denom_safe)
            .rename("SM_index")
            .updateMask(valid)
            .clamp(0, 0.6)
        )

        # -----------------------------
        # 2) SoilGrids θ_res (≈ WP), θ_fc (FC) from official catalog
        # -----------------------------
        depth_band_res = depth_band.replace("-", "_")
        theta_res = (
            ee.Image("ISRIC/SoilGrids250m/v2_0/wv1500")
            .select(f"val_{depth_band_res}_mean")
            .multiply(1e-3)
            .rename("theta_res")
        )

        theta_fc = (
            ee.Image("ISRIC/SoilGrids250m/v2_0/wv0033")
            .select(f"val_{depth_band_res}_mean")
            .multiply(1e-3)
            .rename("theta_fc")
        )

        # -----------------------------
        # 3) θ_sat from bdod (community asset)
        # -----------------------------
        bdod = (
            ee.Image("projects/soilgrids-isric/bdod_mean")
            .select(f"bdod_{depth_band}_mean")
            .rename("bdod")
        )

        rho_b = bdod.divide(100).multiply(1000)  # cg/cm3 -> g/cm3 -> kg/m3
        rho_s = ee.Image.constant(2650)
        theta_sat = (
            ee.Image.constant(1)
            .subtract(rho_b.divide(rho_s))
            .rename("theta_sat")
            .clamp(0.05, 0.75)
        )

        # -----------------------------
        # 4) Scale to volumetric θ
        # θ = SM_index * (θ_sat - θ_res) + θ_res
        # -----------------------------
        theta = (
            sm_index.multiply(theta_sat.subtract(theta_res))
            .add(theta_res)
            .rename("theta")
        )
        theta = theta.max(theta_res).min(theta_sat).rename("theta")

        # -----------------------------
        # 5) Convert to mm and percolation (mm)
        # -----------------------------
        depth_img = ee.Image.constant(depth_mm)
        sm_mm = theta.multiply(depth_img).rename("sm_mm")
        fc_mm = theta_fc.multiply(depth_img).rename("fc_mm")
        perc_mm = sm_mm.subtract(fc_mm).max(ee.Image.constant(0)).rename("perc_mm")

        # -----------------------------
        # 6) Reduce over region (median; you can change to mean if desired)
        # -----------------------------
        reducer = ee.Reducer.median()

        out_dict = ee.Image.cat([theta, theta_fc, sm_mm, fc_mm, perc_mm]).reduceRegion(
            reducer=reducer,
            geometry=region_geometry,
            scale=250,
            maxPixels=1e9
        )

        # pull values client-side
        vals = out_dict.getInfo()

        # Return one row per region (same structure style as your existing code)
        return pd.DataFrame([{
            feature_name: region_id,
            "MedianTheta": vals.get("theta"),
            "MedianThetaFC": vals.get("theta_fc"),
            "MedianSoilMoisture_mm": vals.get("sm_mm"),
            "MedianFieldCapacity_mm": vals.get("fc_mm"),
            "MedianPercolation_mm": vals.get("perc_mm"),
        }])

    except Exception as error:
        logger.error(f"Region ID: {region.get('properties', {}).get(feature_name)} failed. Error: {error}")
        return pd.DataFrame()


def percolation_estimation(config_path) -> None:
    """
    Keeps your original structure (parallel per-region), but updates:
      - Soil moisture index (SM_index) to match the validated GEE logic,
      - Converts to volumetric θ (m3/m3),
      - Assumes depth=100 mm and computes percolation = max(0, SM_mm - FC_mm),
      - Uses SoilGrids wv0033/wv1500 and bdod for θ_fc, θ_res, θ_sat.
    """
    logger = logging.getLogger()
    logger.critical("Started Percolation Estimation")

    script_config = load_yaml_config(config_path)
    secrets_file_path = script_config["Secrets_Path"]["path"]
    secrets = load_yaml_config(rf"{secrets_file_path}")

    gee_service_acc = secrets["GEE_Account"]["username"]
    gee_key_file = secrets["GEE_Account"]["key_file"]
    ensure_ee_initialized(service_account=gee_service_acc, key_file=gee_key_file)

    feature_name = script_config["Irrigation_cmd_area_shapefile"]["feature_name"]

    cores = script_config.get("Multiprocessing", {}).get("cores")
    worker_count = cores if cores is not None else max(1, multiprocessing.cpu_count() - 1)

    irrigation_cmd_area = get_irrigation_cmd_area(config_path)

    start_date = script_config["Date_Running"]["start_date"]
    run_week = script_config["Date_Running"]["run_week"]
    save_data_loc = script_config["Save_Data_Location"]["save_data_loc"]

    # optional config knobs (safe defaults)
    depth_band = script_config.get("SoilGrids", {}).get("depth_band", "0-5cm")
    depth_mm = script_config.get("SoilWater", {}).get("depth_mm", 100)

    for wktime in run_week:
        try:
            if wktime == "currentweek":
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days=10)
            elif wktime == "lastweek":
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days=17)
            else:
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days=10)

            startDate = ee.Date(startdate)
            endDate = ee.Date(enddate)

            logger.critical(f"Running Week: {wktime}")
            logger.critical(f"Start Date: {startdate}")
            logger.critical(f"End Date: {enddate}")
            logger.critical(f"SoilGrids depth_band: {depth_band}, depth_mm: {depth_mm}")

            # Sentinel-1 collection for the window
            s1_collection = (
                ee.ImageCollection("COPERNICUS/S1_GRD")
                .filterBounds(irrigation_cmd_area)
                .filterDate(startDate, endDate)
                .filter(ee.Filter.eq("instrumentMode", "IW"))
                .filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
                .select(["VV"])
            )

            region_list = irrigation_cmd_area.toList(irrigation_cmd_area.size()).getInfo()
            region_dataframes = []

            with ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(
                        estimate_region_soil_moisture,
                        region,
                        s1_collection,
                        feature_name,
                        depth_band,
                        depth_mm
                    ): region for region in region_list
                }

                for future in tqdm(
                    as_completed(futures),
                    total=len(futures),
                    desc="Estimating Soil Moisture + Percolation",
                    unit=" Command Areas"
                ):
                    df = future.result()
                    if not df.empty:
                        region_dataframes.append(df)

            final_df = pd.concat(region_dataframes, ignore_index=True)

            # In this updated approach, each region already returns one row.
            # But to keep your original behavior robust, we can still groupby+mean.
            final_means_df = final_df.groupby(feature_name).mean(numeric_only=True).reset_index()

            os.makedirs(f"{save_data_loc}/percolation", exist_ok=True)
            final_means_df.to_csv(f"{save_data_loc}/percolation/Percolation_{wktime}.csv", index=False)

        except Exception:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            logger.error(
                f"Week: {wktime} encountered an error.\n"
                f"File: {fname}, Line: {exc_tb.tb_lineno}\n"
                f"{traceback.format_exc()}"
            )
            continue

    logger.critical("Percolation Estimation Completed")
