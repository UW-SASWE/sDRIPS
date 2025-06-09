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
from sdrips.utils.ee_initialize import initialize_earth_engine


def estimate_region_soil_moisture(region, s1_collection, field_capacity, feature_name) -> pd.DataFrame:
    """
    Estimate median soil moisture and field capacity for a command area.

    Parameters:
        region (dict): A GeoJSON-like dictionary representing the region feature.
        s1_collection (ee.ImageCollection): Sentinel-1 image collection.
        field_capacity (ee.Image): Field capacity image clipped to irrigation area.
        feature_name (str): Feature property to use for region identification.

    Returns:
        dict: Dictionary with region ID, median soil moisture, and field capacity.
    """
    logger = logging.getLogger()
    try:
        region_feature = ee.Feature(region)
        region_id = region_feature.get(feature_name).getInfo()
        region_geometry = region_feature.geometry()

        def process_image(image):
            smoothed = image.addBands(image.focal_max(30, 'circle', 'meters').rename("Smooth"))
            wet_index = s1_collection.max().select('VV')
            dry_index = s1_collection.min().select('VV')
            sensitivity = wet_index.subtract(dry_index)

            urban_mask = smoothed.select('Smooth').gt(-6)
            water_mask = smoothed.select('Smooth').lt(-17)

            mv = smoothed.select("Smooth").subtract(dry_index).divide(sensitivity)
            mv = mv.updateMask(water_mask.Not()).updateMask(urban_mask.Not())
            mv_upscaled = mv.reduceResolution(reducer=ee.Reducer.mean(), bestEffort=True).reproject(crs=mv.projection(), scale=250)
            mv_clamped = mv_upscaled.clamp(0, 0.6)

            median_ssm = mv_clamped.reduceRegion(
                reducer=ee.Reducer.median(),
                geometry=region_geometry,
                scale=250,
                maxPixels=1e9
            ).get('Smooth')

            median_fc = field_capacity.reduceRegion(
                reducer=ee.Reducer.median(),
                geometry=region_geometry,
                scale=250,
                maxPixels=1e9
            ).get('b10')

            return ee.Feature(None, {
                feature_name: region_id,
                'MedianSoilMoisture': median_ssm,
                'MedianFieldCapacity': median_fc
            })

        features = s1_collection.map(process_image).getInfo()
        data = [{
            feature_name: f['properties'][feature_name],
            'MedianSoilMoisture': f['properties']['MedianSoilMoisture'],
            'MedianFieldCapacity': f['properties']['MedianFieldCapacity']
        } for f in features['features']]

        return pd.DataFrame(data)
    
    except Exception as error:
        logger.error(f"Region ID: {region.get('properties', {}).get(feature_name)} failed. Error: {error}")
        return pd.DataFrame()


def percolation_estimation(config_path) -> None:
    """
    Estimate percolation for command areas using Sentinel-1 data and soil field capacity maps.

    Runs percolation estimates for each week type defined in `run_week`,
    using parallel processing to accelerate per-region computations.
    Saves the weekly percolation results as CSV files.
    Args:
        config_path (str): Path to the main configuration file.
    """
    logger = logging.getLogger()
    logger.critical("Started Percolation Estimation")
    script_config = load_yaml_config(config_path)
    feature_name = script_config['Irrigation_cmd_area_shapefile']['feature_name']
    # print(f"Feature Name: {feature_name}")
    cores = script_config.get("Multiprocessing", {}).get("cores")
    worker_count = cores if cores is not None else multiprocessing.cpu_count() - 1
    irrigation_cmd_area = get_irrigation_cmd_area(config_path)
    start_date = script_config['Date_Running']['start_date']
    run_week = script_config['Date_Running']['run_week']
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    for wktime in run_week:
        try:
            if wktime == "currentweek":
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days=10)
            elif wktime == "lastweek":
                startdate = start_date
                enddate = datetime.datetime.strptime(start_date, "%Y-%m-%d") + datetime.timedelta(days=17)

            startDate = ee.Date(startdate)
            endDate = ee.Date(enddate)

            logger.critical(f"Running Week: {wktime}")
            logger.critical(f"Start Date: {startdate}")
            logger.critical(f"End Date: {enddate}")

            field_capacity = (
                ee.Image("OpenLandMap/SOL/SOL_WATERCONTENT-33KPA_USDA-4B1C_M/v01")
                .select('b10')
                .divide(100)
                .clip(irrigation_cmd_area)
            )

            s1_collection = (
                ee.ImageCollection('COPERNICUS/S1_GRD')
                .filterBounds(irrigation_cmd_area)
                .filterDate(startDate, endDate)
                .filter(ee.Filter.eq('instrumentMode', 'IW'))
                .select(['VV'])
            )

            region_list = irrigation_cmd_area.toList(irrigation_cmd_area.size()).getInfo()
            region_dataframes = []


            with ThreadPoolExecutor(max_workers = worker_count) as executor:
                futures = {
                    executor.submit(
                        estimate_region_soil_moisture,
                        region,
                        s1_collection,
                        field_capacity,
                        feature_name
                    ): region for region in region_list
                }

                for future in tqdm(as_completed(futures), total=len(futures), desc="Estimating Soil Moisture", unit = " Command Areas"):
                    df = future.result()
                    if not df.empty:
                        region_dataframes.append(df)

            final_df = pd.concat(region_dataframes, ignore_index=True)
            final_means_df = final_df.groupby(feature_name).mean().reset_index()
            final_means_df['MedianPercolation'] = (
                final_means_df['MedianSoilMoisture'] - final_means_df['MedianFieldCapacity']
            ).clip(lower=0)

            os.makedirs(f'{save_data_loc}/percolation', exist_ok=True)
            final_means_df.to_csv(f'{save_data_loc}/percolation/Percolation_{wktime}.csv', index=False)

        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            logger.error(
                f"Week: {wktime} encountered an error.\n"
                f"File: {fname}, Line: {exc_tb.tb_lineno}\n"
                f"{traceback.format_exc()}"
            )
            continue
    logger.critical("Percolation Estimation Completed")