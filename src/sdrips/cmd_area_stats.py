import os
import logging
from typing import List, Optional, Tuple
from pathlib import Path
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.mask import mask as riomask
from ruamel.yaml import YAML
from tqdm import tqdm
import multiprocessing

from sdrips.utils.utils import (
    read_region_shapefile, 
    reproject_shapefile,
    get_valid_raster_path,
    get_mean_value_precip,
    initialize_cmd_area_df,
    compute_masked_mean,
    load_yaml_config,
    get_cmd_area_list
)


def process_week(save_data_loc: str, cmd_area_gdf: gpd.GeoDataFrame, week: str, cmd_area_list: list, feature_name: str, precip_condition: bool, percolation_condition: bool, pbar: tqdm) -> gpd.GeoDataFrame:
    """
    Processes command area statistics for a single week.

    Args:
        cmd_area_gdf (gpd.GeoDataFrame): GeoDataFrame with command areas.
        week (str): Week label (e.g., 'currentweek', 'lastweek').
        pbar (tqdm): Progress bar instance.

    Returns:
        gpd.GeoDataFrame: Updated GeoDataFrame with ET and irrigation data.
    """
    is_current = week == 'currentweek'
    is_last = week == 'lastweek'
    day_label = '7Day' if is_current else '14Day'

    cmd_area_gdf[f'{day_label} Irrigation'] = np.nan
    cmd_area_gdf[f'{day_label} Penman ET'] = np.nan
    cmd_area_gdf[f'{day_label} SEBAL ET'] = np.nan

    for i, (cmd_name, *_) in enumerate(cmd_area_list):
        try:
            geometry = cmd_area_gdf.geometry
            penman_path = f"{save_data_loc}/landsat/penman/{week}/penman_eto_{cmd_name.replace(' ', '_').replace('-', '_')}.constant.tif"
            sebal_path = f"{save_data_loc}/landsat/sebal/{week}/sebal_eto_{cmd_name.replace(' ', '_').replace('-', '_')}.eta.tif"

            cmd_area_gdf.loc[i, f'{day_label} Penman ET'] = compute_masked_mean(penman_path, geometry)
            cmd_area_gdf.loc[i, f'{day_label} SEBAL ET'] = compute_masked_mean(sebal_path, geometry)

            if precip_condition and (is_current or is_last):
                ppt_curr_path = f"{save_data_loc}/precip/precip.currentweek.tif"
                ppt_next_path = f"{save_data_loc}/precip/precip.nextweek.tif"
                
                with rio.open(ppt_curr_path) as r:
                    cmd_area_gdf.loc[i, 'Currentweek PPT'] = get_mean_value_precip(r, geometry[i:i+1])
                with rio.open(ppt_next_path) as r:
                    cmd_area_gdf.loc[i, 'Nextweek PPT'] = get_mean_value_precip(r, geometry[i:i+1])
            pbar.update(1)
        except Exception as e:
            logging.error(f"Error processing {cmd_name} for {week}: {e}")
            continue

    cmd_area_gdf[f'{day_label} Irrigation'] = (
        cmd_area_gdf[f'{day_label} SEBAL ET'] - cmd_area_gdf[f'{day_label} Penman ET']
    )

    if percolation_condition:
        percolation_df = pd.read_csv(f'{save_data_loc}/percolation/percolation_{week}.csv')
        cmd_area_gdf = cmd_area_gdf.merge(percolation_df, on=feature_name)

    if is_current or is_last:
        irrigation = cmd_area_gdf[f'{day_label} Irrigation']
        ppt = cmd_area_gdf.get('Currentweek PPT', 0) + cmd_area_gdf.get('Nextweek PPT', 0)
        perc = cmd_area_gdf['MedianPercolation'] * 7 if percolation_condition else 0
        cmd_area_gdf['net_water_req'] = irrigation + ppt - perc

    else:
        cmd_area_gdf['net_water_req'] = cmd_area_gdf[f'{day_label} Irrigation']

    return cmd_area_gdf


def command_area_info(config_path: Path) -> None:
    """
    Orchestrates the generation of irrigation and water requirement statistics
    for all command areas and saves the results as CSV files.
    Args:
        config_path (str): Path to the main configuration file.
    """
    logger = logging.getLogger()
    logger.critical('Started Command Area Info')

    yaml = YAML()
    yaml.preserve_quotes = True
    script_config = load_yaml_config(config_path)
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    irrigation_canals_path = script_config['Irrigation_cmd_area_shapefile']['path']
    feature_name = script_config['Irrigation_cmd_area_shapefile']['feature_name']
    area_column_name = script_config['Irrigation_cmd_area_shapefile']['area_column_name']
    numeric_ID = script_config['Irrigation_cmd_area_shapefile']['numeric_id_name']
    cmd_area_list = get_cmd_area_list(config_path)
    run_week = script_config['Date_Running']['run_week']
    precip_condition = script_config['Precipitation_Config']['consider_preciptation']
    percolation_condition = script_config['Percolation_Config']['consider_percolation']


    cmd_area_gdf = read_region_shapefile(irrigation_canals_path)
    rasterpath = get_valid_raster_path(save_data_loc, cmd_area_list, run_week)

    if not rasterpath:
        raise FileNotFoundError("No valid raster found for CRS extraction.")

    with rio.open(rasterpath) as raster:
        raster_crs = raster.crs
    cmd_area_gdf = reproject_shapefile(cmd_area_gdf, raster_crs)

    cmd_area_gdf['Currentweek PPT'] = np.nan
    cmd_area_gdf['Nextweek PPT'] = np.nan
    cmd_area_gdf = initialize_cmd_area_df(cmd_area_gdf, cmd_area_list, feature_name)

    total_cmd_areas = (
        len(cmd_area_list) + 1 if set(run_week) == {'currentweek', 'lastweek'}
        else len(cmd_area_list) * len(run_week)
    )

    with tqdm(total=total_cmd_areas, desc="Processing Command Area Stats", unit=" Command Areas") as pbar:
        for week in run_week:
            cmd_area_gdf = process_week(save_data_loc = save_data_loc, cmd_area_gdf = cmd_area_gdf, cmd_area_list = cmd_area_list, week = week, feature_name = feature_name, precip_condition = precip_condition, percolation_condition = percolation_condition, pbar = pbar)
            # cmd_area_gdf = cmd_area_gdf.drop(columns=['geometry'])
            cmd_area_gdf['net_water_req_mm'] = cmd_area_gdf['net_water_req'].apply(lambda x: x if x <= 0 else 0) 
            if area_column_name in cmd_area_gdf.columns:
                cmd_area_gdf['net_water_req_m3'] = abs(cmd_area_gdf['net_water_req_mm']/1000) * (cmd_area_gdf[area_column_name]) 
                sorted_cmd_area_gdf = cmd_area_gdf.sort_values(by=f'{numeric_ID}')
                sorted_cmd_area_gdf = sorted_cmd_area_gdf.drop(columns=['geometry'])
                if week == 'currentweek':
                    sorted_cmd_area_gdf.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats.csv",index = False)
                    logging.critical('Finished Command Area Info')
                else:
                    sorted_cmd_area_gdf.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats_{week}.csv",index = False)
                    logging.critical('Finished Last Week Command Area Info')
            else:
                logging.info(f"Area column '{area_column_name}' not found in command area GeoDataFrame. Estimating area from geometry, this might be different than the irrigable area of the command area.")
                if cmd_area_gdf.crs.is_geographic:
                    logging.info("Command area shapefile is in geographic CRS. Estimating UTM CRS for area calculation.")
                    utm_crs = cmd_area_gdf.estimate_utm_crs()
                    cmd_area_gdf = cmd_area_gdf.to_crs(utm_crs)
                cmd_area_gdf['Area_m2'] = cmd_area_gdf.geometry.area
                cmd_area_gdf['net_water_req_m3'] = abs(cmd_area_gdf['net_water_req_mm']/1000) * (cmd_area_gdf['Area_m2']) 
                sorted_cmd_area_gdf = cmd_area_gdf.sort_values(by=f'{numeric_ID}')
                sorted_cmd_area_gdf = sorted_cmd_area_gdf.drop(columns=['geometry'])
                if week == 'currentweek':
                    sorted_cmd_area_gdf.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats.csv",index = False)
                    logging.critical('Finished Command Area Info')
                else:
                    sorted_cmd_area_gdf.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats_{week}.csv",index = False)
                    logging.critical('Finished Last Week Command Area Info')