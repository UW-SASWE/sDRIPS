import os
import logging
from typing import List, Optional
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio as rio
from rasterio.mask import mask as riomask
from tqdm import tqdm

from src.sdrips.utils.utils import (
    read_region_shapefile, 
    reproject_shapefile,
    get_valid_raster_path,
    get_mean_value_precip,
    initialize_cmd_area_df,
    compute_masked_mean
)


def process_week(save_data_loc: str, cmd_area_gdf: gpd.GeoDataFrame, week: str, cmd_area_list: list, feature_name: str, pbar: tqdm) -> gpd.GeoDataFrame:
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
            penman_path = f"{save_data_loc}/landsat/penman/{week}/penman_eto_{cmd_name.replace(' ', '-')}.constant.tif"
            sebal_path = f"{save_data_loc}/landsat/sebal/{week}/sebal_eto_{cmd_name.replace(' ', '-')}.eta.tif"
            
            cmd_area_gdf.loc[i, f'{day_label} Penman ET'] = compute_masked_mean(penman_path, geometry)
            cmd_area_gdf.loc[i, f'{day_label} SEBAL ET'] = compute_masked_mean(sebal_path, geometry)

            if is_current or is_last:
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

    percolation_df = pd.read_csv(f'{save_data_loc}/percolation/percolation_{week}.csv')
    cmd_area_gdf = cmd_area_gdf.merge(percolation_df, on=feature_name)

    if is_current or is_last:
        cmd_area_gdf['net_water_req'] = (
            cmd_area_gdf['Currentweek PPT'] +
            cmd_area_gdf['Nextweek PPT'] +
            cmd_area_gdf[f'{day_label} Irrigation'] -
            cmd_area_gdf['MedianPercolation'] * 7
        )

    return cmd_area_gdf


def command_area_info(save_data_loc: str, irrigation_canals_path: str, run_week: List[str], cmd_area_list: List[List[str]], feature_name: str, numeric_ID: float) -> None:
    """
    Orchestrates the generation of irrigation and water requirement statistics
    for all command areas and saves the results as CSV files.
    """
    logger = logging.getLogger()
    logger.critical('Started Command Area Info')
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
            cmd_area_gdf = process_week(save_data_loc = save_data_loc, cmd_area_gdf = cmd_area_gdf, cmd_area_list = cmd_area_list, week = week, feature_name = feature_name, pbar = pbar)
            cmd_area_gdf = cmd_area_gdf.drop(columns=['geometry'])
            cmd_area_gdf['net_water_req'] = cmd_area_gdf['net_water_req'].apply(lambda x: x if x <= 0 else 0) 
            sorted_cmd_area_gdf = cmd_area_gdf.sort_values(by=f'{numeric_ID}')
            if week == 'currentweek':
                sorted_cmd_area_gdf.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats.csv",index = False)
                logging.critical('Finished Command Area Info')
            else:
                sorted_cmd_area_gdf.to_csv(f"{save_data_loc}/Landsat_Command_Area_Stats_{week}.csv",index = False)
                logging.critical('Finished Last Week Command Area Info')