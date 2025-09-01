
import os
import sys
import logging
import datetime
import geopandas as gpd
import rasterio as rio
import ee
import geemap
import numpy as np
from rasterio.mask import mask as riomask
from ruamel.yaml import YAML
from typing import Tuple, Any, Dict, List, Union, Optional
from tqdm import tqdm
import requests
from pathlib import Path

from sdrips.utils.ee_utils import upload_shapefile_to_ee

yaml = YAML(typ='safe')
yaml.preserve_quotes = True

def get_gdrive_url(file_id: str) -> str:
    """
    Generate direct download URL for Google Drive file
    
    Args:
        file_id (str): Google Drive file ID.

    Returns:
        str: Direct download URL for the file.
    """
    return f"https://drive.google.com/uc?export=download&id={file_id}"


def load_yaml_config(config_path: str) -> Dict:
    """
    Load YAML configuration file.

    Args:
        config_path (str): Path to YAML config file.

    Returns:
        dict: Parsed configuration dictionary.
    """
    with open(config_path, 'r') as f:
        return yaml.load(f)


def read_cmd_area_settings(cmd_area_data: Dict[str, Any], defaults: Dict[str, Any]) -> Dict[str, Any]:
    """
    Read settings for a specific command area with fallback to defaults.

    Args:
        cmd_area_data (Dict[str, Any]): Specific command area configuration.
        defaults (Dict[str, Any]): Default values for planting_date, crop_type, soil_coef.

    Returns:
        Dict[str, Any]: Final resolved settings for the command area.
    """
    if cmd_area_data.get('use_default', False):
        return {
            'planting_date': defaults['planting_date'],
            'crop_type': defaults['crop_type'],
            'soil_coefficient': defaults['soil_coef']
        }
    else:
        return {
            'planting_date': cmd_area_data.get('planting_date', defaults['planting_date']),
            'crop_type': cmd_area_data.get('crop_type', defaults['crop_type']),
            'soil_coefficient': cmd_area_data.get('soil_coef', defaults['soil_coef'])
        }


def read_crop_coefficients(config_path: str, crop: str) -> Dict[Tuple[int, int], Union[float, Tuple]]:
    """
    Read crop coefficient (Kc) values for a given crop from YAML config.

    Args:
        config_path (str): Path to the YAML file containing crop Kc values.
        crop (str): Crop name (case-insensitive).
        - If crop is 'grass' or 'alfa-alfa', returns constant Kc=1.0.

    Returns:
        Dict[Tuple[int, int], Union[float, Tuple]]: Mapping from day ranges to Kc or interpolation info.
    """
    logger = logging.getLogger()
    crop_config = load_yaml_config(config_path)
    normalized_crop_config = {k.lower(): v for k, v in crop_config.items()}
    crop_key = crop.lower()

    if crop.lower() in ['grass', 'alfa-alfa']:
        logger.info(f"Using constant Kc=1.0 for crop '{crop}'")
        return {(0, 999): 1.0}  
    
    if crop_key not in normalized_crop_config:
        logger.error(f"Crop '{crop}' not found in crop config file: {config_path}")
        sys.exit(1)

    crop_data = normalized_crop_config.get(crop_key, {})
    coefficients = {}

    for key, value in crop_data.items():
        try:
            start_day, end_day = map(int, key.split('-'))
        except ValueError:
            logger.error(f"Invalid day range key '{key}' in crop config for crop '{crop}'")
            sys.exit(1)
        if isinstance(value, list) and value[0] == 'linear':
            # Format: ['linear', start_value, end_value, num_days]
            try:
                coefficients[(start_day, end_day)] = (
                    'linear', float(value[1]), float(value[2]), int(value[3])
                )
            except (ValueError, IndexError):
                logger.error(f"Invalid linear format for range {key} in crop '{crop}'")
                sys.exit(1)
        else:
            try:
                coefficients[(start_day, end_day)] = float(value)
            except ValueError:
                logger.error(f"Non-numeric Kc value for range {key} in crop '{crop}'")
                sys.exit(1)

    return coefficients


def get_growth_kc(coefficients: Dict[Tuple[int, int], Union[float, Tuple]], num_days: int) -> Union[float, None]:
    """
    Get crop coefficient (Kc) based on crop growth stage and number of days since planting.

    Args:
        coefficients (dict): Dictionary mapping day ranges to Kc values or linear tuples.
        num_days (int): Days since planting.

    Returns:
        float or None: Crop coefficient for the given day or None if undefined.
    """
    for day_range, kc in coefficients.items():
        if day_range[0] <= num_days <= day_range[1]:
            if isinstance(kc, tuple) and kc[0] == 'linear':
                _, start_val, end_val, days = kc
                days_into_range = num_days - day_range[0]
                delta_per_day = (end_val - start_val) / days
                return start_val + delta_per_day * days_into_range
            return kc
    return None  # If num_days is outside defined ranges


def read_region_shapefile(path: str, target_crs: str = 'EPSG:4326') -> gpd.GeoDataFrame:
    """
    Read a shapefile and ensure it has a valid coordinate reference system (CRS).

    Args:
        path (str): Path to the shapefile.
        target_crs (str, optional): CRS to assign if the file lacks one. Defaults to 'EPSG:4326'.

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with assigned or existing CRS.
    """
    region_gdf = gpd.read_file(path)
    if region_gdf.crs is None:
        region_gdf.crs = target_crs
    return region_gdf


def reproject_shapefile(shapefile: gpd.GeoDataFrame, target_crs) -> gpd.GeoDataFrame:
    """
    Reproject a GeoDataFrame to the specified CRS.

    Args:
        shapefile (gpd.GeoDataFrame): The GeoDataFrame to be reprojected.
        target_crs (str): The target CRS to reproject to (e.g., 'EPSG:32645').

    Returns:
        gpd.GeoDataFrame: Reprojected GeoDataFrame.
    """
    return shapefile.to_crs(target_crs)


def get_valid_raster_path(save_data_loc, cmd_area_list, run_week) -> Optional[str]:
    """
    Returns the first existing irrigation raster path from the canal list.

    Returns:
        Optional[str]: Valid raster file path or None if not found.
    """
    for cmd_area_name, *_ in cmd_area_list:
        path = os.path.join(
            save_data_loc, 'landsat', 'irrigation', run_week[0],
            f"irrigation_{cmd_area_name.replace(' ', '_').replace('-', '_')}.eta.tif"
        )
        if os.path.exists(path):
            return path
    return None


def initialize_cmd_area_df(cmd_area_gdf: gpd.GeoDataFrame, cmd_area_list: List[List[str]], feature_name: str) -> gpd.GeoDataFrame:
    """
    Initializes the command area GeoDataFrame with sorted index and reset geometry.

    Args:
        cmd_area_gdf (gpd.GeoDataFrame): GeoDataFrame of command areas.
        cmd_area_list (List[List[str]]): List of command area names.
        feature_name (str): Column name to be used as index and identifier.

    Returns:
        gpd.GeoDataFrame: Updated and sorted GeoDataFrame.
    """
    simple_cmd_area_names = [item[0] for item in cmd_area_list]
    cmd_area_gdf.set_index(feature_name, inplace=True)
    cmd_area_gdf = cmd_area_gdf.reindex(simple_cmd_area_names)
    cmd_area_gdf.reset_index(inplace=True)
    return cmd_area_gdf


def compute_masked_mean(raster_path: str, geometry) -> float:
    """
    Computes the mean of raster values masked by the given geometry.

    Args:
        raster_path (str): File path to the raster.
        geometry: GeoJSON-like geometry used for masking.

    Returns:
        float: Mean value of the masked raster data.
    """
    with rio.open(raster_path) as src:
        data = src.read(1, masked=True)
        masked_data, _ = riomask(src, geometry, crop=True)
        return np.mean(data)


def get_mean_value_precip(raster, geometry):
    """
    Computes the mean precipitation value for a given geometry.

    Args:
        raster: Opened rasterio dataset object.
        geometry: GeoJSON-like geometry to apply the mask.

    Returns:
        float: Mean precipitation value.
    """
    masked_data, _ = riomask(raster, geometry, crop=True)
    return np.mean(masked_data)


def get_cmd_area_list(config_path: str) -> List[List[str]]:
    """
    Load command area list from the configuration file.
    This function retrieves the command area names from the GEE asset ID or shapefile specified in the configuration.

    Args:
        config_path (str): Path to the main configuration file.
    Returns:
        List[List[str]]: List of command area names.
    """
    yaml = YAML()
    yaml.preserve_quotes = True
    script_config = load_yaml_config(config_path)
    gee_asset_section = script_config.get('GEE_Asset_ID', {})
    feature_name = script_config['Irrigation_cmd_area_shapefile']['feature_name']
    if gee_asset_section.get('id'): 
        gee_asset_id = gee_asset_section['id']
        irrigation_cmd_area = ee.FeatureCollection(gee_asset_id) 
    elif gee_asset_section.get('shp'):
        gee_asset_id = gee_asset_section['shp']
        irrigation_cmd_area = upload_shapefile_to_ee(gee_asset_id)  
    else:
        gee_asset_id = script_config['Irrigation_cmd_area_shapefile']['path']
        secrets_file_path = script_config['Secrets_Path']['path']
        secrets = load_yaml_config(rf'{secrets_file_path}')
        gee_service_acc = secrets['GEE_Account']['username']
        gee_key_file = secrets['GEE_Account']['key_file']
        irrigation_cmd_area = upload_shapefile_to_ee(gee_asset_id, service_account=gee_service_acc, key_file=gee_key_file) 

    # else:
    #     # raise ValueError(
    #     #     "Configuration error: 'GEE_Asset_ID' must contain either 'id' or 'shp'."
    #     #     "Both are missing or empty."
    #     # )
    #     _, gee_asset_id = upload_shapefile_to_ee(gee_asset_id, dry_run = True)  
    # if gee_asset_section.get('id'):
    #     irrigation_cmd_area = ee.FeatureCollection(gee_asset_id) 
    #     # print(f'GEE Asset ID: {gee_asset_id}')
    # else:
    #     irrigation_cmd_area = upload_shapefile_to_ee(gee_asset_id, dry_run = True)  
    cmd_area_list = irrigation_cmd_area.reduceColumns(ee.Reducer.toList(1), [feature_name]).get('list').getInfo()
    return cmd_area_list

def get_irrigation_cmd_area(config_path: str) -> ee.FeatureCollection:
    """
    Load the irrigation command area from the configuration file.
    This function retrieves the command area as an Earth Engine FeatureCollection

    Args:
        config_path (str): Path to the main configuration file.
        
    Returns:
        ee.FeatureCollection: The irrigation command area as a FeatureCollection.
    """
    yaml = YAML()
    yaml.preserve_quotes = True
    script_config = load_yaml_config(config_path)
    gee_asset_section = script_config.get('GEE_Asset_ID', {})
    feature_name = script_config['Irrigation_cmd_area_shapefile']['feature_name']
    if gee_asset_section.get('id'): 
        gee_asset_id = gee_asset_section['id']
        irrigation_cmd_area = ee.FeatureCollection(gee_asset_id) 
    elif gee_asset_section.get('shp'):
        gee_asset_id = gee_asset_section['shp']
        irrigation_cmd_area= upload_shapefile_to_ee(gee_asset_id, dry_run = True) 
    else:
        gee_asset_id = script_config['Irrigation_cmd_area_shapefile']['path']
        secrets_file_path = script_config['Secrets_Path']['path']
        secrets = load_yaml_config(rf'{secrets_file_path}')
        gee_service_acc = secrets['GEE_Account']['username']
        gee_key_file = secrets['GEE_Account']['key_file']
        irrigation_cmd_area = upload_shapefile_to_ee(gee_asset_id, service_account=gee_service_acc, key_file=gee_key_file) 
    # else:
    #     # shp_path= script_config['Irrigation_cmd_area_shapefile']['path']
        
    #     raise ValueError(
    #         "Configuration error: 'GEE_Asset_ID' must contain either 'id' or 'shp'."
    #         "Both are missing or empty."
    #     )
    #     # logger.info(f"Using local shapefile: {gee_asset_id}")
    # if gee_asset_section.get('id'):
    #     irrigation_cmd_area = ee.FeatureCollection(gee_asset_id) 
    #     # print(f'GEE Asset ID: {gee_asset_id}')
    # else:
    #     irrigation_cmd_area= upload_shapefile_to_ee(gee_asset_id, dry_run = True) 
    return irrigation_cmd_area