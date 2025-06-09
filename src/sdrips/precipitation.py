"""
Download and preprocess IMERG precipitation data.

This script fetches daily IMERG precipitation GeoTIFFs using IMERG credentials in secrets,
clips them to a given bounding box, and resamples them using Rasterio.
"""

import logging
import os
import datetime
import requests
import geopandas as gpd
import rasterio as rio
from rasterio.enums import Resampling
from rasterio.warp import calculate_default_transform, reproject
from rasterio.transform import from_origin
from shapely.geometry import box
import rasterio.mask
import numpy as np
from ruamel.yaml import YAML
from pathlib import Path
from tqdm import tqdm
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed

from sdrips.utils.utils import load_yaml_config







def download_imerg_file(date: datetime.date, save_data_loc, imerg_username, imerg_password ) -> Path:
    """
    Download IMERG GeoTIFF file.

    Args:
        date (datetime.date): Date to download.
        save_data_loc (str): Main directory where the sDRIPS is running.
        imerg_username (str): IMERG username for authentication.
        imerg_password (str): IMERG password for authentication.

    Returns:
        Path: Path to downloaded file.
    """
    logger = logging.getLogger()
    datestr = date.strftime('%Y%m%d')
    filename = f'global.precip.imerg.{datestr}.tif'
    out_path = os.path.join(save_data_loc, 'precip', filename)

    year, month = date.year, date.month

    if year > 2024:
        url = (
            f'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/'
            f'{year}/{month:02d}/3B-HHR-E.MS.MRG.3IMERG.{datestr}-S233000-E235959.1410.V07B.1day.tif'
        )
    elif year == 2024:
         url = (
            f'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/'
            f'{year}/{month:02d}/3B-HHR-E.MS.MRG.3IMERG.{datestr}-S233000-E235959.1410.V07B.1day.tif'
        )
    else:
        url = (
            f'https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/'
            f'{year}/{month:02d}/3B-HHR-E.MS.MRG.3IMERG.{datestr}-S233000-E235959.1410.V07B.1day.tif'
        )

    logger.info(f'Downloading {url}')
    response = requests.get(url, auth=(imerg_username, imerg_password))
    response.raise_for_status()

    with open(out_path, 'wb') as f:
        f.write(response.content)

    return out_path


def clip_and_resample_tif(src_path: Path, dst_path: Path, bounds: tuple, resolution: float = 0.1):
    """
    Clip and resample a GeoTIFF using Rasterio.
    - If bounding box < native resolution (0.1°), use all_touched mask clipping.
    - Else, use reproject with calculate_default_transform.

    Args:
        src_path (Path): Path to input GeoTIFF.
        dst_path (Path): Path to output GeoTIFF.
        bounds (tuple): (left, bottom, right, top) bounding box.
        resolution (float): Output resolution in degrees.
    """
    logger = logging.getLogger()
    left, bottom, right, top = bounds
    bbox_width = abs(right - left)
    bbox_height = abs(top - bottom)

    with rio.open(src_path) as src:
        if bbox_width < resolution or bbox_height < resolution:
            logger.info(f"Boudning box {bounds} is smaller than the IMERG's native resolution {resolution}. Using all_touched mask clipping.")
            geom = [box(left, bottom, right, top).__geo_interface__]
            out_image, out_transform = rasterio.mask.mask(
                src, geom, crop=True, all_touched=True, filled=True
            )
            out_meta = src.meta.copy()
            out_meta.update({
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
                "compress": "lzw"
            })
            with rio.open(dst_path, "w", **out_meta) as dest:
                dest.write(out_image)
        else:
            dst_transform, width, height = calculate_default_transform(
                src.crs, src.crs, src.width, src.height, *bounds, resolution=resolution
            )
            kwargs = src.meta.copy()
            kwargs.update({
                'driver': 'GTiff',
                'height': height,
                'width': width,
                'transform': dst_transform,
                'compress': 'lzw'
            })

            with rio.open(dst_path, 'w', **kwargs) as dst:
                for i in range(1, src.count + 1):
                    reproject(
                        source=rio.band(src, i),
                        destination=rio.band(dst, i),
                        src_transform=src.transform,
                        src_crs=src.crs,
                        dst_transform=dst_transform,
                        dst_crs=src.crs,
                        resampling=Resampling.bilinear
                    )

def process_day(date: datetime.date, bounds: tuple, save_data_loc: str, imerg_username: str, imerg_password: str) -> bool:
    """
    Download, clip, and resample IMERG precipitation data for a given day.

    This function downloads the raw IMERG precipitation GeoTIFF file for the specified
    date, clips it to the specified geographic bounds, resamples it if needed, and
    saves the processed file in the provided output location.

    Args:
        date (datetime.date): The date for which to process the IMERG data.
        bounds (tuple): The bounding box to clip the data (minx, miny, maxx, maxy).
        save_data_loc (str): The directory path where the processed data will be saved.
        imerg_username (str): IMERG username for authentication.
        imerg_password (str): IMERG password for authentication.

    Returns:
        bool: True if processing is successful, False otherwise.
    """

    logger = logging.getLogger()
    try:
        raw_tif = download_imerg_file(date, save_data_loc, imerg_username, imerg_password)
        clipped_tif = os.path.join(save_data_loc, 'precip', f'precip.imerg.{date.strftime("%Y%m%d")}.tif')
        clip_and_resample_tif(raw_tif, clipped_tif, bounds)
        logger.info(f'Processed {clipped_tif}')
        return True
    except Exception as e:
        logger.error(f"Failed to process IMERG for {date}: {e}")
        return False

def compute_weekly_cummulative_precip(file_paths:list[Path], output_path: Path) -> Path:
    """
    Compute the weekly cumulative of IMERG rasters and save to output_path.
    
    Args:
        file_paths (list of str): List of file paths to daily GeoTIFF precipitation files.
        output_path (str): File path to save the weekly averaged GeoTIFF.

    Returns:
        str: The file path of the saved weekly averaged raster.
    """
    arrays = []
    meta = None

    for i, file in enumerate(file_paths):
        with rasterio.open(file) as src:
            data = src.read(1)
            data[data < 0] = 0
            arrays.append(data)

            if meta is None:
                meta = src.meta.copy()

    stack = np.stack(arrays)
    weekly_cummulative = np.sum(stack, axis=0) / 10  

    meta.update(dtype=rasterio.float32)

    with rasterio.open(output_path, "w", **meta) as dst:
        dst.write(weekly_cummulative.astype(rasterio.float32), 1)

def imergprecip(config_path: Path):
    """
    Download and process IMERG precipitation data for 7 previous days in parallel.

    Args:
        config_path (Path): Path to the main config file.
    """
    logger = logging.getLogger()
    logger.critical('IMERG Precipitation Download Started')

    yaml = YAML()
    yaml.preserve_quotes = True
    script_config = load_yaml_config(config_path)

    bounds = script_config['Irrigation_cmd_area_shapefile_Bounds']
    left, right = float(bounds['leftlon']), float(bounds['rightlon'])
    top, bottom = float(bounds['toplat']), float(bounds['bottomlat'])
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    start_date_str = script_config['Date_Running']['start_date']
    start_date = datetime.datetime.strptime(start_date_str, '%Y-%m-%d').date()
    default_run_week = script_config['Date_Running']['default_run_week']
    run_week = (
        ["lastweek", "currentweek"]
        if default_run_week
        else script_config['Date_Running'].get('run_week', [])
    )
    cores = script_config.get("Multiprocessing", {}).get("cores")
    worker_count = cores if cores is not None else multiprocessing.cpu_count() - 1

    secrets_file_path = script_config['Secrets_Path']['path']
    secrets = load_yaml_config(rf'{secrets_file_path}')
    imerg_username = secrets['IMERG_Account']['username']
    imerg_password = secrets['IMERG_Account']['password']


    today = datetime.date.today()
    day_diff = abs((start_date - today).days)
    days_back = 7

    all_dates = [
        (start_date - datetime.timedelta(days=day_offset) if day_diff > 7
         else today - datetime.timedelta(days=day_offset))
        for day_offset in range(1, days_back + 1)
    ]

    success_flags = []
    processed_files = [] 

    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = {executor.submit(process_day, date, (left, bottom, right, top), save_data_loc, imerg_username, imerg_password): date for date in all_dates}
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing IMERG Precipitation Days"):
            success = future.result()
            success_flags.append(success)
            if success:
                date = futures[future]
                processed_file = os.path.join(save_data_loc, 'precip', f'precip.imerg.{date.strftime("%Y%m%d")}.tif')
                processed_files.append(processed_file)

    if all(success_flags):
        logger.critical('Finished IMERG precipitation data download successfully.')

        weekly_output_path = os.path.join(save_data_loc, 'precip', f'precip.currentweek.tif')
        
        compute_weekly_cummulative_precip(processed_files, weekly_output_path)
        logger.critical(f'Weekly cumulative precipitation raster saved to: {weekly_output_path}')
    else:
        logger.critical('Error occurred during IMERG precipitation data download or processing.')
