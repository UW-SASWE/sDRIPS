import os
import logging
import zipfile
from tqdm import tqdm
import rasterio as rio
import numpy as np

from sdrips.utils.utils import (
    load_yaml_config,
    get_cmd_area_list
)


def unzip_tiffs(config_path: str) -> None:
    """
    Unzips all .zip files in specified directories organized by sensor, variable, and week.

    Parameters
    ----------
    save_data_loc : str
        Root directory where zipped TIFFs are stored.
    run_weeks : list of str
        List of week identifiers (e.g., ['lastweek', 'currentweek']).

    Returns
    -------
    None
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.critical("Unzipping TIFF files...")
    script_config = load_yaml_config(config_path)
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    run_weeks = script_config['Date_Running']['run_week']
    sensors = ['landsat']
    variables = ['sebal', 'penman','irrigation'] 
    for sensor in sensors:
        for var in variables:
            for week in run_weeks:
                input_dir = os.path.join(save_data_loc, sensor, var, week)
                files_to_process = [f for f in os.listdir(input_dir) if f.endswith(".zip")]
                for filename in tqdm(files_to_process, desc=f"Unzipping {var} for {week}", unit="file"):
                    filepath = os.path.join(input_dir, filename)
                    try:
                        with zipfile.ZipFile(filepath, "r") as zip_file:
                            zip_file.extractall(input_dir)
                    except Exception as err:
                        logger.error(f"Error unzipping {filepath}")
                        logger.exception(err)
                        continue
    logger.critical("Finished unzipping TIFF files.")


def convert_tiffs(config_path: str) -> None:
    """
    Converts and renames TIFF files for upload, setting NoData value and standardizing names.

    Parameters
    ----------
    save_data_loc : str
        Root directory where unzipped TIFFs are stored.
    run_weeks : list of str
        List of week identifiers (e.g., ['lastweek', 'currentweek']).

    Returns
    -------
    None
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.critical("Converting TIFF files...")
    script_config = load_yaml_config(config_path)
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    run_weeks = script_config['Date_Running']['run_week']
    sensors = ['landsat']
    variables = ['sebal', 'penman', 'irrigation']
    for sensor in sensors:
        for week in run_weeks:
            for var in variables:
                input_dir = os.path.join(save_data_loc, sensor, var, week)
                try:
                    files_to_process = [
                            fn for fn in os.listdir(input_dir) if fn.endswith('.tif')
                        ]
                except FileNotFoundError:
                    logger.warning(f"Directory not found: {input_dir}")
                    continue

                for fn in tqdm(files_to_process, desc=f"Processing {var} for {week}", unit=" file"):
                        in_fp = os.path.join(input_dir, fn)
                        cleaned_name = fn.replace('.constant', '') \
                                     .replace('.eta', '') \
                                     .replace('_eto', '') \
                                     .replace('.etr', '')
                        out_fp = os.path.join(save_data_loc, 'uploads', f"{sensor}_{week}_{cleaned_name}")
                        try:
                            with rio.open(in_fp) as src:
                                data = src.read(1).astype('float32') 
                                profile = src.profile.copy()
                                profile.update({
                                    'dtype': 'float32',
                                    'nodata': np.nan,
                                    'driver': 'GTiff'
                                })
                            with rio.open(out_fp, 'w', **profile) as dst:
                                dst.write(data, 1)
                        except Exception as err:
                            logger.error(f"Error processing file: {in_fp}")
                            logger.exception(err)
                            continue
    logger.critical("Finished converting TIFF files.")


def converting_to_eto(config_path:str) -> None:
    """
    Convert Penman-Monteith TIFFs to ETo-compatible TIFFs for upload, renaming accordingly.
    
    Parameters
    ----------  
    config_path : str
        Path to the main configuration file.

    Returns
    -------
    None
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.critical("Starting TIFF conversion for uploading purpose...")
    script_config = load_yaml_config(config_path)
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    run_week = script_config['Date_Running']['run_week']
    cmd_area_list = get_cmd_area_list(config_path)
    for ca_area in tqdm(cmd_area_list, desc="Processing Command Area TIFF to ETo Format", unit=' Command Area'):
        try:
            region_raw = ca_area[0]
            region_id = region_raw.replace(" ", "_").replace("-", "_")

            for sensor in ['landsat']:
                for week in run_week:
                    input_tif = os.path.join(
                        save_data_loc, 'uploads', f"{sensor}_{week}_penman_{region_id}.tif"
                    )
                    output_tif = os.path.join(
                        save_data_loc, 'uploads', f"{sensor}_{week}_penman_eto_{region_id}.tif"
                    )

                    if not os.path.exists(input_tif):
                        logger.warning(f"Missing input TIFF: {input_tif}")
                        continue

                    with rio.open(input_tif) as src:
                        penman_data = src.read(1).astype('float32')
                        metadata = src.meta.copy()
                        metadata.update({
                            'dtype': 'float32',
                            'nodata': np.nan,
                            'driver': 'GTiff'
                        })

                    with rio.open(output_tif, 'w', **metadata) as dst:
                        dst.write(penman_data, 1)

        except Exception as err:
            logger.error("Error while executing converting_to_eto()")
            logger.exception(err)

    logger.info("Successfully saved all TIFF files in the upload folder.")