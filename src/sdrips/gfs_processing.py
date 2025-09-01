import os
import datetime
import logging
import zipfile
import shutil
import requests
import gzip
import numpy as np
import subprocess as sub
from ruamel.yaml import YAML
from tqdm import tqdm
import ee
import urllib.request
import geopandas as gpd
import rasterio as rio
from rasterio.windows import from_bounds
from rasterio.shutil import copy as rio_copy
import tempfile
from typing import Dict, Union

# from sdrips.utils.ee_utils import initialize_earth_engine
from sdrips.utils.ee_utils import ensure_ee_initialized
from sdrips.utils.utils import load_yaml_config




def download_tif_from_ee(image, path, name, bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat):
    """
    Download a GeoTIFF image from Google Earth Engine.

    Args:
        image (ee.Image): Earth Engine image to download.
        path (str): File path where the image will be saved.
        name (str): Name used in the download request.
        bounds_leftlon (float): Western longitude of bounding box.
        bounds_bottomlat (float): Southern latitude of bounding box.
        bounds_rightlon (float): Eastern longitude of bounding box.
        bounds_toplat (float): Northern latitude of bounding box.
    """
    logger = logging.getLogger()
    params = {
        'name':f'{name}',
        'filePerBand': "false",
        'scale': 25000,
        'crs': 'EPSG:4326',
        'fileFormat': 'GeoTIFF',
        'region': ee.Geometry.Rectangle([bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat])
    }
    url = image.getDownloadURL(params)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    urllib.request.urlretrieve(url, path)
    # logger.info(f'Downloaded Successfully From GEE: {path}')



def convert_to_geotiff(input_path, output_path):
    """
    Convert a raster to GeoTIFF format using rasterio.

    Args:
        input_path (str): Path to the input raster file.
        output_path (str): Path to save the GeoTIFF file.
    """
    logger = logging.getLogger()
    with rio.open(input_path) as src:
        profile = src.profile
        profile.update(driver='GTiff')

        with rio.open(output_path, 'w', **profile) as dst:
            dst.write(src.read())
    logger.info(f'Converted to GeoTIFF: {output_path}')


def gfsdata_ee(save_data_loc, start_date, bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat):
    """
    Download and process GFS forecast data from Google Earth Engine.
    
    Args:
        save_data_loc (str): Path to directory for saving processed files.
        start_date (str): Forecast start date in format 'YYYY-MM-DD'.
        bounds_leftlon (float): Western longitude of bounding box.
        bounds_bottomlat (float): Southern latitude of bounding box.
        bounds_rightlon (float): Eastern longitude of bounding box.
        bounds_toplat (float): Northern latitude of bounding box.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    startDate = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    endDate = startDate + datetime.timedelta(days=1)
    datestr = startDate.strftime("%Y%m%d")

    dataset = ee.ImageCollection('NOAA/GFS0P25').filterDate(
        startDate.strftime('%Y-%m-%d'),
        endDate.strftime('%Y-%m-%d')
    ).filterBounds(
        ee.Geometry.Rectangle([bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat])
    )

    apcp_image = dataset.select('total_precipitation_surface') \
                        .filter(ee.Filter.eq('forecast_hours', 168)).first()

    apcp_path = os.path.join(save_data_loc, "precip", f"precip.gfs.{datestr}.nextweek.zip")
    download_tif_from_ee(apcp_image, apcp_path, f"precip_gfs_{datestr}_nextweek", bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat)

    extract_dir = os.path.join(save_data_loc, "precip")
    with zipfile.ZipFile(apcp_path, "r") as zip_input:
        for member in zip_input.namelist():
            extracted_path = zip_input.extract(member, extract_dir)
            new_path = os.path.join(extract_dir, "precip.nextweek.tif")
            os.rename(extracted_path, new_path)
    os.remove(apcp_path)

    weeks = ['currentweek', 'nextweek']
    # params_ee = ['u_component_of_wind_10m_above_ground', 'v_component_of_wind_10m_above_ground', 'temperature_2m_above_ground']
    params_ee = []
    total_iterations = len(weeks) * len(params_ee) * 14 # 14 comes from the range(12, 169, 12)

    if len(params_ee) == 0:
        return

    with tqdm(total=total_iterations, desc="Downloading GFS Data From GEE", unit="file") as pbar:
        for week in weeks:
            week_date = startDate - datetime.timedelta(days=7) if week == 'currentweek' else startDate
            week_date_end = week_date + datetime.timedelta(days=1)
            dataset = ee.ImageCollection('NOAA/GFS0P25').filterDate(
                week_date.strftime('%Y-%m-%d'),
                week_date_end.strftime('%Y-%m-%d')
            ).filterBounds(
                ee.Geometry.Rectangle([bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat])
            )

            for param in params_ee:
                for hours in range(12, 169, 12):
                    image = dataset.select(param).filter(ee.Filter.eq('forecast_hours', hours)).first()
                    folder = 'ugrd' if 'u_component' in param else 'vgrd' if 'v_component' in param else 'avgt'
                    path = os.path.join(save_data_loc, folder, f"{folder}.{week_date.strftime('%Y%m%d')}.{week}.{str(hours).zfill(3)}.zip")
                    download_tif_from_ee(image, path, f"{folder}_{week_date.strftime('%Y%m%d')}_{week}_{str(hours).zfill(3)}", bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat)

                    with zipfile.ZipFile(path, "r") as zip_input:
                        for member in zip_input.namelist():
                            extracted_path = zip_input.extract(member, os.path.join(save_data_loc, folder))
                            gtif_path = os.path.join(
                                save_data_loc, folder,
                                f"{folder}.gfs.{week_date.strftime('%Y%m%d')}.{week}{str(hours).zfill(3)}.tif"
                            )
                            os.rename(extracted_path, gtif_path)

                    final_path = os.path.join(
                        save_data_loc, folder,
                        f"{folder}.{week}.{str(hours).zfill(3)}.tif"
                    )
                    convert_to_geotiff(gtif_path, final_path)
                    logger.info(f"Saved: {final_path}")
                    os.remove(path)
                    pbar.update(1)



def gfsdata_noaa(save_data_loc, start_date, bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat):
    """
    Download and process GFS forecast data from NOAA server.
    
    Args:
        save_data_loc (str): Path to directory for saving processed files.
        start_date (str): Forecast start date in format 'YYYY-MM-DD'.
        bounds_leftlon (float): Western longitude of bounding box.
        bounds_bottomlat (float): Southern latitude of bounding box.
        bounds_rightlon (float): Eastern longitude of bounding box.
        bounds_toplat (float): Northern latitude of bounding box.
    """

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    start_dt = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    datestr = start_dt.strftime("%Y%m%d")

    logger.info("Starting download of GFS forecast data for date %s", datestr)

    precip_url = (
        f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?"
        f"file=gfs.t00z.pgrb2.0p25.f168&var_APCP=on"
        f"&subregion=&leftlon={bounds_leftlon}&rightlon={bounds_rightlon}"
        f"&toplat={bounds_toplat}&bottomlat={bounds_bottomlat}"
        f"&dir=%2Fgfs.{datestr}%2F00%2Fatmos"
    )

    precip_tif_path = os.path.join(save_data_loc, "precip", f"precip.gfs.{datestr}.nextweek.tif")
    precip_asc_path = os.path.join(save_data_loc, "precip", "precip.gfs.nextweek.asc")
    precip_final_tif = os.path.join(save_data_loc, "precip", "precip.nextweek.tif")

    logger.info("GFS download link: %s", precip_url)
    logger.info("Target file: %s", precip_tif_path)

    os.makedirs(os.path.dirname(precip_tif_path), exist_ok=True)

    try:
        sub.run(['curl', '--output', precip_tif_path, precip_url], check=True)
    except sub.CalledProcessError as err:
        logger.error("Failed to download GFS file: %s", err)
        return

    try:
        with rio.open(precip_tif_path) as src:
            data = src.read(2)  # Read band 2
            profile = src.profile
            profile.update(driver='AAIGrid', count=1)
            with rio.open(precip_asc_path, 'w', **profile) as dst:
                dst.write(data[0], 1)
        logger.info("ASC file created: %s", precip_asc_path)
    except Exception as e:
        logger.error("Error during GeoTIFF to ASC conversion: %s", e)

    try:
        with rio.open(precip_asc_path) as src:
            data = src.read(1)
            profile = src.profile
            profile.update(driver='GTiff', count=1)
            with rio.open(precip_final_tif, 'w', **profile) as dst:
                dst.write(data, 1)
        logger.info("Final GeoTIFF created: %s", precip_final_tif)
    except Exception as e:
        logger.error("Error during ASC to GeoTIFF conversion: %s", e)

    # GFS variables
    weeks = ['currentweek', 'nextweek']
    # params = ['tmax', 'tmin', 'ugrd', 'vgrd']
    params = []
    # param_ids = ['TMAX', 'TMIN', 'UGRD', 'VGRD']
    param_ids = []

    total_steps = len(weeks) * len(params) * len(range(12, 169, 12))
    
    if len(params) == 0:
        return 
    
    with tqdm(total=total_steps, desc="Downloading GFS Data From NOAA", unit="file") as pbar:
        for week in weeks:
            for param, param_id in zip(params, param_ids):
                for hour in range(12, 169, 12):
                    date_shift = -7 if week == 'currentweek' else 0
                    date_used = (start_dt + datetime.timedelta(days=date_shift)).strftime("%Y%m%d")

                    url = (
                        f"https://nomads.ncep.noaa.gov/cgi-bin/filter_gfs_0p25.pl?"
                        f"file=gfs.t00z.pgrb2.0p25.f168&var_{param_id}=on"
                        f"&subregion=&leftlon={bounds_leftlon}&rightlon={bounds_rightlon}"
                        f"&toplat={bounds_toplat}&bottomlat={bounds_bottomlat}"
                        f"&dir=%2Fgfs.{date_used}%2F00%2Fatmos"
                    )

                    out_dir = os.path.join(save_data_loc, param)
                    os.makedirs(out_dir, exist_ok=True)

                    gtif_path = os.path.join(out_dir, f"{param}.{date_used}.{week}.{hour:03}.tif")
                    asc_path = os.path.join(out_dir, f"{param}.{week}.{hour:03}.asc")
                    final_tif = os.path.join(out_dir, f"{param}.{week}.{hour:03}.tif")

                    try:
                        sub.run(['curl', '--output', gtif_path, url], check=True)
                    except sub.CalledProcessError:
                        logger.warning("Fallback to base URL for %s", param)
                        fallback_url = url.replace('%2Fatmos', '')
                        sub.run(['curl', '--output', gtif_path, fallback_url], check=True)

                    try:
                        with rio.open(gtif_path) as src:
                            data = src.read(1)
                            profile = src.profile
                            profile.update(driver='AAIGrid', count=1)
                            with rio.open(asc_path, 'w', **profile) as dst:
                                dst.write(data, 1)

                        with rio.open(asc_path) as src:
                            data = src.read(1)
                            profile = src.profile
                            profile.update(driver='GTiff', count=1)
                            with rio.open(final_tif, 'w', **profile) as dst:
                                dst.write(data, 1)

                        logger.info("Processed: %s", final_tif)
                    except Exception as e:
                        logger.error("Error processing %s at hour %d: %s", param, hour, e)

                    pbar.update(1)

def gfsdata(config_path: str)-> None:
    """
    Download and process GFS forecast data based on the date difference.

    If the date difference is greater than 10 days, data is fetched from Google Earth Engine.
    Otherwise, data is fetched from the NOAA server.

    Args:
    - config_path (str): Path to the main configuration file.
    """
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.critical('GFS Data Download Started')


    yaml = YAML()
    yaml.preserve_quotes = True
    script_config = load_yaml_config(config_path)
    secrets_file_path = script_config['Secrets_Path']['path']
    secrets = load_yaml_config(rf'{secrets_file_path}')
    gee_service_acc = secrets['GEE_Account']['username']
    gee_key_file = secrets['GEE_Account']['key_file']
    ensure_ee_initialized(service_account=gee_service_acc, key_file=gee_key_file)

    bounds_leftlon = script_config['Irrigation_cmd_area_shapefile_Bounds']['leftlon']
    bounds_rightlon = script_config['Irrigation_cmd_area_shapefile_Bounds']['rightlon']
    bounds_toplat = script_config['Irrigation_cmd_area_shapefile_Bounds']['toplat']
    bounds_bottomlat = script_config['Irrigation_cmd_area_shapefile_Bounds']['bottomlat']

    if not all([bounds_leftlon, bounds_rightlon, bounds_toplat, bounds_bottomlat]):
        shapefile_path = script_config['Irrigation_cmd_area_shapefile']['path']
        gdf = gpd.read_file(shapefile_path)

        # Ensure CRS is EPSG:4326 for consistency
        if gdf.crs is None:
            logger.error("Shapefile has no CRS defined. Please set it before using bounds.")
            raise ValueError("Shapefile has no CRS defined. Please set it before using bounds.")
        if gdf.crs.to_epsg() != 4326:
            gdf = gdf.to_crs(epsg=4326)

        minx, miny, maxx, maxy = gdf.total_bounds

        bounds_leftlon = minx
        bounds_rightlon = maxx
        bounds_bottomlat = miny
        bounds_toplat = maxy
    else:
        bounds_leftlon = float(bounds_leftlon)
        bounds_rightlon = float(bounds_rightlon)
        bounds_toplat = float(bounds_toplat)
        bounds_bottomlat = float(bounds_bottomlat)

    start_date = script_config['Date_Running']['start_date']
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    start_date_obj = datetime.datetime.strptime(start_date, "%Y-%m-%d")
    date_difference = (datetime.datetime.today() - start_date_obj).days

    if date_difference > 10:
        logger.info("Switching to Google Earth Engine for GFS Forecast Data Download Due To Ten-Day Difference.")
        gfsdata_ee(save_data_loc, start_date, bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat)
        weeks = ['currentweek', 'nextweek']
        # params = ['avgt', 'ugrd', 'vgrd']
        params = []
        letterstr = "ABCDEFGHIJKLMN"

        if len(params) == 0:
            # logger.warning("No parameters found for GFS data processing.")
            return

        total_iterations = len(weeks) * len(params)
        with tqdm(total=total_iterations, desc="Processing GFS Data", unit=" file") as pbar:
            for week in weeks:
                for param in params:
                    datestr = start_date_obj.strftime("%Y%m%d")
                    if week == 'currentweek':
                        datestr = (start_date_obj - datetime.timedelta(days=7)).strftime("%Y%m%d")

                    arrays = []
                    meta = None

                    for idx, hours in enumerate(range(12, 169, 12)):
                        file_path = os.path.join(save_data_loc, param, f"{param}.{week}.{str(hours).zfill(3)}.tif")
                        if os.path.exists(file_path):
                            with rio.open(file_path) as src:
                                if meta is None:
                                    meta = src.meta.copy()
                                arrays.append(src.read(1))
                        else:
                            logger.warning(f"File not found: {file_path}")

                    if arrays and meta:

                        mean_array = np.mean(arrays, axis=0)

                        meta.update(dtype=rio.float32, count=1)

                        out_path = os.path.join(save_data_loc, param, f"{param}.{week}.tif")
                        with rio.open(out_path, 'w', **meta) as dst:
                            dst.write(mean_array.astype(rio.float32), 1)
                        logger.info(f"Saved averaged file: {out_path}")
                    else:
                        logger.warning(f"No data available for parameter: {param} and week: {week}")

                    pbar.update(1)

            for week in weeks:
                ugrd_path = os.path.join(save_data_loc, 'ugrd', f"ugrd.{week}.tif")
                vgrd_path = os.path.join(save_data_loc, 'vgrd', f"vgrd.{week}.tif")
                wind_path = os.path.join(save_data_loc, 'wind', f"wind.{week}.tif")

                if os.path.exists(ugrd_path) and os.path.exists(vgrd_path):
                    with rio.open(ugrd_path) as ugrd_src, rio.open(vgrd_path) as vgrd_src:
                        ugrd = ugrd_src.read(1)
                        vgrd = vgrd_src.read(1)
                        meta = ugrd_src.meta.copy()

                        wind = np.sqrt(ugrd ** 2 + vgrd ** 2)
                        meta.update(dtype=rio.float32, count=1)

                        os.makedirs(os.path.dirname(wind_path), exist_ok=True)
                        with rio.open(wind_path, 'w', **meta) as dst:
                            dst.write(wind.astype(rio.float32), 1)
                        logger.info(f"Saved wind speed file: {wind_path}")
                else:
                    logger.warning(f"Missing ugrd or vgrd files for week: {week}")
    else:
        gfsdata_noaa(save_data_loc, start_date, bounds_leftlon, bounds_bottomlat, bounds_rightlon, bounds_toplat)

        weeks = ['currentweek', 'nextweek']
        # params = ['tmax', 'tmin', 'ugrd', 'vgrd']
        params = []
        letterstr = "ABCDEFGHIJKLMN"

        total_iterations = len(weeks) * len(params)
        with tqdm(total=total_iterations, desc="Processing GFS Data", unit=" file") as pbar:
            for week in weeks:
                for param in params:
                    datestr = start_date_obj.strftime("%Y%m%d")
                    if week == 'currentweek':
                        datestr = (start_date_obj - datetime.timedelta(days=7)).strftime("%Y%m%d")

                    arrays = []
                    meta = None

                    for idx, hours in enumerate(range(12, 169, 12)):
                        file_path = os.path.join(save_data_loc, param, f"{param}.{week}.{str(hours).zfill(3)}.tif")
                        if os.path.exists(file_path):
                            with rio.open(file_path) as src:
                                if meta is None:
                                    meta = src.meta.copy()
                                arrays.append(src.read(1))
                        else:
                            logger.warning(f"File not found: {file_path}")

                    if arrays and meta:
                        mean_array = np.mean(arrays, axis=0)

                        meta.update(dtype=rio.float32, count=1)

                        out_path = os.path.join(save_data_loc, param, f"{param}.{week}.tif")
                        with rio.open(out_path, 'w', **meta) as dst:
                            dst.write(mean_array.astype(rio.float32), 1)
                        logger.info(f"Saved averaged file: {out_path}")
                    else:
                        logger.warning(f"No data available for parameter: {param} and week: {week}")

                    pbar.update(1)

            for week in weeks:
                tmax_path = os.path.join(save_data_loc, 'tmax', f"tmax.{week}.tif")
                tmin_path = os.path.join(save_data_loc, 'tmin', f"tmin.{week}.tif")
                avgt_path = os.path.join(save_data_loc, 'avgt', f"avgt.{week}.tif")

                if os.path.exists(tmax_path) and os.path.exists(tmin_path):
                    with rio.open(tmax_path) as tmax_src, rio.open(tmin_path) as tmin_src:
                        tmax = tmax_src.read(1)
                        tmin = tmin_src.read(1)
                        meta = tmax_src.meta.copy()

                        avgt = (tmax + tmin) / 2
                        meta.update(dtype=rio.float32, count=1)

                        os.makedirs(os.path.dirname(avgt_path), exist_ok=True)
                        with rio.open(avgt_path, 'w', **meta) as dst:
                            dst.write(avgt.astype(rio.float32), 1)
                        logger.info(f"Saved average temperature file: {avgt_path}")
                else:
                    logger.warning(f"Missing tmax or tmin files for week: {week}")

            for week in weeks:
                ugrd_path = os.path.join(save_data_loc, 'ugrd', f"ugrd.{week}.tif")
                vgrd_path = os.path.join(save_data_loc, 'vgrd', f"vgrd.{week}.tif")
                wind_path = os.path.join(save_data_loc, 'wind', f"wind.{week}.tif")

                if os.path.exists(ugrd_path) and os.path.exists(vgrd_path):
                    with rio.open(ugrd_path) as ugrd_src, rio.open(vgrd_path) as vgrd_src:
                        ugrd = ugrd_src.read(1)
                        vgrd = vgrd_src.read(1)
                        meta = ugrd_src.meta.copy()

                        wind = np.sqrt(ugrd ** 2 + vgrd ** 2)
                        meta.update(dtype=rio.float32, count=1)

                        os.makedirs(os.path.dirname(wind_path), exist_ok=True)
                        with rio.open(wind_path, 'w', **meta) as dst:
                            dst.write(wind.astype(rio.float32), 1)
                        logger.info(f"Saved wind speed file: {wind_path}")
                else:
                    logger.warning(f"Missing ugrd or vgrd files for week: {week}")

    logger.critical('Finished GFS Data Download And Processing')