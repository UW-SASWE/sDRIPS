import logging
import re
import ee
from ee.image import Image
from ee.geometry import Geometry
import os, datetime, requests, zipfile, time, math, glob
import pandas as pd
import numpy as np
import geopandas as gpd
from shapely.geometry import Point
from typing import List, Optional, Pattern, Union, Literal
from pathlib import Path
import urllib.request
from bs4 import BeautifulSoup
from ruamel.yaml import YAML
from filelock import FileLock
import json
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
from pykrige.ok import OrdinaryKriging
from geopy.distance import geodesic
import tempfile
import rasterio as rio
from rasterio.transform import from_origin
import warnings
warnings.filterwarnings("ignore")


def download_sensor_data( sensor_date: datetime.date, base_url: str, sensor_data_path: Path, filename_match_regex: Pattern, date_format: str, save_data_loc: Path) -> bool:
    """
    Download the sensor data file for a given date from the server.

    Args:
        sensor_date (datetime.date): Date for which to download the sensor file.
        base_url (str): URL to the sensor data directory.
        sensor_data_path (Path): Local directory to save the file.
        filename_match_regex (Pattern): Regex pattern to find the correct file for the date.
        date_format (str): Date format string (e.g., '%Y_%m_%d') for comparison.
        save_data_loc (Path): Main project working directory.

    Returns:
        bool: True if the file was downloaded successfully, False otherwise.
    """
    logger = logging.getLogger()
    logger.info(f"Downloading sensor data for date: {sensor_date}")
    date_str = sensor_date.strftime("%Y-%m-%d")

    path_of_folder = os.path.join(save_data_loc, sensor_data_path)
    os.makedirs(path_of_folder, exist_ok=True)

    status_json = os.path.join(path_of_folder, "download_log.json")
    status_lock = os.path.join(path_of_folder, "download_log.lock")

    # Step 1: Check download log with lock
    with FileLock(status_lock):
        download_status = {}
        if os.path.exists(status_json):
            with open(status_json, "r") as f:
                download_status = json.load(f)
        if download_status.get(date_str, False):
            logger.info(f"Skipping download for {date_str}; already downloaded this run.")
            return False

    # Step 2: Find the correct file on the server
    try:
        response = requests.get(base_url)
        response.raise_for_status()
    except requests.RequestException as e:
        logger.error(f"Failed to access the given server URL: {e}")
        return False

    soup = BeautifulSoup(response.content, 'html.parser')
    sensor_file_date = sensor_date.strftime(date_format)

    sensor_url = None
    for link in soup.find_all('a'):
        href = link.get('href')
        if href and filename_match_regex.search(href) and sensor_file_date in href:
            sensor_url = base_url + href
            break

    if not sensor_url:
        logger.critical(f"No file found for date: {sensor_file_date}")
        return False

    # Step 3: Download and save the file
    try:
        response = requests.get(sensor_url)
        response.raise_for_status()
    except requests.RequestException as e:
        logger.error(f"Failed to download sensor file: {e}")
        return False

    file_path = os.path.join(path_of_folder, 'Sensor_Data.csv')
    try:
        with open(file_path, "wb") as f:
            f.write(response.content)
        logger.info(f"Sensor data downloaded successfully to {file_path}")
    except IOError as e:
        logger.error(f"Failed to save sensor file: {e}")
        return False

    # Step 4: Update download status 
    with FileLock(status_lock):
        if os.path.exists(status_json):
            with open(status_json, "r") as f:
                download_status = json.load(f)
        else:
            download_status = {}
        download_status[date_str] = True
        with open(status_json, "w") as f:
            json.dump(download_status, f)
        logger.info(f"Updated download log file at {status_json}")

    return True


def get_available_sensor_dates(base_url: str, date_regex: Pattern, date_format: str) -> List[datetime.date]:
    """
    Fetch a sorted list of available sensor dates from the given base URL
    by extracting dates from filenames using a regular expression.

    Args:
        base_url (str): URL to the sensor data directory on the server.
        date_regex (Pattern): Compiled regular expression to extract date (YYYY_MM_DD).
        date_format (str): Format to parse date strings (e.g., '%Y_%m_%d').

    Returns:
        List[datetime.date]: A sorted list of available sensor dates.
    """
    logger = logging.getLogger()
    logger.info(f"Fetching available sensor dates from the base URL: {base_url}")
    try:
        response = requests.get(base_url)
        response.raise_for_status()
    except requests.RequestException as e:
        logger.error(f"Failed to access base URL: {e}")
        return []
    
    soup = BeautifulSoup(response.content, 'html.parser')
    dates = []

    for link in soup.find_all('a'):
        href = link.get('href')
        if href:
            match = date_regex.search(href)
            if match:
                try:
                    date_str = match.group(1)
                    date = datetime.datetime.strptime(date_str, date_format).date()
                    dates.append(date)
                except ValueError:
                    logger.warning(f"Invalid date format in: {href}")

    logger.info(f"Found {len(dates)} available sensor dates.")
    return sorted(dates)


def get_sensor_data_date(landsat_date: datetime.date, base_url: str, date_regex: Pattern, date_format: str, N: float = 12) -> Optional[datetime.date]:
    """
    Find a matching sensor date within N days of the Landsat acquisition date.

    Args:
        landsat_date (datetime.date): The acquisition date of the Landsat image.
        base_url (str): URL to the sensor data directory on the server.
        date_regex (Pattern): Regex pattern to extract date from filenames.
        N (float): The number of days to look back and forward from the Landsat date.

    Returns:
        Optional[datetime.date]: The matching sensor date if found, else None.
    """
    logger = logging.getLogger()
    sensor_dates = get_available_sensor_dates(base_url, date_regex, date_format)
    logger.info(f"Available sensor dates: {sensor_dates}")
    for sensor_date in sensor_dates:
        if landsat_date <= sensor_date and landsat_date >= sensor_date - datetime.timedelta(days=N):
            logging.info(f"Matching sensor date found: {sensor_date}")
            return sensor_date
    logger.warning("No matching sensor date found.")
    return None


def sensor_data_transformation(save_data_loc: Union[str, Path], sensor_data_path: Union[str, Path])-> None:
    """
    Transform raw sensor data from a CSV file into a structured, daily-mean format.

    This function processes:
    - 'Sensor_Data.csv' (raw sensor measurements)
    using metadata provided in:
    - 'Sensor_Types.csv'
    - 'Sensor_Database.csv'

    These two metadata files are required and must be **manually created by the user**.

    Required files and formats:
    ---------------------------
    1. Sensor_Types.csv
       Maps each sensor ID to the type of measurement.

       Example:
           Sensors, Sensor
           101, Water Content
           102, Temperature

    2. Sensor_Database.csv
       Maps each sensor ID to its geographic location.

       Example:
           Sensors, lat, lon
           101, 23.7806, 90.4074
           102, 24.9045, 91.8611

    Args:
        save_data_loc (Union[str, Path]): Main project directory.
        sensor_data_path (Union[str, Path]): The folder (relative to save_data_loc) containing the CSV files.

    Returns:
        None
    """
    logger = logging.getLogger()
    logger.critical(f"Starting sensor data transformation for path: {sensor_data_path}")
    sensor_data_path = os.path.join(save_data_loc, sensor_data_path)
    csv_path, output_path = f'{sensor_data_path}/'+'Sensor_Data.csv',f'{sensor_data_path}/Transformed_Sensor_Data.csv'
    try:
        df = pd.read_csv(csv_path)
    except FileNotFoundError:
        logger.error(f"Sensor_Data.csv not found at: {csv_path}")
        return
    if 'Date' not in df.columns:
        logger.error("'Date' column is missing in Sensor_Data.csv.")
        return
    df['Date'] = pd.to_datetime(df['Date']).dt.date  # Convert to date immediately for daily grouping

    df = df.dropna().reset_index(drop=True)
    data_list = []
    try:
        measure_dict = pd.read_csv(os.path.join(sensor_data_path, 'Sensor_Types.csv')).set_index('Sensors').to_dict(orient='index')
        sensor_dict = pd.read_csv(os.path.join(sensor_data_path, 'Sensor_Database.csv')).set_index('Sensors').to_dict(orient='index')
    except FileNotFoundError as e:
        logger.error(f"Required metadata file missing: {e.filename}")
        return
    # Iterate through each column (excluding 'Date' and 'Line#')
    for column in df.columns.drop(['Line#', 'Date']):
        # Extract the sensor ID and measurement type
        sensor_id = column.split(':')[1].split(")")[0]
        measure_type = column.split(',')[1].strip() 
        # Create a temporary DataFrame with necessary columns
        temp_df = df[['Date', column]].copy()
        temp_df.columns = ['Date', 'Sensor Value']
        temp_df['Sensor ID'] = sensor_id
        temp_df['Measurement Unit'] = measure_type

        # Assign latitude and longitude from the dictionary
        if sensor_id in sensor_dict:
            temp_df['Latitude'] = sensor_dict[sensor_id]['lat']
            temp_df['Longitude'] = sensor_dict[sensor_id]['lon']
        else:
            temp_df['Latitude'] = np.nan  # Assign NaN if no location data is available
            temp_df['Longitude'] = np.nan
        if str(sensor_id) in measure_dict:
            temp_df['Measurement Type'] = list(measure_dict[sensor_id].values())[0]
        else:
            temp_df['Measurement Type'] = 'Unknown'
        # Append to the list
        data_list.append(temp_df)
    if not data_list:
        logger.error("No valid sensor data found for processing.")
        return
    # Concatenate all data into a single DataFrame
    all_sensors_df = pd.concat(data_list)
    # all_sensors_df = all_sensors_df.dropna().reset_index(drop=True)
    # Group by 'Sensor ID', 'Measurement Type', 'Latitude', 'Longitude', and 'Date', then compute daily mean
    daily_means = all_sensors_df.groupby(['Sensor ID', 'Measurement Type', 'Measurement Unit','Latitude', 'Longitude', 'Date']).agg({'Sensor Value': 'mean'})
    daily_means = daily_means.dropna()
    # # Reset index to make the DataFrame more accessible if necessary
    daily_means.reset_index(inplace=True)
    logger.info(f"Saving transformed data to: {output_path}")
    daily_means.to_csv(output_path,index = False)
    logger.info("Sensor data transformation complete.")


def inverse_distance_weighting(
    lats: np.ndarray, lons: np.ndarray, values: np.ndarray, 
    grid_lat: np.ndarray, grid_lon: np.ndarray, power: float = 2.0
) -> np.ndarray:
    """
    Apply IDW interpolation to create a raster over a grid.
    """
    interpolated = np.zeros_like(grid_lat)
    for i in range(grid_lat.shape[0]):
        for j in range(grid_lat.shape[1]):
            dists = np.sqrt((lats - grid_lat[i, j])**2 + (lons - grid_lon[i, j])**2)
            weights = 1 / (dists**power + 1e-12)
            interpolated[i, j] = np.sum(weights * values) / np.sum(weights)
    return interpolated


def interpolate_grid(
    method: str,
    lats: np.ndarray,
    lons: np.ndarray,
    values: np.ndarray,
    lat_grid: np.ndarray,
    lon_grid: np.ndarray,
    power: float = 2.0
) -> np.ndarray:
    """Performs spatial interpolation using various methods."""
    coords = np.column_stack((lons, lats))
    points = np.column_stack((lon_grid.ravel(), lat_grid.ravel()))

    if method == 'idw':
        return inverse_distance_weighting(lats, lons, values, lat_grid, lon_grid, power)

    elif method == 'average':
        return np.full_like(lat_grid, np.mean(values))

    elif method == 'nearest_neighbor':
        tree = cKDTree(coords)
        dist, idx = tree.query(points)
        return values[idx].reshape(lat_grid.shape)

    elif method in ['linear', 'cubic', 'nearest']:
        grid_z = griddata(coords, values, points, method=method, fill_value=np.nan)
        return grid_z.reshape(lat_grid.shape)

    elif method == 'best':
        # Custom logic: use closest station for now
        tree = cKDTree(coords)
        dist, idx = tree.query(points, k=1)
        return values[idx].reshape(lat_grid.shape)

    else:
        raise NotImplementedError(f"Interpolation method '{method}' is not supported.")

# def numpy_to_ee_image(
#     array: np.ndarray,
#     lon_grid: np.ndarray,
#     lat_grid: np.ndarray,
#     crs: str = 'EPSG:4326',
#     scale: float = 25000
# ) -> ee.Image:
#     """
#     Convert a NumPy array and corresponding lat/lon grids to an ee.Image.

#     Args:
#         array (np.ndarray): 2D array with interpolated values.
#         lon_grid (np.ndarray): 2D array of longitudes (same shape as array).
#         lat_grid (np.ndarray): 2D array of latitudes (same shape as array).
#         crs (str): Coordinate reference system (default EPSG:4326).
#         scale (float): Scale of the image in meters (default 25000m).

#     Returns:
#         ee.Image: Earth Engine image with the interpolated values.
#     """
#     if array.shape != lon_grid.shape or array.shape != lat_grid.shape:
#         raise ValueError("Shapes of array, lon_grid, and lat_grid must match.")

#     # Flatten all arrays
#     flat_values = array.flatten().tolist()
#     flat_coords = list(zip(lon_grid.flatten(), lat_grid.flatten()))

#     # Create list of features
#     features = [
#         ee.Feature(ee.Geometry.Point([lon, lat]), {"value": val})
#         for (lon, lat), val in zip(flat_coords, flat_values)
#         if not np.isnan(val)
#     ]

#     fc = ee.FeatureCollection(features)
#     image = fc.reduceToImage(properties=["value"], reducer=ee.Reducer.first())

#     return image


def bias_correction_with_station_data(image: Image, 
    roi: Geometry,
    start_date: Union[str, pd.Timestamp], 
    station_csv_path: Union[str, Path], 
    save_data_loc: Union[str, Path], 
    variable: str = "Temperature", 
    scale: int = 25000, 
    interpolation_method: Literal['idw', 'average', 'nearest_neighbor', 'linear', 'cubic', 'best'] = 'idw', 
    power: float = 2.0, 
    crs: str = 'EPSG:4326'
    ) -> Image:
    """
    Applies bias correction to GFS data (e.g., temperature, wind) using in-situ sensor values.

    Parameters:
        image (ee.Image): The GFS image containing the variable of interest.
        roi (ee.Geometry): The region of interest over which the mean is calculated.
        start_date (str or pd.Timestamp): The date for which correction is applied.
        sensor_data_path (Path): Path to the folder containing sensor data CSV.
        save_data_loc (Path): Main project working directory.
        variable (str): The variable to correct (e.g., "Temperature", "Wind", and "Specific Humidity").
        scale (int): Scale (in meters) for the satellite image resolution.
        interpolation_method (str): Method to interpolate station data.
        power (float): Power parameter for IDW.

    Returns:
        ee.Image: The bias-corrected image.
    """
    logger = logging.getLogger()
    station_csv_path = os.path.join(save_data_loc, station_csv_path)
    logger.info(f'Starting bias correction for variable: {variable}')
    band_map = {
        "Temperature": "temperature_2m_above_ground",
        "Wind Speed": "wind",
        "Specific Humidity": "specific_humidity_2m_above_ground",
        "Pressure": 'pressure'
    }
    band_name = band_map.get(variable)
    if band_name is None:
        logger.error(f"Unsupported variable: {variable}. Supported variables are: {list(band_map.keys())}")
        raise ValueError(f"Unsupported variable: {variable}")
    
    # df_path = os.path.join(sensor_data_path,station_csv_path)
    logger.critical(f'Using station data from: {station_csv_path}')
    try:
        df = pd.read_csv(station_csv_path, parse_dates=['Date'])
    except FileNotFoundError:
        logger.critical(f"Station data file not found at: {station_csv_path}. Please ensure the file exists and the path is correct.")
        raise FileNotFoundError(f"Station data file not found at: {station_csv_path}. Please ensure the file exists and the path is correct.")
    df_on_date = df[df['Date'] == pd.to_datetime(start_date)]
    if df_on_date.empty:
        raise ValueError(f"No station data found for the specified date: {start_date}")

    if variable not in df_on_date.columns or not {'lat', 'lon'}.issubset(df_on_date.columns):
        raise ValueError(f"CSV must include 'lat', 'lon', and '{variable}' columns.  Found columns: {df_on_date.columns.tolist()}")
    

    lats = df_on_date['lat'].to_numpy()
    lons = df_on_date['lon'].to_numpy()
    values = df_on_date[variable].to_numpy()
    logger.info(f"Generating target lat-lon mesh using region at scale {scale}")
    latlon_img = ee.Image.pixelLonLat().reproject(crs=crs, scale=scale)
    latlon_sample = latlon_img.sample(region=roi, scale=scale, geometries=True)
    coords = latlon_sample.aggregate_array('coordinates').getInfo()
    lon_grid = np.array([c[0] for c in coords])
    lat_grid = np.array([c[1] for c in coords])

    lon_mesh, lat_mesh = np.meshgrid(
        np.linspace(min(lon_grid), max(lon_grid), int(np.sqrt(len(coords)))),
        np.linspace(min(lat_grid), max(lat_grid), int(np.sqrt(len(coords))))
    )
    logger.critical(f"Interpolating station data using method: {interpolation_method}")
    interpolated = interpolate_grid(
        method=interpolation_method,
        lats=lats, lons=lons, values=values,
        lat_grid=lat_mesh, lon_grid=lon_mesh,
        power=power
    )
    pixel_distance = geodesic(
        (lat_grid[0, 0], lon_grid[0, 0]),
        (lat_grid[0, 1], lon_grid[0, 1])
    ).meters

    station_interp_image = numpy_to_ee_image(interpolated, lon_mesh, lat_mesh, scale = pixel_distance, crs=crs)

    satellite_band = image.select(variable)
    logger.critical("Calculating bias and correcting image")
    bias = satellite_band.subtract(station_interp_image)
    mean_bias = bias.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=roi,
        scale=scale,
        maxPixels=1e9
    ).get(variable)
    logger.critical(f"Mean bias for {variable} on {start_date}: {mean_bias.getInfo()}")
    corrected = satellite_band.subtract(bias).rename(f"{variable}_corrected")

    return corrected.set({'bias_corrected': True, 'bias_date': start_date})


# def numpy_to_ee_image(array: np.ndarray, lon_grid: np.ndarray, lat_grid: np.ndarray, scale: float, crs: str, roi: ee.Geometry) -> ee.Image:
#     """Converts numpy array and lat/lon grid to Earth Engine image."""
    

#     with tempfile.NamedTemporaryFile(suffix=".tif") as tmp:
#         pixel_x = lon_grid[0, 1] - lon_grid[0, 0] if lon_grid.shape[1] > 1 else scale / 111320
#         pixel_y = lat_grid[0, 0] - lat_grid[1, 0] if lat_grid.shape[0] > 1 else scale / 110540
#         transform = from_origin(
#             lon_grid[0, 0], lat_grid[0, 0], pixel_x, pixel_y
#         )
#         with rio.open(
#             tmp.name,
#             "w",
#             driver="GTiff",
#             height=array.shape[0],
#             width=array.shape[1],
#             count=1,
#             dtype=array.dtype,
#             crs=crs,
#             transform=transform,
#         ) as dst:
#             dst.write(array, 1)
#         return ee.Image.loadGeoTIFF(tmp.name).clip(roi)

def numpy_to_ee_image(array: np.ndarray, lon_grid: np.ndarray, lat_grid: np.ndarray, scale: float, crs: str, roi: ee.Geometry) -> ee.Image:
    """Converts numpy array and lat/lon grid to Earth Engine image."""

    # Create temp file WITHOUT automatic deletion (Windows workaround)
    tmp = tempfile.NamedTemporaryFile(suffix=".tif", delete=False)
    tmp.close()  # Close it so rasterio can access it

    try:
        pixel_x = lon_grid[0, 1] - lon_grid[0, 0] if lon_grid.shape[1] > 1 else scale / 111320
        pixel_y = lat_grid[0, 0] - lat_grid[1, 0] if lat_grid.shape[0] > 1 else scale / 110540
        transform = from_origin(
            lon_grid[0, 0], lat_grid[0, 0], pixel_x, pixel_y
        )

        with rio.open(
            tmp.name,
            "w",
            driver="GTiff",
            height=array.shape[0],
            width=array.shape[1],
            count=1,
            dtype=array.dtype,
            crs=crs,
            transform=transform,
        ) as dst:
            dst.write(array, 1)

        # Load and return the EE image
        return ee.Image(tmp.name).clip(roi)

    finally:
        # Clean up temp file
        os.remove(tmp.name)
    


def interpolate_bias_surface(
    image : ee.Image,    
    roi: ee.Geometry,
    target_date: Union[str, pd.Timestamp],
    station_csv_path: Path,
    save_data_loc: Union[str, Path],
    variable: str = "Temperature",
    scale: float = 25000,
    interpolation_method: Literal['idw', 'average', 'nearest_neighbor', 'linear', 'cubic', 'best'] = 'idw',
    crs : str = "EPSG:4326"
) -> ee.Image:
    """
    Applies bias correction to a GEE image using in-situ station data interpolation.

    Args:
        image (ee.Image): Input image with the specified variable band.
        roi (ee.Geometry): The region of interest to interpolate over.
        target_date (str or pd.Timestamp): Date to filter station data.
        station_csv_path (Path): Path to CSV with station data (lat, lon, value, date).
        save_data_loc (Path): Main project working directory.
        interpolation_method (str): 'idw', 'average', 'nearest_neighbor', 'kriging'.
        scale (float): Pixel size in meters.

    Returns:
        ee.Image: Bias-corrected image with properties set.
    """
    logger = logging.getLogger()
    logger.info(f"Starting bias correction for {variable} on {target_date}")
    band_map = {
        "Temperature": "temperature_2m_above_ground",
        "Wind Speed": "wind",
        "Specific Humidity": "specific_humidity_2m_above_ground",
        "Pressure": 'pressure'
    }
    band_name = band_map.get(variable)
    if band_name is None:
        logger.error(f"Unsupported variable: {variable}. Supported variables are: {list(band_map.keys())}")
        raise ValueError(f"Unsupported variable: {variable}")
    station_csv_path = os.path.join(save_data_loc, station_csv_path)
    logger.info(f"Reading station data from {station_csv_path}")
    df = pd.read_csv(station_csv_path)
    df['Date'] = pd.to_datetime(df['Date']).dt.date
    df = df[df['Date'] == target_date]

    if df.empty:
        logger.warning(f"No station data for target date: {target_date}")
        raise ValueError("No station data for given date.")
    for col in ['lat', 'lon', variable]:
        if col not in df.columns:
            raise ValueError(f"{col} column missing in CSV. Available columns: {df.columns.tolist()}")
    
    # lats = df['lat'].values
    # lons = df['lon'].values
    # vals = df[variable].values

    df = df[['lat', 'lon', variable]].dropna()
    gdf = gpd.GeoDataFrame(df, geometry=gpd.points_from_xy(df.lon, df.lat), crs = crs)

    # Convert ROI bounds to grid
    roi_coords = roi.bounds().getInfo()['coordinates'][0]
    lons = [c[0] for c in roi_coords]
    lats = [c[1] for c in roi_coords]
    min_lon, max_lon = min(lons), max(lons)
    min_lat, max_lat = min(lats), max(lats)

    # Compute number of pixels
    lon_res = scale / 111320  # rough meters per degree longitude
    lat_res = scale / 110540  # rough meters per degree latitude
    n_cols = max(2, int((max_lon - min_lon) / lon_res))
    n_rows = max(2, int((max_lat - min_lat) / lat_res))

    if n_cols < 2 or n_rows < 2:
        logger.warning("ROI is smaller than scale; fallback to nearest station interpolation.")
        interpolation_method = "nearest_neighbor"

    grid_lon = np.linspace(min_lon, max_lon, n_cols)
    grid_lat = np.linspace(max_lat, min_lat, n_rows)
    lon_grid, lat_grid = np.meshgrid(grid_lon, grid_lat)

    if interpolation_method == "idw":
        tree = cKDTree(gdf[["lon", "lat"]].values)
        dists, idxs = tree.query(np.c_[lon_grid.ravel(), lat_grid.ravel()], k=4)
        weights = 1 / (dists + 1e-6)
        values_arr = gdf[variable].values  # shape (n_stations,)
        selected_vals = values_arr[idxs]    # shape (n_points, k)
        interpolated = np.sum(weights * selected_vals, axis=1) / np.sum(weights, axis=1)
        interpolated = interpolated.reshape(lat_grid.shape)

    elif interpolation_method == "nearest_neighbor":
        tree = cKDTree(gdf[["lon", "lat"]].values)
        dists, idxs = tree.query(np.c_[lon_grid.ravel(), lat_grid.ravel()], k=1)
        interpolated = gdf.iloc[idxs][variable].values.reshape(lat_grid.shape)

    elif interpolation_method == "average":
        interpolated = np.full_like(lon_grid, gdf[{variable}].mean())

    elif interpolation_method == "kriging":
        OK = OrdinaryKriging(
            gdf["lon"].values,
            gdf["lat"].values,
            gdf[{variable}].values,
            variogram_model="linear",
            verbose=False,
            enable_plotting=False,
        )
        interpolated, _ = OK.execute("grid", grid_lon, grid_lat)
        interpolated = np.flipud(interpolated)

    else:
        raise ValueError(f"Unsupported interpolation method: {interpolation_method}")
        # Compute pixel size robustly
    if lon_grid.shape[1] > 1 and lat_grid.shape[0] > 1:
        pixel_x = lon_grid[0, 1] - lon_grid[0, 0]
        pixel_y = lat_grid[0, 0] - lat_grid[1, 0]
    else:
        logger.warning("Grid is too small; using estimated pixel sizes")
        pixel_x = lon_res
        pixel_y = lat_res

    logger.info("Converting interpolated surface to EE image")
    pixel_distance = geodesic(
        (lat_grid[0, 0], lon_grid[0, 0]),
        (lat_grid[0, 0] - pixel_y, lon_grid[0, 0] + pixel_x),
    ).meters
    # ee_img = numpy_to_ee_image(interpolated.astype(np.float32), lon_grid, lat_grid, scale = pixel_distance, roi = roi, crs = crs)
    logger.info("Computing bias-corrected image...")
    satellite_band = image.select(band_name)
    mean_interp_value = float(np.nanmean(interpolated))
    bias = satellite_band.subtract(mean_interp_value)
    mean_bias = bias.reduceRegion(
        reducer=ee.Reducer.mean(),
        geometry=roi,
        scale=scale,
        maxPixels=1e9
    ).get(band_name)
    logger.info(f"Mean bias for {variable}: {mean_bias.getInfo()}")
    corrected = satellite_band.subtract(bias)
    corrected = corrected.copyProperties(
        source=image,
        exclude=[]  # leave empty list if you want all properties
    )
    corrected = ee.Image(corrected).set(
        {
            'bias_corrected': True,
            'bias_date': str(target_date),
        }
    )
    return corrected