import logging
import re
import ee
from ee.image import Image
from ee.geometry import Geometry
import os, datetime, requests, zipfile, time, math, glob
import pandas as pd
import numpy as np
from typing import List, Optional, Pattern, Union
from pathlib import Path
import urllib.request
from bs4 import BeautifulSoup
from ruamel.yaml import YAML
from filelock import FileLock
import json
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


def bias_correction(image: Image, roi: Geometry, start_date: Union[str, pd.Timestamp], sensor_data_path: Union[str, Path], save_data_loc: Union[str, Path], variable: str = "Temperature"):
    """
    Applies bias correction to GFS data (e.g., temperature, wind) using in-situ sensor values.

    Parameters:
        image (ee.Image): The GFS image containing the variable of interest.
        roi (ee.Geometry): The region of interest over which the mean is calculated.
        start_date (str or pd.Timestamp): The date for which correction is applied.
        sensor_data_path (Path): Path to the folder containing sensor data CSV.
        save_data_loc (Path): Main project working directory.
        variable (str): The variable to correct (e.g., "Temperature", "Wind", and "Specific Humidity").

    Returns:
        ee.Image: The bias-corrected image.
    """
    logger = logging.getLogger()
    sensor_data_path = os.path.join(save_data_loc, sensor_data_path)
    logger.info(f'Starting bias correction for variable: {variable}')
    band_map = {
        "Temperature": "temperature_2m_above_ground",
        "Wind": "wind",
        "Specific Humidity": "specific_humidity_2m_above_ground",
    }
    band_name = band_map.get(variable)
    if band_name is None:
        logger.error(f"Unsupported variable: {variable}. Supported variables are: {list(band_map.keys())}")
        raise ValueError(f"Unsupported variable: {variable}")
    
    df_path = f'{sensor_data_path}/Transformed_Sensor_Data.csv'
    logger.info(f'Using sensor data path: {df_path}')
    df = pd.read_csv(df_path)
    df['Date'] = pd.to_datetime(df['Date']).dt.date
    logger.info(f"Bias correction for {variable} on {start_date}")
    start_date = pd.to_datetime(start_date).date() 
    specific_date_df = df[df['Date'] ==  start_date]
    df = specific_date_df[specific_date_df['Measurement Type'] == variable]
    if df.empty:
        logger.critical(f"No sensor data for {start_date} and variable '{variable}'. Skipping correction.")
    def extract_value(row: pd.Series) -> float:
        point = ee.Geometry.Point([row['Longitude'], row['Latitude']])
        selscale = 30
        value = image.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=point,
            scale = selscale
        ).get(band_name)
        try:
            return value.getInfo() if value else np.nan
        except Exception as e:
            logger.warning(f"Error getting model value for point {point}: {e}")
            return np.nan

    df['Model Value'] = df.apply(extract_value, axis=1)
    df.dropna(subset=['Model Value'], inplace=True)

    if df.empty:
        logger.critical(f"No valid satellite data for {start_date}. Skipping correction.")

    mean_sensor = df['Sensor Value'].mean()
    mean_model = df['Model Value'].mean()
    bias = mean_sensor - mean_model
    logger.info(f"Mean sensor value: {mean_sensor:.2f}, Mean model value: {mean_model:.2f}, Bias: {bias:.2f}")
    corrected_image = image.add(bias)
    logger.info(f'Bias correction for variable {variable} completed successfully.')

    return corrected_image

