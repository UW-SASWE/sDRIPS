import os
import ee
import geemap
import time
import logging

from sdrips.utils.utils import load_yaml_config
from sdrips.utils.ee_initialize import initialize_earth_engine


def initialize_earth_engine():
    """
    Authenticate and initialize the Earth Engine API.

    If the user is not authenticated, this function will prompt for authentication.
    Once authenticated, it initializes the Earth Engine library.

    Raises:
        Exception: If initialization fails after authentication.
    """
    logger = logging.getLogger()

    try:
        ee.Initialize()
        logger.info("Earth Engine initialized successfully.")
    except ee.ee_exception.EEException as init_error:
        logger.critcal("Earth Engine not initialized. Attempting authentication...")
        try:
            ee.Authenticate()
            ee.Initialize()
            logger.info("Earth Engine authenticated and initialized successfully.")
        except ee.ee_exception.EEException as auth_error:
            logger.error(f"Earth Engine authentication or initialization failed: {auth_error}")
            raise
    except Exception as e:
        logger.error(f"An unexpected error occurred during Earth Engine initialization: {e}")
        raise Exception("Failed to initialize Earth Engine after authentication.") from e



def create_ee_asset_folder(asset_folder_id: str) -> None:
    """
    Create an Earth Engine asset folder if it does not already exist.

    Args:
        asset_folder_id (str): Full asset path for the folder
            (e.g., 'users/username/foldername').

    Raises:
        ee.EEException: If an error occurs during the creation.
    """
    try:
        logger = logging.getLogger()
        ee.data.getAsset(asset_folder_id)
        logger.critical(f"Asset folder already exists: {asset_folder_id}")
    except ee.EEException:
        try:
            parent_path = "/".join(asset_folder_id.split("/")[:-1])
            folder_name = asset_folder_id.split("/")[-1]
            logger.critical(f"Creating asset folder: {asset_folder_id}")
            ee.data.createAsset(
                {'type': 'Folder'}, parent_path, folder_name
            )
            logger.critical(f"Created folder: {asset_folder_id}")
        except ee.EEException as e:
            logger.error(f"Failed to create asset folder: {asset_folder_id}")
            raise e


def upload_shapefile_to_ee(shp_path: str, asset_folder: str = "sdrips/") -> ee.FeatureCollection:
    """
    Uploads a local shapefile to the specified Earth Engine asset folder and returns it as an ee.FeatureCollection.
    
    Args:
        shp_path (str): Path to the local shapefile.
        asset_folder (str): Folder in the user's EE assets where the shapefile will be uploaded.
    
    Returns:
        ee.FeatureCollection: The uploaded FeatureCollection reference.
    """
    initialize_earth_engine()
    logger = logging.getLogger()
    shp_name = os.path.splitext(os.path.basename(shp_path))[0]
    # logger.critical(f'shp_name: {shp_name}')
    root = ee.data.getAssetRoots()[0]
    username= root["id"].replace("projects/earthengine-legacy/assets/", "")
    # logger.critical(f'username: {username}')
    asset_folder_id = f"{username}/{asset_folder}"
    # logger.critical(f'asset_folder_id: {asset_folder_id}')
    asset_id = f"{asset_folder_id}{shp_name}"
    # logger.critical(f'asset_id: {asset_id}')
    create_ee_asset_folder(asset_folder_id)
    try:
        ee.data.getAsset(asset_id)
        logger.critical(f"Asset already exists: {asset_id}")
    except ee.EEException:
        logger.critical(f"Uploading new asset to EE: {asset_id}")
        # logger.critical(f"shp_path: {shp_path}")
        ee_fc = geemap.shp_to_ee(shp_path)
        task = ee.batch.Export.table.toAsset(
            collection=ee_fc,
            description=f"upload_{shp_name}",
            assetId=asset_id
        )
        task.start()

        logger.critical("Waiting for export to complete...")
        max_wait_time = 60 * 10  # 10 minutes in seconds
        poll_interval = 10  # seconds
        elapsed_time = 0

        while task.active():
            if elapsed_time >= max_wait_time:
                logger.error("Export task timed out after 30 minutes.")
                raise TimeoutError("Export task did not finish within 30 minutes.")
            logger.critical("Still exporting...")
            time.sleep(poll_interval)
            elapsed_time += poll_interval
        if task.status()['state'] == 'COMPLETED':
            logger.critical(f"Export completed successfully: {asset_id}")
        else:
            error_message = task.status().get("error_message", "Unknown error")
            logger.error(f"Export failed: {error_message}")
            raise RuntimeError(f"Export task failed: {error_message}")

    return ee.FeatureCollection(asset_id)
