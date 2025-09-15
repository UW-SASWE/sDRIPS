import os
import subprocess
import ee
import geemap
import time
from datetime import datetime
import logging
from google.cloud.resourcemanager import ProjectsClient

# from sdrips.utils.utils import load_yaml_config
# from sdrips.utils.ee_initialize import initialize_earth_engine



def initialize_earth_engine(service_account: str = None, key_file: str = None):
    """
    Authenticate and initialize the Earth Engine API.

    If the user is not authenticated, this function will prompt for authentication.
    Once authenticated, it initializes the Earth Engine library.

    Args:
        service_account (str, optional): Service account email.
        key_file (str, optional): Path to the service account key JSON file.

    Raises:
        Exception: If initialization fails after authentication.
    """
    logger = logging.getLogger()

    if service_account and key_file:
        logger.info(f"Initializing Earth Engine with service account: {service_account}")
        try:
            credentials = ee.ServiceAccountCredentials(service_account, key_file)
            ee.Initialize(credentials)
            logger.info(f"Earth Engine initialized with service account: {service_account}")
            return
        except Exception as sa_error:
            logger.warning(f"Service account initialization failed: {sa_error}")

    try:
        client = ProjectsClient()
        for project in client.search_projects():
            # Check if the project has the 'earth-engine' label
            if "earth-engine" in project.labels and "default project" in project.display_name.lower():
                default_project = project.project_id
                break
        logger.info(f"Default Earth Engine project identified: {default_project}")
        # default_project = ee.data.getAssetRoots()[0]["id"].split("/")[1]
        ee.Initialize(project = default_project)
        logger.info("Earth Engine initialized successfully.")
    except ee.EEException:
        logger.critical("Earth Engine not initialized. Attempting authentication...")
        try:
            ee.Authenticate()
            client = ProjectsClient()
            for project in client.search_projects():
                # Check if the project has the 'earth-engine' label
                if "earth-engine" in project.labels and "default project" in project.display_name.lower():
                    default_project = project.project_id
                    break
            logger.info(f"Default Earth Engine project identified: {default_project}")
            # default_project = ee.data.getAssetRoots()[0]["id"].split("/")[1]
            ee.Initialize(project = default_project)
            logger.info("Earth Engine authenticated and initialized successfully.")
        except ee.ee_exception.EEException as auth_error:
            try:
                ee.Authenticate()
                ee.Initialize()
            except Exception as e:
                logger.error(f"Earth Engine authentication or initialization failed: {auth_error}")
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


# def upload_shapefile_to_ee(shp_path: str, asset_folder: str = "sdrips/") -> ee.FeatureCollection:
#     """
#     Uploads a local shapefile to the specified Earth Engine asset folder and returns it as an ee.FeatureCollection.
    
#     Args:
#         shp_path (str): Path to the local shapefile.
#         asset_folder (str): Folder in the user's EE assets where the shapefile will be uploaded.
    
#     Returns:
#         ee.FeatureCollection: The uploaded FeatureCollection reference.
#     """
#     initialize_earth_engine()
#     logger = logging.getLogger()

#     if not os.path.exists(shp_path):
#         raise FileNotFoundError(f"Shapefile not found: {shp_path}")
    
#     root = ee.data.getAssetRoots()[0]["id"]
#     logger.info(f"EE root: {root}")

#     if root.startswith("projects/"):
#         idx = root.find("/assets")
#         base = root[:idx + len("/assets")]  # ensures base ends with /assets
#     else:
#         base = root

#     asset_folder_id = f"{base.rstrip('/')}/{asset_folder.strip('/')}"

#     parts = asset_folder_id.split("/")
#     if parts[0] == "projects":
#         start_idx = 3  # ['projects', '<id>', 'assets']
#     elif parts[0] == "users":
#         start_idx = 2
#     else:
#         raise ValueError(f"Unknown root format: {parts[0]}")

#     for i in range(start_idx, len(parts)+1):
#         prefix = "/".join(parts[:i])
#         try:
#             ee.data.getAsset(prefix)
#             logger.info(f"Folder exists: {prefix}")
#         except ee.EEException:
#             logger.info(f"Creating folder: {prefix}")
#             ee.data.createFolder(prefix)
#             logger.info(f"Created folder: {prefix}")

#     shp_name = os.path.splitext(os.path.basename(shp_path))[0]
#     asset_id = f"{asset_folder_id}/{shp_name}"
#     logger.info(f"Final asset ID: {asset_id}")

#     try:
#         ee.data.getAsset(asset_id)
#         logger.info(f"Asset already exists: {asset_id}")
#     except ee.EEException:
#         logger.info(f"Uploading new asset to EE: {asset_id}")
#         ee_fc = geemap.shp_to_ee(shp_path)
#         task = ee.batch.Export.table.toAsset(
#             collection=ee_fc,
#             description=f"upload_{shp_name}",
#             assetId=asset_id
#         )
#         task.start()
#         logger.info("Started export task. Waiting for completion...")

#         max_wait_time = 60 * 10  # 10 minutes
#         poll_interval = 10  # seconds
#         elapsed = 0
#         while task.active():
#             if elapsed >= max_wait_time:
#                 logger.error("Export task timed out after 10 minutes.")
#                 raise TimeoutError("Export task did not finish within 10 minutes.")
#             logger.info("Still exporting...")
#             time.sleep(poll_interval)
#             elapsed += poll_interval

#         status = task.status()
#         if status.get("state") == "COMPLETED":
#             logger.info("Export completed successfully.")
#         else:
#             error_message = status.get("error_message", "Unknown error")
#             logger.error(f"Export failed: {error_message}")
#             raise RuntimeError(f"Export task failed: {error_message}")

#     return ee.FeatureCollection(asset_id)


def upload_shapefile_to_ee(shp_path: str, asset_folder: str = "sdrips_folder", dry_run: bool = False, service_account: str = None, key_file: str = None) -> ee.FeatureCollection:
    """
    Upload a local shapefile to the specified EE asset folder, creating folders if needed.

    Args:
        shp_path (str): Path to local shapefile.
        asset_folder (str): Folder name for EE assets.
        dry_run (bool): If True, just compute and return asset_id as ee.FeatureCollection without exporting.

    Returns:
        ee.FeatureCollection: Reference to the uploaded EE asset.
    """
    assert os.path.exists(shp_path), f"Shapefile not found: {shp_path}"
    ensure_ee_initialized(service_account=service_account, key_file=key_file)
    logger = logging.getLogger(__name__)
    
    shp_name = os.path.splitext(os.path.basename(shp_path))[0]

    roots = ee.data.getAssetRoots()
    if roots == []:
        client = ProjectsClient()
        for project in client.search_projects():
            if "earth-engine" in project.labels and "default project" in project.display_name.lower():
                default_project = project.project_id
                break  # once condition is met come out of the for loop

        project_id = default_project
        shp_name = os.path.splitext(os.path.basename(shp_path))[0]

        asset_folder = f"projects/{project_id}/assets/{asset_folder}"
        result = subprocess.run(
            ["earthengine", "ls", asset_folder],
            capture_output=True,
            text=True
        )
        logger.debug(f'result: {result.stdout}')
        if "not found" in result.stdout.lower():
            logger.info(f"Folder does not exist. Creating: {asset_folder}")
            subprocess.run(["earthengine", "create", "folder", asset_folder])
        else:
            logger.info(f"Folder exists: {asset_folder}")
        base_folder = asset_folder

    else:
        root = roots[0]["id"]
        if asset_folder not in roots:
            base_folder = f"{root}/{asset_folder}"
        else:
            base_folder = asset_folder
        try:
            ee.data.getAsset(base_folder)
            logger.critical(f"Folder already exists: {base_folder}")
        except ee.EEException:
            logger.critical(f"Creating folder: {base_folder}")
            ee.data.createFolder(base_folder)

    asset_id = f"{base_folder}/{shp_name}"

    try:
        ee.data.getAsset(asset_id)
        asset_id = f"{base_folder}/{shp_name}"
        logger.critical(f"Asset exists. Returning: {asset_id}")
        return ee.FeatureCollection(asset_id)
    except ee.EEException:
        logger.critical(f"Uploading new asset: {asset_id}")

    if dry_run:
        logger.info(f"Dry run - asset ID would be: {asset_id}")
        return ee.FeatureCollection(asset_id)

    # Convert shapefile to EE FeatureCollection
    ee_fc = geemap.shp_to_ee(shp_path)
    task = ee.batch.Export.table.toAsset(
        collection=ee_fc,
        description=f"upload_{shp_name}",
        assetId=asset_id
    )
    task.start()

    logger.critical("Waiting for export to complete...")
    max_wait_time = 10 * 60  # 10 minutes
    poll_interval = 10
    elapsed = 0
    while task.active():
        if elapsed >= max_wait_time:
            logger.error("Export timed out after 10 minutes")
            raise TimeoutError("Export did not finish within 10 minutes")
        logger.critical("Still exporting...")
        time.sleep(poll_interval)
        elapsed += poll_interval

    status = task.status()
    if status["state"] == "COMPLETED":
        logger.critical(f"Export completed: {asset_id}")
    else:
        error_message = status.get("error_message", "Unknown error")
        logger.error(f"Export failed: {error_message}")
        raise RuntimeError(f"Export task failed: {error_message}")

    return ee.FeatureCollection(asset_id)


def ensure_ee_initialized(service_account: str = None, key_file: str = None):
    try:
        ee.Number(1).getInfo()
    except Exception:
        initialize_earth_engine(service_account=service_account, key_file=key_file)
