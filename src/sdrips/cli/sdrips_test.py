# import argparse
# import requests
# from pathlib import Path
# import sys
# from typing import Dict
# from ruamel.yaml import YAML
# import multiprocessing
# from pathlib import Path
# import pandas as pd
# import numpy as np
# import rasterio

# from utils.utils import get_gdrive_url
# from sdrips.run_sdrips import run_sdrips

# # Google Drive File IDs 
# CONFIG_SOURCES = {
#     "ca_config.yaml": "1DI_e47mdscV_VZBgSFiXGvShwm07g_zk", 
#     "rupasi.shp": "1DI_e47mdscV_VZBgSFiXGvShwm07g_zk",
#     "folder": "1DI_e47mdscV_VZBgSFiXGvShwm07g_zk" # should not be downloaded as zip in download_files function. Also should be downloaded in the folder of expected_output
# }


# # Download config files from Google Drive
# def download_files():
#     for file_name, file_id in CONFIG_SOURCES.items():
#         gdrive_url = get_gdrive_url(file_id)
#         response = requests.get(gdrive_url)
#         if response.status_code == 200:
#             if file_name.endswith(".yaml"):
#                 # save in test folder in config_files, create test folder directory if not exists
#                 (Path("config_files/test")).mkdir(parents=True, exist_ok=True)  
#                 # use joinpath with test and config files
#                 with open(Path("config_files").joinpath("test", file_name), 'wb') as f:
#                     f.write(response.content)
#                 print(f"Downloaded {file_name}")
#             if file_name == "folder":
#                 # save in the expected_output folder, create directory if not exists. this expected_output folder should be in test_rupasi
#                 (Path("expected_output/test_rupasi")).mkdir(parents=True, exist_ok=True)
#                 with open(Path("expected_output/test_rupasi").joinpath(file_name), 'wb') as f:
#                     f.write(response.content)
#             else:
#                 # save in Shapefiles/Rupasi, create Rupasi directory if not exists
#                 (Path("Shapefiles/Rupasi")).mkdir(parents=True, exist_ok=True)  
#                 with open(Path("Shapefiles").joinpath("Rupasi", file_name), 'wb') as f:
#                     f.write(response.content)
#                 print(f"Downloaded {file_name}")
#         else:
#             print(f"Failed to download {file_name}")

# def update_config_file_test(config_path, save_loc=None):
#     """
#     Updates multiprocessing cores and save_data_loc in the YAML config.

#     Args:
#         config_path (str or Path): Path to the sdrips_config.yaml file.
#         save_loc (str or Path, optional): Location to save data. Defaults to current directory.

#     """
#     config_path = Path(config_path)
#     save_loc = save_loc or str(Path.cwd())  

#     yaml = YAML()
#     yaml.preserve_quotes = True  

#     with open(config_path, 'r') as f:
#         config = yaml.load(f)

#     # Update multiprocessing cores
#     max_cores = max(1, multiprocessing.cpu_count() - 1)
#     if "Multiprocessing" in config:
#         config["Multiprocessing"]["cores"] = max_cores
#     else:
#         raise ValueError("Multiprocessing section not found in config.")

#     # Update save_data_loc, i want to add test directory in the save_data_loc,pls do it for me
#     if "Save_Data_Location" in config:
#         config["Save_Data_Location"]["save_data_loc"] = str(Path(save_loc).joinpath("test_rupasi"))
#     else:
#         config["Save_Data_Location"] = {"save_data_loc": str(Path(save_loc).joinpath("test_rupasi"))}

#     # Update Shapefiles path
#     if "Irrigation_cmd_area_shapefile" in config:
#         config["Irrigation_cmd_area_shapefile"]["path"] = str(Path(save_loc).joinpath("Shapefiles", "Rupasi", "Rupasi.shp"))
#         config["Irrigation_cmd_area_shapefile"]["feature_name"] = "ADM3_EN"
#         config["Irrigation_cmd_area_shapefile"]["numeric_id_name"] = "ID"
#         config["Irrigation_cmd_area_shapefile"]["area_column_name"] = "Shape_Area"
#     else:
#         raise ValueError("Irrigation_cmd_area_shapefile section not found in config.")
    
#     if "GEE_Asset_ID" in config:
#         config["GEE_Asset_ID"]["shp"] = str(Path(save_loc).joinpath("Shapefiles", "Rupasi", "Rupasi.shp"))
#     else:
#         raise ValueError("GEE_Asset_ID section not found in config.")
        
#     if "Date_Running" in config:
#         config["Date_Running"]["start_date"] = "ADM3_EN"
#         config["Date_Running"]["default_run_week"] = "false"
#         config["Date_Running"]["run_week"] = ["currentweek"]
#     else:
#         raise ValueError("Date_Running or Cmd_Area_Config section not found in config.")

#     if "Cmd_Area_Config" in config:
#         config["Cmd_Area_Config"]["path"] = str(Path(save_loc).joinpath("config_files", "test", "ca_config.yaml"))
#     else:
#         raise ValueError("Cmd_Area_Config section not found in config.")

#     with open(config_path, 'w') as f:
#         yaml.dump(config, f)

# run_sdrips('config_files/test/ca_config.yaml')




# def compare_rasters(raster1_path: Path, raster2_path: Path, n_samples: int = 5, tol: float = 1e-2, random_seed: int = 47):
#     """
#     Compares two rasters by shape, non-NODATA values, and randomly sampled locations,
#     rounding values to 2 decimal places.

#     Args:
#         raster1_path (str or Path): Path to the first raster (.tif)
#         raster2_path (str or Path): Path to the second raster (.tif)
#         n_samples (int): Number of random points to sample for comparison
#         tol (float): Tolerance for comparing raster values (after rounding)
#         random_seed (int): Seed for reproducibility

#     Returns:
#         bool: True if all checks pass, False otherwise
#     """
#     results = {}

#     # Read rasters
#     with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
#         arr1 = src1.read(1).astype(float)
#         arr2 = src2.read(1).astype(float)

#         nodata1 = src1.nodata
#         nodata2 = src2.nodata

#     # Check shape
#     if arr1.shape != arr2.shape:
#         return False

#     # Mask NODATA values
#     mask = np.ones(arr1.shape, dtype=bool)
#     if nodata1 is not None:
#         mask &= arr1 != nodata1
#     if nodata2 is not None:
#         mask &= arr2 != nodata2

#     # Round arrays to 2 decimals and compare non-NODATA values
#     arr1_rounded = np.round(arr1[mask], 2)
#     arr2_rounded = np.round(arr2[mask], 2)
#     if not np.all(np.isclose(arr1_rounded, arr2_rounded, atol=tol)):
#         return False

#     # Randomly sample points and compare
#     np.random.seed(random_seed)
#     valid_indices = np.argwhere(mask)
#     n_samples = min(n_samples, len(valid_indices))
#     sampled_indices = valid_indices[np.random.choice(len(valid_indices), n_samples, replace=False)]

#     for idx in sampled_indices:
#         i, j = idx
#         if not np.isclose(round(arr1[i, j], 2), round(arr2[i, j], 2), atol=tol):
#             return False

#     return True


# def check_raster_exists_and_nonempty(raster_path: Path):
#     """
#     Checks if a raster file exists and is not empty.

#     Args:
#         raster_path (str or Path): Path to the raster (.tif) file

#     Returns:
#         bool: True if file exists and is non-empty, False otherwise
#     """
#     raster_path = Path(raster_path)

#     # Check existence
#     if not raster_path.exists():
#         print(f"Raster file does not exist: {raster_path}")
#         return False

#     # Check non-empty (file size > 0)
#     if raster_path.stat().st_size == 0:
#         print(f"Raster file is empty: {raster_path}")
#         return False

#     return True

# def compare_csv_expected(csv_expected_path: Path, csv_observed_path: Path, float_tol: float = 1e-6):
#     """
#     Compares an observed CSV against an expected CSV.
#     All expected columns must exist in the observed CSV and match in values.
#     Extra columns in observed CSV are ignored.
    
#     Args:
#         csv_expected_path (str or Path): Path to the expected CSV file
#         csv_observed_path (str or Path): Path to the observed CSV file
#         float_tol (float): Tolerance for comparing float values
    
#     Returns:
#         bool: True if all expected columns match, False otherwise
#     """
#     df_expected = pd.read_csv(csv_expected_path)
#     df_observed = pd.read_csv(csv_observed_path)
    
#     expected_cols = df_expected.columns
#     observed_cols = df_observed.columns
    
#     # Check if all expected columns exist in observed CSV
#     missing_cols = [col for col in expected_cols if col not in observed_cols]
#     if missing_cols:
#         print(f"Expected columns missing in observed CSV: {missing_cols}")
#         return False
    
#     # Compare values for expected columns only
#     for col in expected_cols:
#         series_exp = df_expected[col]
#         series_obs = df_observed[col]
        
#         if pd.api.types.is_numeric_dtype(series_exp):
#             if not np.allclose(series_exp.fillna(np.nan), series_obs.fillna(np.nan), atol=float_tol, equal_nan=True):
#                 print(f"Values do not match in column: {col}")
#                 return False
#         else:
#             if not series_exp.fillna("").equals(series_obs.fillna("")):
#                 print(f"Values do not match in column: {col}")
#                 return False
    
#     return True


# ## Comparing testing outputs with expected outputs

# ### 1. Comparing penman-monteith output
# actual_raster = 'outputs/penman_monteith.tif' # the directory should be with the test directory <project_root>/data/landsat/penman/currentweek/penman_eto_Phulpur.constant.tif could you pls change it?
# expected_raster = 'expected/penman_monteith.tif' # the directory should be the test
# assert compare_rasters(actual_raster, expected_raster), "Penman-Monteith output does not match expected."

# ### 2. Comparing sebal output
# actual_raster = 'outputs/sebal.tif'
# expected_raster = 'expected/sebal.tif'
# assert compare_rasters(actual_raster, expected_raster), "SEBAL output does not match expected."

# ### 3. Comparing irrigation output
# actual_raster = 'outputs/irrigation.tif'
# expected_raster = 'expected/irrigation.tif'
# assert compare_rasters(actual_raster, expected_raster), "Irrigation output does not match expected."

# ### 4. Comparing percolation output
# actual_raster = 'outputs/percolation.tif'
# expected_raster = 'expected/percolation.tif'
# assert compare_rasters(actual_raster, expected_raster), "Percolation output does not match expected."

# ### 5. Testing if IMERG Precipitation is working fine. Not comparing with expected output as 
# file_path = 'outputs/precipitation.tif'
# assert check_raster_exists_and_nonempty(file_path), "IMERG Precipitation output is not valid."

# ### 6. Testing if GFS Precipitation is working fine. Not comparing with expected output as 
# file_path = 'outputs/gfs_precipitation.tif'
# assert check_raster_exists_and_nonempty(file_path), "GFS Precipitation output is not valid."

# ### 7. Comparing command stats CSV file
# actual_csv = 'outputs/command_stats.csv'
# expected_csv = 'expected/command_stats.csv'
# assert compare_csv_expected(expected_csv, actual_csv), "Command stats CSV output does not match expected."

import sys, os
from pathlib import Path
import requests
import multiprocessing
import pandas as pd
import numpy as np
import rasterio
from ruamel.yaml import YAML
import shutil
import argparse


from sdrips.run_sdrips import run_sdrips
from sdrips.utils.utils import get_gdrive_url


CONFIG_SOURCES = {
    "ca_config.yaml": "1yyqTOThok0btGbiXhdBWNdHX3Rs0ir7O",
    "Rupasi.shx": "1i_TZHaB3qCZv48RymaufQV6-M1RjY5O3",
    "Rupasi.shp": "1hg-UAmFITlKWQxCG4oN65B70JQODP0TZ", 
    "Rupasi.dbf": "10bILvNFqTE398WJouzxdePVbKfo1Z8Og", 
    "Rupasi.prj": "1q-zLkgyCmd0MTmESlO1RpdbsjZt498rO", 
    "Rupasi.cpg": "1a6MORP0CFMGU9x1IrbNzgj-Agov_iyQ6",
    "penman_eto_Rupasi.constant.tif": "146rvhSpHpJhLgXzB4UlwjzKJysN3K0Z4",
    "sebal_eto_Rupasi.eta.tif": "1lbOdSIq87PeOHoUYHcQClkoEntRwntlO",
    "irrigation_Rupasi.eta.tif": "1d2GO9Eh1DAJKcRJ9ZhAbPVdqJEVUYa_L",
    "Percolation_currentweek.csv": "1C2BluSX2J3lRB4jp0yiM94UpxKhpQJ8E",
    "Landsat_Command_Area_Stats.csv": "1KkTHKbVTyATiYmsrNkVbITeyl61LewaD"
}

# ------------------------------
# Download utility
# ------------------------------
def download_test_files(test_data_dir, config_dir, shapefiles_dir, outputs_dir, csv_dir):
    """
    Downloads all expected outputs and config files from Google Drive
    into their respective test_data subfolders.
    """
    for file_name, file_id in CONFIG_SOURCES.items():
        try:
            gdrive_url = get_gdrive_url(file_id)
            response = requests.get(gdrive_url)
            if response.status_code != 200:
                print(f"Failed to download {file_name}")
                continue

            # Determine save directory
            if file_name.endswith(".yaml"):
                save_dir = config_dir
            elif file_name.endswith((".shp", ".shx", ".dbf", ".prj", ".cpg", ".geojson")):
                save_dir = shapefiles_dir
            elif file_name.endswith(".tif") or 'Percolation' in file_name:
                save_dir = outputs_dir
            elif file_name.endswith(".csv") and not 'Percolation' in file_name:
                save_dir = csv_dir
            else:
                save_dir = test_data_dir.joinpath("misc")

            save_dir.mkdir(parents=True, exist_ok=True)
            save_path = save_dir.joinpath(file_name)
            with open(save_path, "wb") as f:
                f.write(response.content)
            print(f"Downloaded {file_name} -> {save_path}")
        except Exception as e:
            print(f"Error downloading {file_name}: {e}")

# ------------------------------
# Copy utility
# ------------------------------

def copy_config_files_from_project(config_dir):
    """
    Copy all config files from the main project directory to test config directory
    """
    project_config_dir = Path("config_files")
    if not project_config_dir.exists():
        print(f"Project config directory not found: {project_config_dir}")
        return False
    
    # Copy all files from project config to test config
    copied_files = []
    for config_file in project_config_dir.glob("*"):
        if config_file.is_file():
            dest_path = config_dir.joinpath(config_file.name)
            shutil.copy2(config_file, dest_path)
            copied_files.append(config_file.name)
            
    return True

def copy_config_files_from_project(project_root_dir: Path, config_dir: Path):
    """
    Copy all config files from the main project directory to test config directory
    
    Args:
        project_root_dir (Path): Path to the main project root directory
        config_dir (Path): Path to the test config directory where files should be copied
    """
    project_config_dir = project_root_dir.joinpath("config_files")
    if not project_config_dir.exists():
        print(f"Project config directory not found: {project_config_dir}")
        return False
    
    # copied_files = []
    for config_file in project_config_dir.glob("*"):
        if config_file.is_file():
            dest_path = config_dir.joinpath(config_file.name)
            shutil.copy2(config_file, dest_path)
            # copied_files.append(config_file.name)
            print(f"Copied {config_file.name} to test config directory")
    
    # print(f"Copied {len(copied_files)} config files")
    return True

# ------------------------------
# YAML config updater
# ------------------------------
def update_config_file_for_test(config_path: Path, save_outputs_dir: Path, test_data_dir: Path):
    """
    Updates the config file for testing:
    - Sets save_data_loc to an external directory (outside expected results)
    - Updates config file paths to use test data from expected results folder
    
    Args:
        config_path (Path): Path to the config file to modify
        save_outputs_dir (Path): External directory where outputs should be saved
        test_data_dir (Path): Test data directory (expected results folder) containing config files
    """
    yaml = YAML()
    yaml.preserve_quotes = True

    # Define subdirectories within the test data directory
    config_dir = test_data_dir.joinpath("config_files")
    shapefiles_dir = test_data_dir.joinpath("Shapefiles")

    with open(config_path, "r") as f:
        config = yaml.load(f)

    # Save data location - point to EXTERNAL directory
    config["Save_Data_Location"] = {"save_data_loc": str(save_outputs_dir)}

    # path to roi in shapefiles folder
    roi_path = shapefiles_dir.joinpath("Rupasi.shp")
    if "Irrigation_cmd_area_shapefile" in config:
        config["Irrigation_cmd_area_shapefile"]["path"] = str(roi_path)
        config["Irrigation_cmd_area_shapefile"]["feature_name"] = "ADM4_EN"
        config["Irrigation_cmd_area_shapefile"]["numeric_id_name"] = "FID"
        config["Irrigation_cmd_area_shapefile"]["area_column_name"] = "Shape_Area"

    if "GEE_Asset_ID" in config:
        config["GEE_Asset_ID"]["shp"] = str(roi_path)

    ca_config_file_path = config_dir.joinpath("ca_config.yaml")
    if "Cmd_Area_Config" in config:
        config["Cmd_Area_Config"]["path"] = str(ca_config_file_path)

    crop_config_file_path = config_dir.joinpath("crop_config.yaml")
    if "Crop_Config" in config:
        config["Crop_Config"]["path"] = str(crop_config_file_path)

    config_links_file_path = config_dir.joinpath("config_links.yaml")
    if "Config_links" in config:
        config["Config_links"]["path"] = str(config_links_file_path)

    secrets_file_path = config_dir.joinpath("secrets.yaml")
    if "Secrets_Path" in config:
        config['Secrets_Path']["path"] = str(secrets_file_path)

    # Date Running
    if "Date_Running" in config:
        config['Date_Running']["start_date"] = "2023-02-14"
        config['Date_Running']["default_run_week"] = False
        config['Date_Running']["run_week"] = ["currentweek"] 

    # High-Level Controls:
    config["Command_Area_Net_Water_Requirement"] = True
    config["Canal_water_allotment"] = False
    config["Insitu_Sensor_integration"] = False
    config["Weather_station_integration"] = False

    with open(config_path, "w") as f:
        yaml.dump(config, f)
    
    return True


# ------------------------------
# Raster & CSV comparison utils
# ------------------------------
def compare_rasters(raster1_path: Path, raster2_path: Path, tol: float = 1e-2, n_samples: int = 5, random_seed: int = 47):
    """
    Compare two rasters by shape, non-NODATA values, and random sample points.

    Args:
        raster1_path (Path): Path to the first raster file.
        raster2_path (Path): Path to the second raster file.
        tol (float): Tolerance for comparing raster values.
        n_samples (int): Number of random samples to compare.
        random_seed (int): Seed for random number generation.

    Returns:
        bool: True if rasters are similar, False otherwise.
    """

    # Check if files exist first
    if not raster1_path.exists():
        raise FileNotFoundError(f"Actual raster file not found: {raster1_path}")
    if not raster2_path.exists():
        raise FileNotFoundError(f"Expected raster file not found: {raster2_path}")
    
    if not check_raster_exists_and_nonempty(raster1_path):
        raise ValueError(f"Actual raster file is empty or unreadable: {raster1_path}")
    if not check_raster_exists_and_nonempty(raster2_path):
        raise ValueError(f"Expected raster file is empty or unreadable: {raster2_path}")
    
    # Check if files are readable
    try:
        with rasterio.open(raster1_path) as src1, rasterio.open(raster2_path) as src2:
            arr1, arr2 = src1.read(1).astype(float), src2.read(1).astype(float)
            nodata1, nodata2 = src1.nodata, src2.nodata

        if arr1.shape != arr2.shape:
            print(f"Shape mismatch: {arr1.shape} vs {arr2.shape}")
            return False

        mask = np.ones(arr1.shape, dtype=bool)
        if nodata1 is not None:
            mask &= arr1 != nodata1
        if nodata2 is not None:
            mask &= arr2 != nodata2

        if not np.allclose(np.round(arr1[mask], 2), np.round(arr2[mask], 2), atol=tol):
            print(f"Value mismatch beyond tolerance {tol}")
            return False

        # Random sampling
        valid_idx = np.argwhere(mask)
        
        if len(valid_idx) == 0:
            print("No valid data points to compare")
            return False
        
        np.random.seed(random_seed)
        sampled_idx = valid_idx[np.random.choice(len(valid_idx), min(n_samples, len(valid_idx)), replace=False)]
        for i, j in sampled_idx:
            if not np.isclose(round(arr1[i, j], 2), round(arr2[i, j], 2), atol=tol):
                print(f"Sample mismatch at position ({i}, {j}): {arr1[i, j]} vs {arr2[i, j]}")
                return False
        return True
    
    except Exception as e:
        print(f"Error comparing rasters {raster1_path} and {raster2_path}: {e}")
        return False

def check_raster_exists_and_nonempty(raster_path: Path):
    """
    Checks if a raster file exists and is not empty.

    Args:
        raster_path (str or Path): Path to the raster (.tif) file

    Returns:
        bool: True if file exists and is non-empty, False otherwise
    """
    raster_path = Path(raster_path)

    # Check existence
    if not raster_path.exists():
        print(f"Raster file does not exist: {raster_path}")
        return False

    # Check non-empty (file size > 0)
    if raster_path.stat().st_size == 0:
        print(f"Raster file is empty: {raster_path}")
        return False

    return True

def compare_csv(csv_expected: Path, csv_observed: Path, tol: float = 1e-6):
    """
    Compares an observed CSV file against an expected CSV file.
    All columns in the expected CSV must exist in the observed CSV and match in values.
    Extra columns in the observed CSV are ignored.

    Columns containing 'PPT' or the net_water_req columns are excluded from comparison. 
    This is because IMERG data links updates frequently, which may change the test pipeline; 
    the comparison checks the IMERG module is functioning correctly without directly checking its values. 
    Similarly, net water requirement columns are excluded, as precipitation directly affects these values.

    Args:
        csv_expected_path (str or Path): Path to the expected CSV file
        csv_observed_path (str or Path): Path to the observed CSV file
        float_tol (float): Tolerance for comparing float values
    
    Returns:
        bool: True if all expected columns match, False otherwise
    """
    try:
        df_exp = pd.read_csv(csv_expected)
        df_obs = pd.read_csv(csv_observed)

        exclude_cols = ["net_water_req", "net_water_req_mm", "net_water_req_m3"]
        cols_to_check = [
            col for col in df_exp.columns
            if "PPT" not in col and col not in exclude_cols
        ]
        for col in cols_to_check:
            if pd.api.types.is_numeric_dtype(df_exp[col]):
                if not np.allclose(
                    df_exp[col].fillna(np.nan),
                    df_obs[col].fillna(np.nan),
                    atol=tol,
                    equal_nan=True
                ):
                    print(f"Mismatch in column {col}")
                    return False
            else:
                if not df_exp[col].fillna("").equals(df_obs[col].fillna("")):
                    print(f"Mismatch in column {col}")
                    return False
        return True
    
    except Exception as e:
        print(f"Error comparing CSV files {csv_expected} and {csv_observed}: {e}")
        return False

# ------------------------------
# Main test runner
# ------------------------------
def run_tests(test_dir: Path | str = "./tests") -> None:
    """
    Run sDRIPS verification tests to validate installation or for developer testing.

    The test workflow excludes comparisons of precipitation outputs and net water 
    requirement columns in the stats CSV. IMERG data links update frequently, which may 
    change the test pipeline; the comparison therefore verifies that the IMERG module 
    is functioning correctly without directly checking its values. Similarly, net water 
    requirement columns are excluded, as they are directly affected by precipitation.

    Args:
        test_dir (Path): Path to the test directory within the project (default: "./tests").

    Workflow:
        1. Downloads test files.
        2. Copies configuration files from the project.
        3. Updates configuration for testing.
        4. Runs sDRIPS.
        5. Compares raster and CSV outputs against expected results.
    """
    test_dir = Path(test_dir)

    if not test_dir.exists():
        # raise FileNotFoundError(f"Test directory does not exist: {test_dir}")
        os.makedirs(test_dir, exist_ok=True)
    
    # ------------------------------
    # Constants for test data
    # ------------------------------
    TEST_DATA_DIR = test_dir.joinpath("test_expected_data")
    CONFIG_DIR = TEST_DATA_DIR.joinpath("config_files")
    SHAPEFILES_DIR = TEST_DATA_DIR.joinpath("Shapefiles")
    OUTPUTS_DIR = TEST_DATA_DIR.joinpath("expected_outputs")
    CSV_DIR = TEST_DATA_DIR.joinpath("expected_field_stats")

    print("\nDownloading test files...")
    download_test_files(TEST_DATA_DIR, CONFIG_DIR, SHAPEFILES_DIR, OUTPUTS_DIR, CSV_DIR)

    print("\nCopying config files from project...")
    project_root = test_dir.parent
    if not copy_config_files_from_project(project_root, CONFIG_DIR):
        print("Warning: Could not copy all config files from project")
        return False
    
    test_config = CONFIG_DIR.joinpath("sdrips_config.yaml")
    # if not update_config_file_for_test(config_path = test_config, ):
    #     print("Failed to update config file")
    #     return False
    test_model_data_dir = test_dir.joinpath("Data")
    update_config_file_for_test(
        config_path=test_config,
        save_outputs_dir = test_model_data_dir,  
        test_data_dir=TEST_DATA_DIR 
    )

    
    print("\nRunning sDRIPS for testing...")
    run_sdrips(test_config)

    print("\nComparing outputs with expected...")
    raster_tests = [
        {
            "actual_filename": "penman_eto_Rupasi.constant.tif",
            "expected_filename": "penman_eto_Rupasi.constant.tif", 
            "actual_path": (project_root.joinpath('Data/landsat/penman/currentweek/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "sebal_eto_Rupasi.eta.tif",
            "expected_filename": "sebal_eto_Rupasi.eta.tif",
            "actual_path": (project_root.joinpath('Data/landsat/sebal/currentweek/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "irrigation_Rupasi.eta.tif", 
            "expected_filename": "irrigation_Rupasi.eta.tif",
            "actual_path": (project_root.joinpath('Data/landsat/irrigation/currentweek/')),
            "expected_path": OUTPUTS_DIR
        }
    ]

    for test in raster_tests:
        actual = test["actual_path"].joinpath(test["actual_filename"])
        expected = test["expected_path"].joinpath(test["expected_filename"])
        assert compare_rasters(actual, expected), f"Raster mismatch: {test['actual_filename']}"
    
    print("\nTests (3 Tests) for ET based raster outputs completed successfully.")
    
    # Check IMERG and GFS outputs exist and non-empty
    imerg_output = (project_root.joinpath('Data/precip/precip.currentweek.tif'))
    gfs_output = (project_root.joinpath('Data/precip/precip.nextweek.tif'))
    assert check_raster_exists_and_nonempty(imerg_output), "IMERG output is missing or empty"
    assert check_raster_exists_and_nonempty(gfs_output), "GFS output is missing or empty"

    print("\nChecks (2 Checks) for precipitation raster outputs completed successfully.")

    # Check CSV outputs
    csv_tests = [
        {
            "actual_filename": "Percolation_currentweek.csv",
            "expected_filename": "Percolation_currentweek.csv",
            "actual_path": (project_root.joinpath('Data/percolation/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "Landsat_Command_Area_Stats.csv", 
            "expected_filename": "Landsat_Command_Area_Stats.csv",
            "actual_path": (project_root.joinpath('Data/')),
            "expected_path": CSV_DIR
        }
    ]
    
    for test in csv_tests:
        actual = test["actual_path"].joinpath(test["actual_filename"])
        expected = test["expected_path"].joinpath(test["expected_filename"])
        assert compare_csv(expected, actual), f"CSV mismatch: {test['actual_filename']}"

    print("\nTests (2 Tests) for CSV outputs completed successfully.")

    print("\nAll tests passed!")

# ------------------------------
# CLI / Notebook entry
# ------------------------------
if __name__ == "__main__":
    # run_tests()
    parser = argparse.ArgumentParser(description="Run sDRIPS verification tests")
    parser.add_argument(
        "test_dir",
        type=Path,
        nargs="?", 
        default=Path("./tests"),
        help="Path to the test directory inside the project (e.g. ./tests)"
    )
    args = parser.parse_args()
    run_tests(args.test_dir)