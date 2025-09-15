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


CONFIG_SOURCES_rupasi = {
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

CONFIG_SOURCES_sensor = {
    "ca_config.yaml": "1B5XvmWzDlky__6ziOh_6nKjuP4qJGkp3",
    "UW_ROI.shx": "1h_joLyeDj-l9wRNYRotIvbJWJF1lTA1c",
    "UW_ROI.shp": "13xkxSSbaV2vX6BGIGbe_N5G-fS4NUDU3", 
    "UW_ROI.dbf": "1JXd-2MSkiPY0opI7CEKaqh4UOQN2hK29", 
    "UW_ROI.prj": "13JOwvKBMu3Dd_B5mzhntPKGzaehcwh8g", 
    "UW_ROI.cpg": "1WH0_pdAIOTmydNwVPB78Th60zXTEPaKS",
    "penman_eto_Seattle.constant.tif": "1YOgv_jCQPrXqPwkx3A-w7Yl69hheQj_V",
    "sebal_eto_Seattle.eta.tif": "1mO_mvivKoaryowN4zrO6meCa8i8Xrw9d",
    "irrigation_Seattle.eta.tif": "1tQnMlba67sXl5c7ZX9dIq4c1wq33epnW",
    "Percolation_currentweek.csv": "1Bx8n9cDMpS00XazXVGSAfSPnhaZQK_VE",
    "Landsat_Command_Area_Stats.csv": "1PZBpLr5sJFhh09HEDGa1ae-s23Ti0wjz",
    "Sensor_Types.csv": "1Dvo5ZV0j4dmNUA0R9wwYBYYAmwdiCxTu",
    "Sensor_Database.csv": "1-qkLE0ugIxBhugUSZnZxcTFM9G2SxoQW"
}

# ------------------------------
# Download utility
# ------------------------------
def download_test_files(test_data_dir, config_dir, shapefiles_dir, outputs_dir, csv_dir, sensor_test=False):
    """
    Downloads all expected outputs and config files from Google Drive
    into their respective test_data subfolders.
    """

    dictionary_to_use = CONFIG_SOURCES_sensor if sensor_test else CONFIG_SOURCES_rupasi
    for file_name, file_id in dictionary_to_use.items():
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

def move_sensor_files(csv_dir: Path, test_data_dir: Path):
    """
    Move test sensor meta data files from the downloaded directory to test data directory.
    """
    if not csv_dir.exists():
        print(f"Expected Field Stats directory not found: {csv_dir}")
        return False

    move_files = ["Sensor_Types.csv", "Sensor_Database.csv"]
    for sensor_file in csv_dir.glob("*"):
        if sensor_file.is_file() and sensor_file.name in move_files:
            dest_path = test_data_dir.joinpath("Sensor_Data",sensor_file.name)
            dest_path.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(sensor_file, dest_path)
            print(f"Moved {sensor_file.name} to test data directory")
            
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
def update_config_file_for_test(config_path: Path, save_outputs_dir: Path, test_data_dir: Path, sensor_test: bool = False):
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
    if not sensor_test:
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

    elif sensor_test:
        # Save data location - point to EXTERNAL directory
        config["Save_Data_Location"] = {"save_data_loc": str(save_outputs_dir)}

        # path to roi in shapefiles folder
        roi_path = shapefiles_dir.joinpath("UW_ROI.shp")
        if "Irrigation_cmd_area_shapefile" in config:
            config["Irrigation_cmd_area_shapefile"]["path"] = str(roi_path)
            config["Irrigation_cmd_area_shapefile"]["feature_name"] = "ADM3_EN"
            config["Irrigation_cmd_area_shapefile"]["numeric_id_name"] = "id"
            config["Irrigation_cmd_area_shapefile"]["area_column_name"] = "Irr_Area"

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
            config['Date_Running']["start_date"] = "2024-06-07"
            config['Date_Running']["default_run_week"] = False
            config['Date_Running']["run_week"] = ["currentweek"] 

        # LUCL condition
        config["GLCC_Mask"]['glcc_mask'] = False
        config["GLCC_Mask"]['glcc_id'] = "users/saswegee/shahzaib/Lidar_NDVI_Map_glcc"

        # High-Level Controls:
        config["Command_Area_Net_Water_Requirement"] = True
        config["Canal_water_allotment"] = False
        config["Insitu_Sensor_integration"] = True
        config["Weather_station_integration"] = False

        # Sensor details
        config['Insitu_Sensor_integration_config']["air_temperature_sensor"] = True
        config['Insitu_Sensor_integration_config']["air_temperature_variable"] = 'Temperature'
        config['Insitu_Sensor_integration_config']["soil_moisture_sensor"] = False
        config['Insitu_Sensor_integration_config']["wind_speed_sensor"] = True
        config['Insitu_Sensor_integration_config']["wind_speed_variable"] = 'Wind'
        config['Insitu_Sensor_integration_config']["specific_humidity_sensor"] = True
        config['Insitu_Sensor_integration_config']["specific_humidity_variable"] = 'Specific Humidity'
        config['Insitu_Sensor_integration_config']["hub_url"] = "https://depts.washington.edu/saswe/Sensor_Cronjob/"
        config['Insitu_Sensor_integration_config']["sensor_data_path"] = "Sensor_Data"


    with open(config_path, "w") as f:
        yaml.dump(config, f)
    
    return True

def update_config_links_file_for_test(config_path: Path):
    """
    Updates the precipitation related data links in config_links file for testing:
    
    Args:
        config_path (Path): Path to the config file to modify.
    """
    yaml = YAML()
    yaml.preserve_quotes = True

    with open(config_path, "r") as f:
        config = yaml.load(f)
    
    config["precipitation"]["base_urls"]["before"] = "https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/2024/"
    config["precipitation"]["file_pattern"]["prefix"] = "3B-HHR-L.MS.MRG.3IMERG."
    
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
        raise FileNotFoundError(f"Actual raster file not found (full path): {raster1_path.absolute()}")
    if not raster2_path.exists():
        raise FileNotFoundError(f"Expected raster file not found (full path): {raster2_path.absolute()}")

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
                    difference = df_exp[col].fillna(np.nan) - df_obs[col].fillna(np.nan)
                    print(f"Mismatch in column {col}. Differences: {difference}")
                    return False
            else:
                if not df_exp[col].fillna("").equals(df_obs[col].fillna("")):
                    print(f"Mismatch in column {col}.")
                    return False
        return True
    
    except Exception as e:
        print(f"Error comparing CSV files {csv_expected} and {csv_observed}: {e}")
        return False

# ------------------------------
# Main test runner
# ------------------------------
def run_tests(test_dir: Path | str = "./tests", sensor_test: bool = False) -> None:
    """
    Run sDRIPS verification tests to validate installation or for developer testing.

    The test workflow excludes comparisons of precipitation outputs and net water 
    requirement columns in the stats CSV. IMERG data links update frequently, which may 
    change the test pipeline; the comparison therefore verifies that the IMERG module 
    is functioning correctly without directly checking its values. Similarly, net water 
    requirement columns are excluded, as they are directly affected by precipitation.

    Args:
        test_dir (Path): Path to the test directory within the project (default: "./tests").
        sensor_test (bool): If True, runs the sensor test case; otherwise runs the Rupasi (satellite only) test case.

    Workflow:
        1. Downloads test files.
        2. Copies configuration files from the project.
        3. Updates configuration for testing.
        4. Runs sDRIPS.
        5. Compares raster and CSV outputs against expected results.
    """
    test_dir = Path(test_dir)
    print('Test directory:', test_dir.resolve())
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
    download_test_files(TEST_DATA_DIR, CONFIG_DIR, SHAPEFILES_DIR, OUTPUTS_DIR, CSV_DIR, sensor_test = sensor_test)

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
        test_data_dir=TEST_DATA_DIR,
        sensor_test = sensor_test
    )
    if sensor_test:
        update_config_links_file_for_test(config_path=CONFIG_DIR.joinpath("config_links.yaml"))
    else:
        pass
    print("\nMoving sensor meta files to test model data directory...")
    if not move_sensor_files(CSV_DIR, test_model_data_dir):
        print("Warning: Could not move all sensor files")
        return False
    
    print("\nRunning sDRIPS for testing...")
    run_sdrips(test_config)

    print("\nComparing outputs with expected...")
    raster_tests_rupasi = [
        {
            "actual_filename": "penman_eto_Rupasi.constant.tif",
            "expected_filename": "penman_eto_Rupasi.constant.tif", 
            "actual_path": (test_dir.joinpath('Data/landsat/penman/currentweek/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "sebal_eto_Rupasi.eta.tif",
            "expected_filename": "sebal_eto_Rupasi.eta.tif",
            "actual_path": (test_dir.joinpath('Data/landsat/sebal/currentweek/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "irrigation_Rupasi.eta.tif", 
            "expected_filename": "irrigation_Rupasi.eta.tif",
            "actual_path": (test_dir.joinpath('Data/landsat/irrigation/currentweek/')),
            "expected_path": OUTPUTS_DIR
        }
    ]

    raster_tests_uw = [
        {
            "actual_filename": "penman_eto_Seattle.constant.tif",
            "expected_filename": "penman_eto_Seattle.constant.tif", 
            "actual_path": (test_dir.joinpath('Data/landsat/penman/currentweek/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "sebal_eto_Seattle.eta.tif",
            "expected_filename": "sebal_eto_Seattle.eta.tif",
            "actual_path": (test_dir.joinpath('Data/landsat/sebal/currentweek/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "irrigation_Seattle.eta.tif", 
            "expected_filename": "irrigation_Seattle.eta.tif",
            "actual_path": (test_dir.joinpath('Data/landsat/irrigation/currentweek/')),
            "expected_path": OUTPUTS_DIR
        }
    ]
    raster_tests = raster_tests_rupasi if not sensor_test else raster_tests_uw
    for test in raster_tests:
        actual = test["actual_path"].joinpath(test["actual_filename"])
        expected = test["expected_path"].joinpath(test["expected_filename"])
        assert compare_rasters(actual, expected), f"Raster mismatch: {test['actual_filename']}"
    
    if not sensor_test:
        print("\nTests (3 Tests) for ET based raster outputs completed successfully.")
    else:
        print("\nTests (3 Tests) for sensor corrected ET based raster outputs completed successfully.")
    
    # Check IMERG and GFS outputs exist and non-empty
    imerg_output = (test_dir.joinpath('Data/precip/precip.currentweek.tif'))
    gfs_output = (test_dir.joinpath('Data/precip/precip.nextweek.tif'))
    assert check_raster_exists_and_nonempty(imerg_output), "IMERG output is missing or empty"
    assert check_raster_exists_and_nonempty(gfs_output), "GFS output is missing or empty"

    print("\nChecks (2 Checks) for precipitation raster outputs completed successfully.")

    # Check CSV outputs
    csv_tests = [
        {
            "actual_filename": "Percolation_currentweek.csv",
            "expected_filename": "Percolation_currentweek.csv",
            "actual_path": (test_dir.joinpath('Data/percolation/')),
            "expected_path": OUTPUTS_DIR
        },
        {
            "actual_filename": "Landsat_Command_Area_Stats.csv", 
            "expected_filename": "Landsat_Command_Area_Stats.csv",
            "actual_path": (test_dir.joinpath('Data/')),
            "expected_path": CSV_DIR
        }
    ]
    
    for test in csv_tests:
        actual = test["actual_path"].joinpath(test["actual_filename"])
        expected = test["expected_path"].joinpath(test["expected_filename"])
        assert compare_csv(expected, actual), f"CSV mismatch: {test['actual_filename']}"

    if not sensor_test:
        print("\nTests (2 Tests) for CSV outputs completed successfully.")
    else:
        print("\nTests (2 Tests) for sensor corrected CSV outputs completed successfully.")

    print("\nAll tests passed!")

# ------------------------------
# CLI / Notebook entry
# ------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sDRIPS verification tests")
    parser.add_argument(
        "test_dir",
        type=Path,
        nargs="?", 
        default=Path("./tests"),
        help="Path to the test directory inside the project (e.g. ./tests)"
    )
    parser.add_argument(
        "--sensor",
        action="store_true",
        help="Enable sensor test mode"
    )
    
    args = parser.parse_args()
    run_tests(args.test_dir, args.sensor)