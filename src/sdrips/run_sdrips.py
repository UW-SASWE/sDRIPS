"""
Main script to run sD.R.I.P.S framework.

This script reads configuration from a YAML file, sets up logging (with multiprocessing support),
and processes all defined command areas based on the configuration file flags.

"""

import os
import logging
import re
import datetime
import time
from ruamel.yaml import YAML
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from logging.handlers import QueueHandler, QueueListener
import multiprocessing
from multiprocessing.queues import Queue
from typing import Tuple, Dict, List

from sdrips.utils.clean import clear_data
from sdrips.utils.make_directory import make_directory
from sdrips.evapotranspiration import (
    process_cmd_area_parallel
)
from sdrips.utils.utils import (
    load_yaml_config,
    get_cmd_area_list,
    get_irrigation_cmd_area    
)
from sdrips.utils.logging_utils import (
    setup_logger_with_queue,
    worker_logger_setup,
    worker_init
)
from sdrips.tiff_preprocess import (
    unzip_tiffs, 
    convert_tiffs, 
    converting_to_eto
)
from sdrips.precipitation import imergprecip
from sdrips.gfs_processing import gfsdata
from sdrips.percolation import percolation_estimation
from sdrips.cmd_area_stats import command_area_info
from sdrips.utils.initialize import load_config
from sdrips.canal_water_distribution import calculate_canal_cumulative_discharge
from sdrips.sen_corr_evapotranspiration import (
    process_cmd_area_sensor_parallel
)
from sdrips.stn_corr_evapotranspiration import (
    process_cmd_area_stn_parallel
)


def run_et_for_ca(config_path: str, log_queue: Queue, max_workers: float, cmd_area_nwr: bool = True, Insitu_Sensor_integration: bool = False, Weather_station_integration: bool = False) -> None:
    """
    Run ET estimation for all command areas.

    Args:
        save_data_loc (str): Directory to save data and logs.
        log_queue (multiprocessing.Queue): Queue for logging messages.
        max_workers (float): Maximum number of worker processes to use.

    Returns:
        None
    """
    worker_init(log_queue)
    logger = logging.getLogger(__name__)
    try:
        logger.info("Starting ET estimation...")
        logger.info(f"Configuration - Insitu : {Insitu_Sensor_integration}, Weather Station : {Weather_station_integration}, Command Area NWR : {cmd_area_nwr}")
        if Insitu_Sensor_integration == True and cmd_area_nwr == True:
            process_cmd_area_sensor_parallel(config_path, logger, log_queue, max_workers)
        elif Weather_station_integration == True and cmd_area_nwr == True:
            process_cmd_area_stn_parallel(config_path, logger, log_queue, max_workers)
        elif cmd_area_nwr == True :
            process_cmd_area_parallel(config_path, logger, log_queue, max_workers)
        

    except Exception as e:
        logger.exception("Unhandled exception in main process")

def parse_args():
    """
    Parse command-line arguments using argparse.

    Returns:
        argparse.Namespace: Parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Runs the sDRIPS framework for surface water irrigation optimization.\n\n"
            "This script supports both full and modular runs (with multiprocessing support) depending on the flags "
            "specified in the sdrips configuration YAML file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '-c', '--config',
        required=True,
        help='Path to the script configuration YAML file.'
    )
    return parser.parse_args()


def run_sdrips(config_path: str):
    """
    Main function to orchestrate sDRIPS execution.

    Args:
        config_path (str): Path to the main configuration file.

    Returns:
        None
    """
    print('Reading Configurations and Setting Up Logging...')
    config = load_config(config_path)

    save_data_loc = config.Save_Data_Location.save_data_loc
    run_week = config.Date_Running.run_week
    cores = config.Multiprocessing.cores if config.Multiprocessing.cores is not None else multiprocessing.cpu_count() - 1
    clear_condition = config.Clean_Directory.clear_directory_condition
    run_et = config.Run_ET_Estimation.et_estimation
    run_precip = config.Precipitation_Config.consider_preciptation
    run_weather = config.Weather_Config.consider_forecasted_weather
    run_soil_moisture = config.Percolation_Config.consider_percolation
    run_region_stats = config.Region_stats.estimate_region_stats
    cmd_area_nwr = config.Command_Area_Net_Water_Requirement
    canal_water_allotment = config.Canal_water_allotment
    Insitu_Sensor_integration = config.Insitu_Sensor_integration
    Weather_station_integration = config.Weather_station_integration

    log_queue, queue_listener, log_file_path = setup_logger_with_queue(config_path)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(QueueHandler(log_queue))

    logger.info("===== Starting sDRIPS Framework Execution =====")
    logger.info(f"Configuration loaded from: {config_path}")
    logger.info(f"Data directory: {save_data_loc}")
    logger.info(f"Number of Cores Used: {cores}")
    try:
        make_directory(save_data_loc, run_week)
        if clear_condition:
            logger.info("Initiating directory cleaning...")
            clear_data(
                save_data_loc,
                clear_et=True,
                clear_precip=True,
                clear_weather=True,
                clear_uploads=True
            )

        if run_et:
            start_time = time.time()
            logger.info("Running ET module for all command areas...")
            run_et_for_ca(config_path, log_queue, cores, cmd_area_nwr, Insitu_Sensor_integration, Weather_station_integration)
            logger.info(f"ET module finished in {time.time() - start_time:.2f} seconds")

            start_time = time.time()
            unzip_tiffs(config_path)
            logger.info(f"Tiff unzipping finished in {time.time() - start_time:.2f} seconds")

            start_time = time.time()
            convert_tiffs(config_path)
            logger.info(f"Tiff conversion finished in {time.time() - start_time:.2f} seconds")

            start_time = time.time()
            converting_to_eto(config_path)
            logger.info(f"Conversion to ETo finished in {time.time() - start_time:.2f} seconds")

        if run_precip:
            start_time = time.time()
            logger.info("Running Precipitation module...")
            imergprecip(config_path)
            logger.info(f"Precipitation module finished in {time.time() - start_time:.2f} seconds")

        if run_weather:
            start_time = time.time()
            logger.info("Running GFS module...")
            gfsdata(config_path)
            logger.info(f"GFS module finished in {time.time() - start_time:.2f} seconds")
            
        if run_soil_moisture:
            start_time = time.time()
            logger.info("Running Percolation module...")
            percolation_estimation(config_path)
            logger.info(f"Percolation module finished in {time.time() - start_time:.2f} seconds")

        if run_region_stats:
            start_time = time.time()
            logger.info("Running Command Area Statistics module...")
            command_area_info(config_path)
            logger.info(f"Command Area Statistics module finished in {time.time() - start_time:.2f} seconds")

        if canal_water_allotment:
            start_time = time.time()
            logger.info("Running Canal Water Allotment module...")
            calculate_canal_cumulative_discharge(config_path)
            logger.info(f"Canal Water Allotment module finished in {time.time() - start_time:.2f} seconds")
        
        logger.info("===== sDRIPS Framework Finished =====")
        print(f'sDRIPS execution completed. Check logs for details. Address of the log: {log_file_path}')

    except Exception as e:
        logger.exception("Unhandled exception during sDRIPS execution")
        print(f'Error occurred during sDRIPS execution. \nCheck logs for details. Address of the log: {log_file_path}')

    finally:
        # logger.info("===== sDRIPS Framework Finished =====")
        queue_listener.stop()
        
        
def main():
    args = parse_args()
    run_sdrips(args.config)


if __name__ == '__main__':
    main()