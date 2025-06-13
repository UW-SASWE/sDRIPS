"""
Main script to run sD.R.I.P.S framework.

This script reads configuration from a YAML file, sets up logging (with multiprocessing support),
and processes all defined command areas based on the configuration file flags.

"""

import os
import logging
import datetime
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
from sdrips.initialize import load_config
from sdrips.canal_water_distribution import calculate_canal_cumulative_discharge


def run_et_for_ca(config_path: str, log_queue: Queue, max_workers: float) -> None:
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
            "specified in the script configuration YAML file.",
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
    canal_water_allotment = config.Canal_water_allotment

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
            logger.info("Running ET module for all command areas...")
            run_et_for_ca(config_path, log_queue, cores)
            unzip_tiffs(config_path)
            convert_tiffs(config_path)
            converting_to_eto(config_path)
        if run_precip:
            logger.info("Running Precipitation module...")
            imergprecip(config_path)
        if run_weather:
            logger.info("Running GFS module...")
            gfsdata(config_path)
        if run_soil_moisture:
            logger.info("Running Percolation module...")
            percolation_estimation(config_path)
        if run_region_stats:
            logger.info("Running Command Area Statistics module...")
            command_area_info(config_path)
        if canal_water_allotment:
            logger.info("Running Canal Water Allotment module...")
            calculate_canal_cumulative_discharge(config_path)

    except Exception as e:
        logger.exception("Unhandled exception during sDRIPS execution")

    finally:
        logger.info("===== sDRIPS Framework Finished =====")
        queue_listener.stop()
        
def main():
    args = parse_args()
    run_sdrips(args.config)


if __name__ == '__main__':
    main()