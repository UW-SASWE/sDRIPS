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
    process_cmd_area_parallel,
    get_cmd_area_list)
from sdrips.utils.utils import load_yaml_config
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


def run_et_for_ca(save_data_loc: str, log_queue) -> Dict[str, List[float]]:
    """
    Run ET estimation for all command areas.

    Args:
        save_data_loc (str): Directory to save data and logs.
    """
    worker_init(log_queue)
    logger = logging.getLogger(__name__)
    try:
        logger.info("Starting ET estimation...")

        penman_mean_values,sebal_mean_values,irr_mean_values = process_cmd_area_parallel(logger, log_queue)
        # cmd_area_list = process_cmd_area_parallel(save_data_loc, log_queue)
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


def main():
    """
    Main function to orchestrate sDRIPS execution.
    """
    args = parse_args()
    config = load_yaml_config(args.config)

    save_data_loc = config['Save_Data_Location']['save_data_loc']
    run_week = config['Date_Running']['run_week']
    clear_condition = config['Clean_Directory'].get('clear_directory_condition', False)
    run_et = config['Run_ET_Estimation'].get('et_estimation', False)
    cmd_area_list = get_cmd_area_list()

    log_queue, queue_listener, log_file_path = setup_logger_with_queue(save_data_loc)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.addHandler(QueueHandler(log_queue))

    logger.info("===== Starting sDRIPS Framework Execution =====")
    logger.info(f"Configuration loaded from: {args.config}")
    logger.info(f"Data directory: {save_data_loc}")

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
            cmd_area_list = run_et_for_ca(save_data_loc, log_queue)
            unzip_tiffs(save_data_loc, run_week)
            convert_tiffs(save_data_loc, run_week)
            converting_to_eto(cmd_area_list, save_data_loc, run_week)

    except Exception as e:
        logger.exception("Unhandled exception during sDRIPS execution")

    finally:
        logger.info("===== sDRIPS Framework Finished =====")
        queue_listener.stop()
        

if __name__ == '__main__':
    main()