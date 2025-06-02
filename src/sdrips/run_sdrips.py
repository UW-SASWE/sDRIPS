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

from src.sdrips.evapotranspiration import process_cmd_area_parallel 
from src.sdrips.utils.utils import load_yaml_config
from src.sdrips.utils.logging_utils import (
    setup_logger_with_queue,
    worker_logger_setup,
    worker_init
)


def run_et_for_ca(save_data_loc: str) -> Dict[str, List[float]]:
    """
    Run ET estimation for all command areas.

    Args:
        save_data_loc (str): Directory to save data and logs.
    """
    log_queue, queue_listener, log_file = setup_logger_with_queue(save_data_loc)
    main_logger = logging.getLogger()
    main_logger.setLevel(logging.DEBUG)
    main_logger.addHandler(QueueHandler(log_queue))

    try:
        main_logger.info("Location of the Data Folder: " + save_data_loc)
        main_logger.info("Starting ET estimation...")

        penman_mean_values,sebal_mean_values,irr_mean_values = process_cmd_area_parallel(main_logger, log_queue)
    except Exception as e:
        main_logger.exception("Unhandled exception in main process")

    finally:
        queue_listener.stop()

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


if __name__ == '__main__':
    args = parse_args()
    config = load_yaml_config(args.config)

    save_data_loc = config['Save_Data_Location']['save_data_loc']
    ET_estimation = config['Run_ET_Estimation'].get('et_estimation', False)

    if ET_estimation:
        run_et_for_ca(save_data_loc)