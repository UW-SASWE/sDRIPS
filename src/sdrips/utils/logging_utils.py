import logging
import logging.handlers
from logging.handlers import QueueHandler, QueueListener
import multiprocessing
import os
import datetime
from ruamel.yaml import YAML
from multiprocessing.queues import Queue
from typing import Tuple, Dict, List

from sdrips.utils.utils import load_yaml_config


class ExcludeModulesFilter(logging.Filter):
    def __init__(self, excluded_modules):
        super().__init__()
        self.excluded_modules = excluded_modules 

    def filter(self, record):
        return not any(record.pathname.endswith(mod) for mod in self.excluded_modules)

class IncludeModulesFilter(logging.Filter):
    def __init__(self, included_modules):
        super().__init__()
        self.included_modules = included_modules
    
    def filter(self, record):
        # Check if the filename is in included modules
        filename = record.filename
        return filename in self.included_modules 
    
def setup_logger_with_queue(config_path: str) -> Tuple[Queue, QueueListener, str]:
    """
    Setup logging configuration with a multiprocessing queue for log messages.
    This function initializes a logger that writes to a file, excluding logs from specific modules.
    Args:
        config_path (str): Path to the main configuration file.

    Returns:
        Tuple[Queue, QueueListener, str]: A tuple containing the log queue, queue listener, and the log file path.
    """
    script_config = load_yaml_config(config_path)
    save_data_loc = script_config['Save_Data_Location']['save_data_loc']
    os.makedirs(f'{save_data_loc}/logs/', exist_ok=True)
    dt_fmt = '%Y%m%d_%H%M%S'
    log_file = os.path.abspath(f'{save_data_loc}/logs/{datetime.datetime.today().strftime(dt_fmt)}.log')

    log_queue = multiprocessing.Queue()

    file_handler = logging.FileHandler(log_file)
    formatter = logging.Formatter('%(asctime)s - %(filename)s:%(lineno)d - %(message)s')
    file_handler.setFormatter(formatter)

    # excluded_modules = ['discovery.py', 'connectionpool.py', 'env.py', 
    # '__init__.py', 'warp.py', 'mask.py', 
    # 'utils.py', 'features.py', 'collection.py',
    # 'collection.py', 'geodataframe.py', 'retry.py'
    # 'google_auth_httplib2.py', 'session.py', 'font_manager.py'
    # 'file.py', 'google_auth_httplib2.py', 'api.py'
    # ]
    included_modules = [
        'sdrips_init.py',
        'sdrips_test.py',              
        'ee_utils.py',   
        'utils.py',
        'clean.py',
        'make_directory.py',
        'sensor_utils.py',
        'evapotranspiration.py',          
        'gfs_processing.py',             
        'run_sdrips.py',
        'canal_water_dist.py',
        'canal_water_distribution.py',
        'cmd_config.py',
        'percolation.py',
        'precipitation.py',
        'stn_corr_evapotranspiration.py',
        'sen_corr_evapotranspiration.py',
        'tiff_preprocess.py',
        'cmd_area_stats.py',
        'initialize.py',
        'logging_utils.py',
        'weather_station_utils.py',
        'command_area_creation.py'

        # Add any other files you want to monitor
    ]
    # file_handler.addFilter(ExcludeModulesFilter(excluded_modules))
    file_handler.addFilter(IncludeModulesFilter(included_modules))

    queue_listener = QueueListener(log_queue, file_handler)
    queue_listener.start()

    return log_queue, queue_listener, log_file
    
def worker_logger_setup(log_queue: multiprocessing.Queue) -> None:
    """
    Setup logger for worker processes using a multiprocessing queue.
    This function configures the logger to use a QueueHandler that sends log messages
    to a shared queue for processing by the main process.

    Args:
        log_queue (multiprocessing.Queue): Shared log queue for multiprocessing.

    Returns:
        None
    """
    qh = logging.handlers.QueueHandler(log_queue)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers = []
    logger.addHandler(qh)

def worker_init(log_queue: multiprocessing.Queue) -> None:
    """
    Initialize logger for worker processes.

    Args:
        log_queue (multiprocessing.Queue): Shared log queue for multiprocessing.
    
    Returns:
        None
    """
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    logger.handlers = []
    logger.addHandler(QueueHandler(log_queue))
