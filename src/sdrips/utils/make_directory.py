import os
import logging
from typing import List


def make_directory(save_data_loc: str, run_week: List[str]) -> None:
    """
    Create required directories for sDRIPS framework if they do not exist.

    Args:
        save_data_loc (str): Base path to create data folders under.
        run_week (List[str]): List of week identifiers (e.g., ['lastweek', 'currentweek']).
    """
    logger = logging.getLogger()
    folders_list = ['landsat', 'uploads', 'precip', 'logs']#, 'avgt', 'tmax', 'tmin', 'ugrd', 'vgrd', 'wind']
    et_params = ['penman', 'sebal', 'irrigation']

    try:
        for folder in folders_list:
            if folder == 'landsat':
                for param in et_params:
                    for week in run_week:
                        path = os.path.join(save_data_loc, folder, param, week)
                        os.makedirs(path, exist_ok=True)
            else:
                path = os.path.join(save_data_loc, folder)
                os.makedirs(path, exist_ok=True)

        logger.info("All necessary directories created (if not already present).")

    except Exception as e:
        logger.error("Error while executing make_directory():")
        logger.exception(e)
