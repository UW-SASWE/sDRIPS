
import os
import sys
import logging
import datetime
from ruamel.yaml import YAML
from typing import Tuple, Any, Dict, List, Union

yaml = YAML(typ='safe')
yaml.preserve_quotes = True

def load_yaml_config(config_path: str) -> Dict:
    """
    Load YAML configuration file.

    Args:
        config_path (str): Path to YAML config file.

    Returns:
        dict: Parsed configuration dictionary.
    """
    with open(config_path, 'r') as f:
        return yaml.load(f)

def read_cmd_area_settings(cmd_area_data: Dict[str, Any], defaults: Dict[str, Any]) -> Dict[str, Any]:
    """
    Read settings for a specific command area with fallback to defaults.

    Args:
        cmd_area_data (Dict[str, Any]): Specific command area configuration.
        defaults (Dict[str, Any]): Default values for planting_date, crop_type, soil_coef.

    Returns:
        Dict[str, Any]: Final resolved settings for the command area.
    """
    if cmd_area_data.get('use_default', False):
        return {
            'planting_date': defaults['planting_date'],
            'crop_type': defaults['crop_type'],
            'soil_coefficient': defaults['soil_coef']
        }
    else:
        return {
            'planting_date': cmd_area_data.get('planting_date', defaults['planting_date']),
            'crop_type': cmd_area_data.get('crop_type', defaults['crop_type']),
            'soil_coefficient': cmd_area_data.get('soil_coef', defaults['soil_coef'])
        }


def read_crop_coefficients(config_path: str, crop: str) -> Dict[Tuple[int, int], Union[float, Tuple]]:
    """
    Read crop coefficient (Kc) values for a given crop from YAML config.

    Args:
        config_path (str): Path to the YAML file containing crop Kc values.
        crop (str): Crop name (case-insensitive).

    Returns:
        Dict[Tuple[int, int], Union[float, Tuple]]: Mapping from day ranges to Kc or interpolation info.
    """
    crop_config = load_yaml_config(config_path)
    normalized_crop_config = {k.lower(): v for k, v in crop_config.items()}
    crop_key = crop.lower()

    if crop.lower() in ['grass', 'alfa-alfa']:
        logging.info(f"Using constant Kc=1.0 for crop '{crop}'")
        return {(0, 999): 1.0}  
    
    if crop_key not in normalized_crop_config:
        logging.error(f"Crop '{crop}' not found in crop config file: {config_path}")
        sys.exit(1)

    crop_data = normalized_crop_config.get(crop_key, {})
    coefficients = {}

    for key, value in crop_data.items():
        try:
            start_day, end_day = map(int, key.split('-'))
        except ValueError:
            logging.error(f"Invalid day range key '{key}' in crop config for crop '{crop}'")
            sys.exit(1)
        if isinstance(value, list) and value[0] == 'linear':
            # Format: ['linear', start_value, end_value, num_days]
            try:
                coefficients[(start_day, end_day)] = (
                    'linear', float(value[1]), float(value[2]), int(value[3])
                )
            except (ValueError, IndexError):
                logging.error(f"Invalid linear format for range {key} in crop '{crop}'")
                sys.exit(1)
        else:
            try:
                coefficients[(start_day, end_day)] = float(value)
            except ValueError:
                logging.error(f"Non-numeric Kc value for range {key} in crop '{crop}'")
                sys.exit(1)

    return coefficients


def get_growth_kc(coefficients: Dict[Tuple[int, int], Union[float, Tuple]], num_days: int) -> Union[float, None]:
    """
    Get crop coefficient (Kc) based on crop growth stage and number of days since planting.

    Args:
        coefficients (dict): Dictionary mapping day ranges to Kc values or linear tuples.
        num_days (int): Days since planting.

    Returns:
        float or None: Crop coefficient for the given day or None if undefined.
    """
    for day_range, kc in coefficients.items():
        if day_range[0] <= num_days <= day_range[1]:
            if isinstance(kc, tuple) and kc[0] == 'linear':
                _, start_val, end_val, days = kc
                days_into_range = num_days - day_range[0]
                delta_per_day = (end_val - start_val) / days
                return start_val + delta_per_day * days_into_range
            return kc
    return None  # If num_days is outside defined ranges