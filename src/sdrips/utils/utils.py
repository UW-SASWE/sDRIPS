
import os
import logging
import datetime
from ruamel.yaml import YAML
from typing import Tuple, Dict, List

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

