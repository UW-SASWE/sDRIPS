"""
Initialize configuration variables from the YAML config file
and return them as a SimpleNamespace object for easy attribute access.
"""
from ruamel.yaml import YAML
from types import SimpleNamespace
from typing import Any


yaml = YAML(typ='safe')
yaml.preserve_quotes = True

def load_config(config_path: str = "script_config.yaml") -> SimpleNamespace:
    """
    Load and parse the YAML configuration file into a SimpleNamespace.

    Args:
        config_path (str): Path to the YAML config file.

    Returns:
        SimpleNamespace: A namespace object containing config parameters as attributes.
    """
    with open(config_path, "r") as file:
        raw_config: dict[str, Any] = yaml.load(file)

    # Convert nested dictionaries into SimpleNamespaces recursively
    def dict_to_namespace(obj: Any) -> Any:
        if isinstance(obj, dict):
            return SimpleNamespace(**{k: dict_to_namespace(v) for k, v in obj.items()})
        elif isinstance(obj, list):
            return [dict_to_namespace(item) for item in obj]
        else:
            return obj

    return dict_to_namespace(raw_config)
