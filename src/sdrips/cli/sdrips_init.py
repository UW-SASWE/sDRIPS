import argparse
import requests
from pathlib import Path
import sys
from typing import Dict
from ruamel.yaml import YAML
import multiprocessing
from pathlib import Path

from sdrips.cli.sdrips_test import run_tests
from sdrips.run_sdrips import run_sdrips


# GitHub links
CONFIG_SOURCES = {
    "config_links.yaml": "https://raw.githubusercontent.com/UW-SASWE/sDRIPS/refs/heads/main/config_files_template/config_links.yaml", 
    "crop_config.yaml": "https://raw.githubusercontent.com/UW-SASWE/sDRIPS/refs/heads/main/config_files_template/crop_config.yaml",
    "sdrips_config.yaml": "https://raw.githubusercontent.com/UW-SASWE/sDRIPS/refs/heads/main/config_files_template/sdrips_config.yaml",
    "secrets.yaml": "https://raw.githubusercontent.com/UW-SASWE/sDRIPS/refs/heads/main/config_files_template/secrets.yaml"
}


def init_project_structure(project_path: Path) -> None:
    """Create sDRIPS directory structure"""
    project_path.joinpath("Data").mkdir(parents=True, exist_ok=True)
    project_path.joinpath("Shapefiles").mkdir(exist_ok=True)
    project_path.joinpath("config_files").mkdir(exist_ok=True)

def download_config(raw_url: str, dest_path: Path, force: bool = False) -> bool:
    """
    Download config from Github sDRIPS repository.

    Args:
        raw_url (str): The raw URL of the config file.
        dest_path (Path): The local path to save the downloaded config file.
        force (bool): Whether to force re-download the file if it exists.

    Returns:
        bool: True if download was successful, False otherwise.
    """
    if dest_path.exists() and not force:
        print(f" {dest_path.name} (use --force to update)")
        return True

    try:
        response = requests.get(raw_url, timeout=15)
        response.raise_for_status()

        dest_path.write_bytes(response.content)
        return True
        
    except Exception as e:
        print(f" Download failed: {str(e)}", file=sys.stderr)
        return False

def update_config_file(config_path: str | Path, save_loc=None):
    """
    Updates multiprocessing cores and save_data_loc in the YAML config.

    Args:
        config_path (str or Path): Path to the sdrips_config.yaml file.
        save_loc (str or Path, optional): Location to save data. Defaults to current directory.
    """
    config_path = Path(config_path)
    # save_loc = save_loc or str(Path.cwd())  

    yaml = YAML()
    yaml.preserve_quotes = True  
    try:

        with open(config_path, 'r') as f:
            config = yaml.load(f)

        # Update multiprocessing cores
        max_cores = max(1, multiprocessing.cpu_count() - 1)
        if "Multiprocessing" in config:
            config["Multiprocessing"]["cores"] = max_cores
        else:
            raise ValueError("Multiprocessing section not found in config.")

        # Update save_data_loc
        if "Save_Data_Location" in config:
            data_path = Path(str(save_loc)).joinpath('Data')
            config["Save_Data_Location"]["save_data_loc"] = str(data_path)
        else:
            config["Save_Data_Location"] = {"save_data_loc": str(data_path)}

        with open(config_path, 'w') as f:
            yaml.dump(config, f)

        print(f"sDRIPS Configurations updated: cores={max_cores}, save_data_loc={save_loc}")

    except Exception as e:
        print(f"Error updating config file: {e}.")

    
def initialize_project(project_dir: str = ".", force: bool = False) -> bool:
    """
    Core initialization workflow

    Args:
        project_dir (str): The project directory.
        force (bool): Whether to force re-initialization.
    """
    project_path = Path(project_dir).expanduser().resolve()
    
    try:
        print(f" Initializing project at: {project_path}")
        init_project_structure(project_path)
        
        print("\n Fetching latest configurations:")
        config_dir = project_path.joinpath("config_files")
        
        for filename, url in CONFIG_SOURCES.items():
            print(f"• {filename}...", end=" ", flush=True)
            dest_path = config_dir.joinpath(filename)
            success = download_config(url, dest_path, force)
            print("Configuration File has been downloaded successfully." if success else "Configuration File download failed.")

        update_config_file(config_dir.joinpath("sdrips_config.yaml"), save_loc=project_path)

        print(f"\n Project ready at: {project_path}")
        return True
        
    except Exception as e:
        print(f"\n Initialization failed: {str(e)}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="sDRIPS CLI")
    subparsers = parser.add_subparsers(dest="command")

    # init subcommand
    init_parser = subparsers.add_parser("init", help="Initialize a sDRIPS project",
                                        description=initialize_project.__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter
                                        )
    init_parser.add_argument(
        "--dir", "-d",
        default=".",
        help="Project directory (default: current)"
    )
    init_parser.add_argument(
        "--force", "-f",
        action="store_true",
        help="Force update all config files"
    )

    # ------------------------------
    # test subcommand
    # ------------------------------
    test_parser = subparsers.add_parser("test", help="Run sDRIPS verification tests",
                                        description=run_tests.__doc__, 
                                        formatter_class=argparse.RawDescriptionHelpFormatter
                                        )
    test_parser.add_argument(
        "--dir", "-d",
        type=Path,
        default=Path("./tests"),
        help="Path to the test directory (default: ./tests)"
    )
    test_parser.add_argument(
        "--sensor","-s",
        action="store_true",
        help="Enable sensor test mode"
    )

    # ------------------------------
    # run subcommand
    # ------------------------------
    run_parser = subparsers.add_parser("run", help="Run sDRIPS",
                                        description=run_sdrips.__doc__,
                                        formatter_class=argparse.RawDescriptionHelpFormatter
                                        )
    run_parser.add_argument(
        "--config", "-c",
        type=Path,
        default=Path("./config_files/sdrips_config.yaml"),
        help="Path to the config file (default: ./config_files/sdrips_config.yaml)"
    )

    args = parser.parse_args()

    if args.command == "init":
        sys.exit(0 if initialize_project(args.dir, args.force) else 1)
    
    elif args.command == "test":
        run_tests(args.dir, args.sensor)  

    elif args.command == "run":
        run_sdrips(args.config)

    else:
        parser.print_help()