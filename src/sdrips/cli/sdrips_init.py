import argparse
import requests
from pathlib import Path
import sys
from typing import Dict

# Google Drive File IDs 
CONFIG_SOURCES = {
    "config_links.yaml": "1MuT0Wa1O-PXI1GtkyNqeEkJfGVMWpVEE", 
    "crop_config.yaml": "1C6zQvCJYaED35aBeQhV51eyHmp_8EQLB",
    "sdrips_config.yaml": "1bdeUYhK9tQZx-suZpqfAV-4OrBwgIqdB",
    "secrets.yaml": "1FDXnPmHQxD71CmjiZwCJ8oSd4NLtlD-Y"
}

def get_gdrive_url(file_id: str) -> str:
    """Generate direct download URL for Google Drive file"""
    return f"https://drive.google.com/uc?export=download&id={file_id}"

def init_project_structure(project_path: Path) -> None:
    """Create sDRIPS directory structure"""
    project_path.joinpath("Data").mkdir(parents=True, exist_ok=True)
    project_path.joinpath("Shapefiles").mkdir(exist_ok=True)
    project_path.joinpath("config_files").mkdir(exist_ok=True)

def download_config(file_id: str, dest_path: Path, force: bool = False) -> bool:
    """Download config from Google Drive with virus scan handling (if prompts)"""
    if dest_path.exists() and not force:
        print(f" {dest_path.name} (use --force to update)")
        return True

    try:
        # First request (may get virus warning page)
        url = get_gdrive_url(file_id)
        session = requests.Session()
        response = session.get(url, timeout=15)
        response.raise_for_status()
        
        # Handle Google Drive's virus scan warning
        if "virus scan warning" in response.text:
            confirm_link = f"https://drive.google.com/uc?export=download&confirm=t&id={file_id}"
            response = session.get(confirm_link, timeout=15)
        
        with open(dest_path, "wb") as f:
            f.write(response.content)
        return True
    except Exception as e:
        print(f" Download failed: {str(e)}", file=sys.stderr)
        return False

def initialize_project(project_dir: str = ".", force: bool = False) -> bool:
    """Core initialization workflow"""
    project_path = Path(project_dir).expanduser().resolve()
    
    try:
        print(f" Initializing project at: {project_path}")
        init_project_structure(project_path)
        
        print("\n Fetching latest configurations:")
        config_dir = project_path.joinpath("config_files")
        
        for filename, file_id in CONFIG_SOURCES.items():
            print(f"• {filename}...", end=" ", flush=True)
            dest_path = config_dir.joinpath(filename)
            success = download_config(file_id, dest_path, force)
            print("Configuration Files have been downloaded successfully." if success else "Configuration Files download failed.")
        
        print(f"\n Project ready at: {project_path}")
        return True
        
    except Exception as e:
        print(f"\n Initialization failed: {str(e)}", file=sys.stderr)
        return False

def main():
    parser = argparse.ArgumentParser(description="sDRIPS CLI")
    subparsers = parser.add_subparsers(dest="command")

    # init subcommand
    init_parser = subparsers.add_parser("init", help="Initialize a sDRIPS project")
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

    args = parser.parse_args()

    if args.command == "init":
        sys.exit(0 if initialize_project(args.dir, args.force) else 1)
    else:
        parser.print_help()