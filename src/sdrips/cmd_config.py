import argparse
from pathlib import Path
from typing import List, Union
import geopandas as gpd
import yaml


def get_ca_ids(shapefile_path: Union[str, Path], column_name: str) -> List[str]:
    """Read shapefile and return unique IDs from the specified column."""
    gdf = gpd.read_file(shapefile_path)
    
    if column_name not in gdf.columns:
        raise ValueError(f"Column '{column_name}' not found in the shapefile.\n Available columns: {gdf.columns.tolist()}")
    
    return gdf[column_name].astype(str).unique().tolist()


def create_yaml_file(
    ca_ids: List[str],
    yaml_file_path: Union[str, Path],
    default_planting_date: str,
    default_crop_type: str,
    default_soil_coef: float,
    default_distribution_unif: float
) -> None:
    """Create a YAML config file for the command areas with defaults."""
    config = {
        'DEFAULT': {
            'use_default': True,
            'planting_date': default_planting_date,
            'crop_type': default_crop_type,
            'soil_coef': default_soil_coef,
            'distribution_unif': default_distribution_unif,
        }
    }

    for ca_id in ca_ids:
        config[ca_id] = {
            'use_default': False,
            'planting_date': default_planting_date,
            'crop_type': default_crop_type,
            'soil_coef': default_soil_coef,
            'distribution_unif': default_distribution_unif,
        }

    with open(yaml_file_path, 'w') as yaml_file:
        yaml.dump(config, yaml_file, sort_keys=False, default_flow_style=False)


def main():
    parser = argparse.ArgumentParser(
        description="Generate command area config YAML using the shapefile."
    )
    parser.add_argument(
        "-s", "--shp_path", required=True,
        help="Path to the shapefile. Use quotes if the path has spaces."
    )
    parser.add_argument(
        "-c", "--column_name", required=True,
        help="Column name in shapefile containing unique command area IDs."
    )
    parser.add_argument(
        "-d", "--default_planting_date", default="2023-04-01",
        help="Default planting date in YYYY-MM-DD format."
    )
    parser.add_argument(
        "-cc", "--default_crop_type", default="Rice",
        help="Default crop type (e.g., Rice, Wheat, Corn)."
    )
    parser.add_argument(
        "-sc", "--default_soil_coef", type=float, default=0.5,
        help="Default soil coefficient (e.g., 0.3, 0.5)."
    )
    parser.add_argument(
        "-du", "--default_distribution_unif", type=float, default=1.0,
        help="Default distribution uniformity (e.g., 0.7, 1.0)."
    )
    parser.add_argument(
        "-o", "--output_path", default=None,
        help="Optional output path for the YAML config file."
    )


    args = parser.parse_args()

    # Set output YAML path and ensure the folder exists
    if args.output_path:
        output_path = Path(args.output_path)
    else:
        output_path = Path("config_files/ca_config.yaml")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Get unique command area IDs
    ca_ids = get_ca_ids(args.shp_path, args.column_name)

    # Create the YAML file
    create_yaml_file(
        ca_ids=ca_ids,
        yaml_file_path=output_path,
        default_planting_date=args.default_planting_date,
        default_crop_type=args.default_crop_type,
        default_soil_coef=args.default_soil_coef,
        default_distribution_unif=args.default_distribution_unif,
    )

    print(f"Command area config file created at: {output_path.resolve()}")


if __name__ == "__main__":
    main()
