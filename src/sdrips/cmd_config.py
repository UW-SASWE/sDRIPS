import argparse
from pathlib import Path
from typing import List, Union, Optional
import geopandas as gpd
import yaml
from ruamel.yaml import YAML


def get_ca_ids(
    data_path: Union[str, Path],
    column_name: str,
    layer: Optional[Union[str, int]] = None,
) -> List[str]:
    """
    Read a vector dataset and return unique IDs from the specified column.

    Supports any format readable by GeoPandas, including:
      - ESRI Shapefile (.shp) and zipped shapefiles (.zip)
      - GeoJSON (.geojson, .json)
      - GeoPackage (.gpkg)  [may require 'layer']
      - ESRI File Geodatabase (.gdb)  [requires 'layer']
      - Many others depending on your GDAL build
    Note on KML/KMZ: Reading KML/KMZ depends on your GDAL build and drivers. If you run into issues, it’s often easier to convert to GeoPackage or GeoJSON first.

    Parameters
    ----------
    data_path : str | Path
        Path to the vector dataset (file path, a .zip of a shapefile, or a
        directory like a .gdb). For multi‑layer datasets, provide `layer`.
    column_name : str
        Name of the attribute column containing unique IDs.
    layer : str | int | None, optional
        Layer name or index for multi‑layer sources (e.g., GeoPackage/FileGDB).
        Ignored for single‑layer files like .shp or .geojson.

    Returns
    -------
    List[str]
        Unique IDs as strings, in arbitrary order.

    Raises
    ------
    ValueError
        If the column isn't found, or the layer is missing/invalid.
    """
    try:
        gdf = gpd.read_file(data_path, layer=layer) if layer is not None else gpd.read_file(data_path)
    except Exception as e:
        # Provide a friendlier hint for common multi-layer/driver issues
        msg = (
            f"Failed to read vector data at '{data_path}'. "
            "If this is a multi‑layer dataset (e.g., .gpkg or .gdb), pass the 'layer' name or index. "
            f"Original error: {e}"
        )
        raise ValueError(msg) from e

    if column_name not in gdf.columns:
        raise ValueError(
            f"Column '{column_name}' not found. "
            f"Available columns: {gdf.columns.tolist()}"
        )

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
    
    yaml = YAML()
    yaml.indent(sequence=4, offset=2)
    yaml.default_flow_style = False

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
        # yaml.dump(config, yaml_file, sort_keys=False, default_flow_style=False)
        for key, value in config.items():
            yaml.dump({key: value}, yaml_file)
            yaml_file.write("\n")

def run_cmd_config(
    shp_path: Union[str, Path],
    column_name: str,
    default_planting_date: str = "2023-04-01",
    default_crop_type: str = "Rice",
    default_soil_coef: float = 0.5,
    default_distribution_unif: float = 1.0,
    output_path: Union[str, Path, None] = None,
    layer: Optional[Union[str, int]] = None
) -> Path:
    """
    Generate a YAML configuration file for command areas found in a vector dataset.

    This function reads a vector dataset (e.g., Shapefile, GeoJSON, GeoPackage, FileGDB)
    using GeoPandas, extracts the unique command‑area identifiers from the specified
    attribute column, and writes a YAML config with one section per ID plus a DEFAULT
    section containing baseline parameters. Here Shapefile format is highly recommended as it 
    can be integrated with Google Earth Engine easily. 

    Parameters
    ----------
    shp_path : str or pathlib.Path
        Path to the vector dataset. Can be a file (e.g., ``.shp``, ``.geojson``,
        ``.gpkg``) or a directory/container supported by your GDAL/Fiona build
        (e.g., a File Geodatabase ``.gdb``). Zipped shapefiles (``.zip``) are also
        supported by GeoPandas. Here ``.shp`` format is highly recommended as it 
        can be integrated with Google Earth Engine easily.
    column_name : str
        Name of the attribute column that contains the unique command‑area IDs to
        index the YAML sections (e.g., ``"CA_ID"``).
    default_planting_date : str, optional
        Default planting date to write into both the DEFAULT section and each
        per‑ID section (ISO format ``YYYY-MM-DD``). Default is ``"2023-04-01"``.
    default_crop_type : str, optional
        Default crop type (e.g., ``"Rice"``, ``"Wheat"``, ``"Corn"``). Written to
        the DEFAULT section and each per‑ID section. Default is ``"Rice"``.
    default_soil_coef : float, optional
        Default soil coefficient value (dimensionless) to include in the config.
        Default is ``0.5``.
    default_distribution_unif : float, optional
        Default distribution uniformity (DU) value (dimensionless, e.g., 0.7–1.0)
        representing irrigation application efficiency. Default is ``1.0``.
    output_path : str, pathlib.Path, or None, optional
        Where to write the YAML file. If ``None``, the file is written to
        ``config_files/ca_config.yaml``. Parent directories are created if missing.

    Returns
    -------
    pathlib.Path
        Path to the written YAML file.

    """
    output_path = Path(output_path) if output_path else Path("config_files/ca_config.yaml")
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ca_ids = get_ca_ids(shp_path, column_name, layer=layer)
    create_yaml_file(
        ca_ids=ca_ids,
        yaml_file_path=output_path,
        default_planting_date=default_planting_date,
        default_crop_type=default_crop_type,
        default_soil_coef=default_soil_coef,
        default_distribution_unif=default_distribution_unif,
    )
    return output_path


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
        help="Column name in shapefile containing unique command area names. It is recommended to use a column with string values and no spaces in between."
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
        help="[Optional] Default soil coefficient (e.g., 0.3, 0.5)."
    )
    parser.add_argument(
        "-du", "--default_distribution_unif", type=float, default=1.0,
        help="[Optional] Default distribution uniformity (e.g., 0.7, 1.0)."
    )
    parser.add_argument(
        "-o", "--output_path", default=None,
        help="[Optional] Output path for the YAML config file. Default is config_files/ca_config.yaml"
    )


    args = parser.parse_args()

    # # Set output YAML path and ensure the folder exists
    # if args.output_path:
    #     output_path = Path(args.output_path)
    # else:
    #     output_path = Path("config_files/ca_config.yaml")
    # output_path.parent.mkdir(parents=True, exist_ok=True)

    # # Get unique command area IDs
    # ca_ids = get_ca_ids(args.shp_path, args.column_name)

    # # Create the YAML file
    # create_yaml_file(
    #     ca_ids=ca_ids,
    #     yaml_file_path=output_path,
    #     default_planting_date=args.default_planting_date,
    #     default_crop_type=args.default_crop_type,
    #     default_soil_coef=args.default_soil_coef,
    #     default_distribution_unif=args.default_distribution_unif,
    # )

    output = run_cmd_config(
        shp_path=args.shp_path,
        column_name=args.column_name,
        default_planting_date=args.default_planting_date,
        default_crop_type=args.default_crop_type,
        default_soil_coef=args.default_soil_coef,
        default_distribution_unif=args.default_distribution_unif,
        output_path=args.output_path
    )
    print(f"Command area config file created at: {output.resolve()}")


if __name__ == "__main__":
    main()
