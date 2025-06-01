import tempfile
import shutil
from pathlib import Path
import geopandas as gpd
import pytest
from shapely.geometry import Polygon

from sdrips.cmd_config import get_ca_ids, create_yaml_file

@pytest.fixture
def temp_shapefile():
    """Creates a temporary shapefile with fake CA_IDs"""
    temp_dir = tempfile.mkdtemp()
    shapefile_path = Path(temp_dir) / "test_ca.shp"

    data = {
        "CA_ID": ["CA101", "CA102", "CA103"],
        "geometry": [
            Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
            Polygon([(2, 0), (3, 0), (3, 1), (2, 1)]),
            Polygon([(4, 0), (5, 0), (5, 1), (4, 1)])
        ]
    }

    gdf = gpd.GeoDataFrame(data, crs="EPSG:4326")
    gdf.to_file(shapefile_path)

    yield shapefile_path

    shutil.rmtree(temp_dir)  

def test_create_yaml_file(temp_shapefile):
    # Setup
    ca_ids = get_ca_ids(temp_shapefile, "CA_ID")
    output_yaml = Path(tempfile.mkdtemp()) / "ca_config.yaml"

    # Run
    create_yaml_file(
        ca_ids=ca_ids,
        yaml_file_path=output_yaml,
        default_planting_date="2023-04-01",
        default_crop_type="Rice",
        default_soil_coef=0.5,
        default_distribution_unif=1.0
    )

    # Assert output YAML exists and contains CA IDs
    assert output_yaml.exists()
    content = output_yaml.read_text()
    assert "CA101" in content
    assert "crop_type: Rice" in content
