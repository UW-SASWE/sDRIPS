# Jupyter Notebook Interface of sDRIPS

This page shows the notebook-friendly version of the sDRIPS command workflows. Use these functions when you want to run one step at a time, inspect intermediate files, or combine sDRIPS with your own analysis code.

!!! tip "Notebook setup"
    Start Jupyter from the root of your sDRIPS project directory. This keeps relative paths consistent with the command-line examples.

## Workflow Map

| Task | Import | Function |
| --- | --- | --- |
| Initialize a project | `sdrips.cli.sdrips_init` | `initialize_project()` |
| Run verification tests | `sdrips.cli.sdrips_test` | `run_tests()` |
| Run sDRIPS | `sdrips.run_sdrips` | `run_sdrips()` |
| Build command-area config | `sdrips.cmd_config` | `run_cmd_config()` |
| Estimate command-area polygons | `sdrips.cmd_area_creation` | `create_cmd_area()` |

## Recommended Imports

```python
from pathlib import Path

from sdrips.cli.sdrips_init import initialize_project
from sdrips.cli.sdrips_test import run_tests
from sdrips.run_sdrips import run_sdrips
from sdrips.cmd_config import run_cmd_config
from sdrips.cmd_area_creation import create_cmd_area
```

!!! note "After installation or updates"
    If you installed or upgraded sDRIPS while Jupyter was already open, restart the kernel before importing the package.

## Initialize a Project

Create the standard project folders and default configuration files:

```python
initialize_project(project_dir=".", force=False)
```

The function creates:

```text
Data/
Shapefiles/
config_files/
```

and downloads:

```text
config_links.yaml
crop_config.yaml
sdrips_config.yaml
secrets.yaml
```

!!! warning "Force refresh"
    Set `force=True` only when you intentionally want to refresh existing configuration files from the default templates.

## Create Command-Area Polygons

Use this step only when command-area boundaries are not already available. The function reads a canal or line network and creates estimated command-area polygons using the target area column.

```python
created_path = create_cmd_area(
    shape_path="Shapefiles/canal_network/Irrigation_Canal_Networks.shp",
    column_name="CMD_Ar_m2",
    output_path="Shapefiles/created_cmd_area/created_cmd_area.shp",
)

created_path
```

| Parameter | Description |
| --- | --- |
| `shape_path` | Input canal or line-network vector file. |
| `column_name` | Column containing target area values, usually in square meters. |
| `output_path` | Output shapefile path. |

!!! warning "Use measured boundaries when available"
    `create_cmd_area()` estimates command-area polygons. If measured command-area or farm-field boundaries are available, those are preferred.

## Generate `ca_config.yaml`

Use `run_cmd_config()` to create command-area crop settings from a vector file.

```python
ca_config_path = run_cmd_config(
    shp_path="Shapefiles/cmd_area/Teesta_Command_Areas.shp",
    column_name="CNLNM",
    default_planting_date="2025-07-07",
    default_crop_type="Rice",
    default_soil_coef=0.5,
    default_distribution_unif=1.0,
    output_path="config_files/ca_config.yaml",
)

ca_config_path
```

| Parameter | Description |
| --- | --- |
| `shp_path` | Command-area vector file. Shapefile format is recommended for Google Earth Engine integration. |
| `column_name` | Column containing unique command-area IDs. |
| `default_planting_date` | Default planting date in `YYYY-MM-DD` format. |
| `default_crop_type` | Crop type available in `crop_config.yaml`, such as `Rice`, `Wheat`, or `Corn`. |
| `default_soil_coef` | Soil coefficient used for each generated command area unless edited later. |
| `default_distribution_unif` | Distribution uniformity used for each generated command area unless edited later. |
| `output_path` | Output YAML path. Default behavior writes to `config_files/ca_config.yaml`. |

Preview the generated file:

```python
print(Path(ca_config_path).read_text())
```

!!! tip "Edit per command area"
    The generated file includes a `DEFAULT` section and one section for each command area. You can manually adjust planting date, crop type, soil coefficient, or distribution uniformity before running sDRIPS.

## Inspect Configuration Files

Two utility loaders are useful in notebooks:

```python
from sdrips.utils.initialize import load_config
from sdrips.utils.utils import load_yaml_config
```

Use `load_config()` when you want attribute-style access:

```python
config = load_config("config_files/sdrips_config.yaml")

config.Save_Data_Location.save_data_loc
config.Date_Running.run_week
config.Multiprocessing.cores
```

Use `load_yaml_config()` when you want a standard dictionary:

```python
config_dict = load_yaml_config("config_files/sdrips_config.yaml")

config_dict["Save_Data_Location"]["save_data_loc"]
```

!!! note "Utility functions"
    The `sdrips.utils` package does not currently expose separate terminal commands. Its functions are best used from notebooks, scripts, or internal modules.

## Check Crop Coefficients

You can inspect the crop coefficient values that sDRIPS will use for a crop.

```python
from sdrips.utils.utils import read_crop_coefficients, get_growth_kc

coefficients = read_crop_coefficients("config_files/crop_config.yaml", crop="Rice")
kc_day_35 = get_growth_kc(coefficients, num_days=35)

kc_day_35
```

This is useful when checking whether a planting date and crop type are producing reasonable growth-stage values.

## Earth Engine Utilities

Most users can let the main sDRIPS workflow initialize Earth Engine automatically. For notebooks, you can also initialize Earth Engine explicitly.

```python
from sdrips.utils.ee_utils import initialize_earth_engine

initialize_earth_engine(
    service_account="your-service-account@project.iam.gserviceaccount.com",
    key_file="config_files/path_to_key.json",
)
```

!!! tip "Credential setup"
    Store credential paths and account information in `secrets.yaml` for routine sDRIPS runs. Direct initialization in notebooks is mainly useful for debugging or exploratory work.

## Run Verification Tests

Run the standard verification workflow:

```python
run_tests(test_dir=Path("./tests"), sensor_test=False)
```

Run the sensor-corrected test case:

```python
run_tests(test_dir=Path("./tests"), sensor_test=True)
```

The test workflow downloads expected files, prepares test configuration, runs sDRIPS, and compares selected raster and CSV outputs.

## Run sDRIPS

After editing `sdrips_config.yaml`, run:

```python
run_sdrips("config_files/sdrips_config.yaml")
```

The function uses the module controls in `sdrips_config.yaml`, including ET estimation, precipitation, weather, percolation, command-area statistics, sensor correction, weather-station correction, and canal water allotment.

Typical outputs are written under the configured `Save_Data_Location`:

```text
Data/
  Landsat/
  percolation/
  logs/
  Landsat_Command_Area_Stats.csv
```

!!! tip "Track the run"
    `run_sdrips()` prints the log-file path at the end of the workflow. Use that file to diagnose missing credentials, failed downloads, or module-specific errors.

## Notebook Help

Jupyter can show function signatures and docstrings directly:

```python
run_cmd_config?
create_cmd_area?
run_sdrips?
```

For a fuller walkthrough, see the notebook tutorials:

| Tutorial | Use it when |
| --- | --- |
| [Creating command areas](/en/latest/tutorials/cmd_area_creation/Creating_CMD_Config/) | You need to estimate command-area polygons from canal data. |
| [Creating command-area config files](/en/latest/tutorials/config_files/Creating_CMD_Config/) | You need to generate `ca_config.yaml` interactively. |
