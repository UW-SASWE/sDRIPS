# Command Area Configuration File
The **command area configuration file** (`ca_config.yaml`) defines **region-specific management and cropping conditions**. It records **crop type, planting date, soil stress coefficient and water distribution efficiency** within each field boundary. This enables sDRIPS to model multiple fields or command areas in parallel, even when they differ in cropping patterns or water delivery characteristics. A representative snippet of the command area configuration file is shown below.

```yaml
DEFAULT:
  use_default: true
  planting_date: '2023-01-15'
  crop_type: Rice
  soil_coef: 0.5
  distribution_unif: 1.0
Teesta_Canal_A_1:
  use_default: false
  planting_date: '2023-01-15'
  crop_type: Rice
  soil_coef: 0.5
  distribution_unif: 1.0
T1T_2:
  use_default: false
  planting_date: '2023-01-15'
  crop_type: Rice
  soil_coef: 0.5
  distribution_unif: 1.0
```
## Parameter Definitions
- `planting_date`: The date when seeds were sown in the field.
- `crop_type`: The crop cultivated in the specific field or command area.
- `soil_coef`: The soil stress coefficient, which defaults to 0.5 as per the study by [Bose et al. (2021)](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2020WR028654)
- `distribution_unif`: The water distribution efficiency in the field or command area, with a default of 1.0.

## How `use_default` Works

- **DEFAULT section**  

    - If `use_default` in the `DEFAULT` section is set to `true`, sDRIPS **ignores parameters specified for individual command areas** and applies the default settings **globally**.  
    - If set to `false`, sDRIPS considers parameters for each command area separately.  

- **Individual command area**  

    - If `use_default` is set to `true` for a specific command area, sDRIPS **ignores the parameters in that command area** and applies the values from the `DEFAULT` section.  
    - If `false`, sDRIPS uses the parameters specified for that command area.  

This mechanism reduces the burden on the user by minimizing manual updates across multiple command areas, especially for large-scale runs.  

!!! tip_note "Tip"

    - To apply a particular setting across all command areas, **modify the parameter in the `DEFAULT` section** instead of updating each command area individually.
    - A **step-by-step guide** for **creating a command area configuration file** is available in the tutorial section using [***CLI***](/en/latest/tutorials/config_files/cli/) and [***Jupyter Notebook***](/en/latest/tutorials/config_files/Creating_CMD_Config/).
