# Design Philosophy

sDRIPS was designed not only as an irrigation advisory model, but as a practical, reusable, and adaptable decision-support framework for satellite-informed surface-water irrigation management.

Earlier versions of the framework were developed as customized implementations for different regions, users, and operational settings. These included farmer-focused advisories, canal-operator advisories, sensor-integrated landscape irrigation, and farm-scale irrigation scheduling. While these versions demonstrated the flexibility of the approach, they also highlighted the need for a single centralized package that could reduce duplicated effort, simplify maintenance, and allow users to configure the system for their own region of interest.

The current sDRIPS package addresses this need through a design philosophy built around openness, modularity, scalability, cloud computing, data adaptability, minimal user input, and fail-safe operation.

![Key software design principles of sDRIPS](../images/Model/Design_Philosophy.png)

*Key software design principles of sDRIPS. The framework emphasizes open-source development, open-access data use, cloud computing, modularity, scalability, customizability, data adaptability, minimal user input, parallelization, efficiency, flexibility, and fail-safe operation.*

## Open-source and open-access design

sDRIPS is developed as an open-source Python package so that users can inspect, reproduce, modify, and extend the workflow. The framework relies on publicly available Earth observation and weather datasets where possible, reducing dependence on costly proprietary data sources. This design makes sDRIPS more accessible to researchers, irrigation managers, farmers, and decision-support practitioners working in both data-rich and resource-constrained environments.

sDRIPS is also designed to be easy to understand and maintain by following standard Python documentation practices, including [PEP 257](https://peps.python.org/pep-0257/) guidelines for clear and consistent docstrings.


## Cloud-based computation

A major design choice in sDRIPS is the use of cloud-based geospatial processing. Instead of requiring users to download and process large volumes of satellite imagery locally, sDRIPS uses cloud-hosted datasets and cloud computing where possible. Computationally intensive steps, such as satellite image preprocessing and evapotranspiration estimation, can be handled through cloud-based platforms. This lowers the demand on the user's local machine and allows sDRIPS to run on modest computing resources.

## Minimal user input

sDRIPS is designed to reduce the amount of information that users must provide manually. In a typical setup, the main user-provided inputs include:

| Required input | Purpose |
|---|---|
| Google Earth Engine account | Access cloud-based geospatial datasets |
| PPS account | Download IMERG precipitation data |
| Boundary polygon | Define the region of interest |
| Crop type and planting date | Estimate crop water demand |
| Sensor or weather-station data | Optional; used only when local data integration is enabled |

This design allows users to focus on the irrigation-management question rather than repetitive satellite-data preprocessing.

## Configuration-driven workflow

Most user interaction with sDRIPS occurs through configuration files rather than direct modification of the source code. This makes the package easier to use, reproduce, and adapt across regions. The configuration system controls high-level model behavior, including whether the user wants to generate farmer advisories, command-area net water requirement estimates, canal water allotment outputs, in situ sensor integration, or local weather-station integration. Lower-level settings define file paths, project directories, time windows, and execution details.

## Modular architecture

sDRIPS separates the irrigation advisory workflow into modular components:

1. Configuration setup
2. Evapotranspiration estimation
3. Percolation estimation
4. Precipitation integration
5. Advisory generation

Because these components are modular, users can activate, deactivate, or modify selected parts of the workflow without rebuilding the entire system. For example, one user may run sDRIPS using only satellite and global model data, while another user may enable local sensor or weather-station integration.

## Data adaptability

Different regions have different levels of data availability. Some users may have local weather stations or in situ soil moisture sensors, while others may only have a boundary polygon and crop information. sDRIPS was designed to operate across this range of data conditions. When local data are available, they can supplement or bias-correct satellite and model-based estimates. When local data are unavailable, sDRIPS can still operate using globally available satellite and numerical weather model products.

## Tested on operational-scale irrigation data

sDRIPS is designed for both small field-level applications and larger irrigation-management systems. The framework was benchmarked using the Teesta Barrage Project in Bangladesh as an operational-scale case study. For more details, refer to the HESS paper by [Khan et al. (2026)](https://doi.org/10.5194/hess-30-3675-2026).

| Benchmark item | Value |
|---|---:|
| Number of command areas | 153 |
| Total analysis area | 1548.65 km² |
| Landsat pixels processed | 1,720,721 |
| Sentinel-1 pixels processed | 15,486,485 |
| Example run date | 28 February 2023 |

This example demonstrates why scalability, cloud computing, and parallelization are central to the design of sDRIPS. The package is intended not only for individual plots, but also for command areas and canal systems where many spatial units must be processed repeatedly.

## Scalable and efficient execution

sDRIPS supports parallel processing so that multiple farms, command areas, or spatial units can be processed more efficiently. If the number of cores is not specified, sDRIPS uses one fewer than the maximum number of available cores. In this benchmark, that means 7 cores on the Windows laptop with 8 available cores, and 63 cores on the UNIX cluster with 64 available cores.

In the Teesta Barrage Project benchmark, total runtime decreased from 66.6 minutes to 10.4 minutes on a Windows laptop when parallel processing was used. On a UNIX-based cluster, runtime decreased from 46.4 minutes to 6 minutes.

| Workflow step | Windows laptop, 1 core | Windows laptop, max (8 - 1) cores | UNIX cluster, 1 core | UNIX cluster, max (64 - 1) cores |
|---|---:|---:|---:|---:|
| Evapotranspiration estimation (Penman-Monteith and SEBAL) | 49.2 mins | 6.3 mins | 35.8 mins | 4.6 mins |
| Raster preparation (unzipping and renaming) | 2.5 mins | 32.4 sec | 5.9 sec | 5.8 sec |
| Precipitation data acquisition and processing (nowcast and forecast) | 25.5 sec | 12.9 sec | 5.5 sec | 1.7 sec |
| Percolation estimation | 11.2 mins | 1.5 mins | 9.8 mins | 37.6 sec |
| Net water requirement estimation | 3.2 mins | 1.9 mins | 35.5 sec | 37.3 sec |
| **Total run time** | **66.6 mins** | **10.4 mins** | **46.4 mins** | **6.0 mins** |

These results show that sDRIPS can be used on standard computers while also benefiting from higher-performance systems when available.

## Fail-fast and fail-safe operation

sDRIPS follows a fail-fast and fail-safe design philosophy.

**Fail-fast** means that errors are detected early and reported clearly. This helps users identify missing files, incorrect paths, credential problems, configuration issues, or data-access errors before a long computational run fails.

**Fail-safe** means that the system is designed to continue operating when optional data streams are unavailable. For example, in a sensor-integrated setup, sDRIPS can fall back to a satellite-only workflow if sensor data are missing or interrupted.

This is important for operational irrigation advisory systems, where users still need outputs even when some data sources fail.

## Flexible interaction modes

sDRIPS supports both command-line and notebook-based workflows.

The command-line interface is useful for routine operational runs, including scheduled weekly execution. Jupyter Notebook tutorials support learning, exploration, and research-oriented applications.

This dual interaction mode makes sDRIPS usable by different types of users, including practitioners, researchers, students, and developers.

## From model to decision-support framework

The design philosophy of sDRIPS is centered on translating satellite remote sensing and weather-model data into usable irrigation advisories.

Its value comes not only from the underlying evapotranspiration and water-balance calculations, but also from the way the package makes those calculations accessible, configurable, scalable, and operationally robust.

In this sense, sDRIPS is more than a computational model. It is a flexible decision-support framework for satellite-informed surface-water irrigation management.
