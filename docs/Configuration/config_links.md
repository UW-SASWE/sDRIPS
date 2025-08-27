# Links Configuration File
The **links configuration file** (`config_links.yaml`) centralizes and manages URLs for external data services, such as **IMERG precipitation datasets**. Maintaining these links in a configuration file simplifies updates and ensures that the pipeline remains compatible with evolving data access protocols. A representative snippet of the configuration file is shown below.
```yaml
precipitation: 
  base_urls:
    before: "https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/"
    year:   "https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/"
    after:  "https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early/"

  file_pattern:
    prefix: "3B-HHR-E.MS.MRG.3IMERG."
    suffix: "-S233000-E235959.1410.V07B.1day.tif"

```

## How to update the link
To ensure that the pipeline accesses the latest data, follow these steps to update the configuration:  

<span class="preparation_step">Step 1: Check the latest dataset version</span> <br>
Visit the official IMERG data directory through [link](https://jsimpsonhttps.pps.eosdis.nasa.gov/imerg/gis/early).  

<span class="preparation_step">Step 2: Identify the prefix and suffix</span> <br>
View the prefix and suffix by hovering over the file names, or by copying a file link and pasting it into the browser to check the prefix and suffix details.

- For example, lets say that prefix on the file link is `3B-HHR-E.MS.MRG.3IMERG.`
- Suffix for the file link is `-S213000-E235959.1410.V07D.3day.tif`
- In the suffix, pay particular attention to the **version segment** (`V07D`) part as this indicates the latest release.  

<span class="preparation_step">Step 3: Update the configuration file</span> <br>
Replace the old prefix or suffix with the updated values:

- suffix: `-S233000-E235959.1410.V07D.1day.tif`. 
- prefix: `3B-HHR-E.MS.MRG.3IMERG.`

!!! tip_note "Tip"

    - Only the **version** (`V07D`) needs to be updated; time (`233000`) and interval (`1day`) remain unchanged.
    - Ensure the structure of the file names, including periods (.) and other delimiters, matches the format in the dataset.