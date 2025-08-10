# Conceptual Model of sDRIPS
sDRIPS is a cloud-based, open-source irrigation advisory system designed to optimize surface water use by integrating Earth observation datasets and numerical weather models. It leverages publicly available data from multiple sources, including:

1. **Landsat** and **Sentinel-1** satellites for surface monitoring
2. **Global Forecast System (GFS)** for meteorological forecasts
3. **Global Precipitation Measurement (GPM)** – IMERG for precipitation estimates. 

At its core, sDRIPS estimates evapotranspiration (ET) using two models:

1. **Penman–Monteith** <a href="https://www.fao.org/4/x0490e/x0490e00.htm" target="_blank">(Allen et al., 1998)</a> method – used as a proxy for potential crop water demand under optimal (non-stressed) conditions.

2. **Surface Energy Balance Algorithm for Land (SEBAL)** <a href="https://doi.org/10.1016/S0022-1694(98)00253-4" target="_blank">(Bastiaanssen et al., 1998)</a> – used to estimate actual ET under prevailing field conditions.

By comparing these two ET estimates, sDRIPS quantifies the gap between water needed by crops and water actually supplied. This information is then combined with - nowcast and forecast precipitation data and estimated percolation losses to calculate the net irrigation requirement.

sDRIPS applies both energy balance and water balance principles to generate actionable irrigation advisories at multiple spatial scales — from individual farms to large irrigation command areas. 

![sDRIPS Concept Illustration](../images/sDRIPS_Illustration.png)
<p style="text-align:center"> <sup><i>Illustration of how sDRIPS works. Green represents the regions where surplus irrigation has been provided, while red regions are deficit regions needing more irrigation.</i></sup> </p>


## Evapotranspiration 
### Penman-Monteith ET
The **Penman–Monteith** evapotranspiration (\(ETo\)) represents the **potential water demand** for a well-watered reference crop. In sDRIPS, \(ETo\) is computed for each Landsat overpass date and across all regions of interest, following the methodology of <a href="https://www.fao.org/4/x0490e/x0490e00.htm" target="_blank">(Allen et al., 1998)</a>.  

The **FAO-56 Penman–Monteith** equation is:  

$$ ETo = \frac{0.408\,\Delta\,(R_n - G) + \gamma \,\frac{900}{T + 273} \, u_2 \,(e_s - e_a)}{\Delta + \gamma \,(1 + 0.34\,u_2)}$$
    
Where:  

- $ETo$ — Reference evapotranspiration (mm/day)  
- $R_n$ — Net radiation at the crop surface (MJ m$^{-2}$ day$^{-1}$)  
- $G$ — Soil heat flux density (MJ m$^{-2}$ day$^{-1}$)  
- $T$ — Mean daily air temperature at 2 m height (°C)  
- $u_2$ — Wind speed at 2 m height (m/s)  
- $e_s$ — Saturation vapor pressure (kPa)  
- $e_a$ — Actual vapor pressure (kPa)  
- $(e_s - e_a)$ — Vapor pressure deficit (kPa)  
- $\Delta$ — Slope of the vapor pressure curve (kPa °C$^{-1}$)  
- $\gamma$ — Psychrometric constant (kPa °C$^{-1}$)  

The $ETo$ computed above represents the water use of a reference crop under optimal (non-stressed) conditions. To estimate crop-specific potential evapotranspiration ($ET$), sDRIPS applies crop and water stress coefficients:  

$$
ET = ETo \times K_c \times K_s
$$

Where:  

- $ET$ — Potential evapotranspiration for the specific crop (mm/day)  
- $K_c$ — Crop coefficient (varies by crop type and growth stage)  
- $K_s$ — Soil water stress coefficient (accounts for reduced transpiration under water-limited conditions)  

### SEBAL ET
The SEBAL algorithm, developed by <a href="https://doi.org/10.1016/S0022-1694(98)00253-4" target="_blank">Bastiaanssen et al. (1998)</a>, has been successfully implemented in more than 30 countries <a href="https://doi.org/10.1061/(ASCE)0733-9437(2005)131:1(85)" target="_blank">(Bastiaanssen et al., 2005)</a>. It has demonstrated high accuracy in estimating evapotranspiration.

SEBAL solves the surface energy balance equation to estimate evapotranspiration using satellite imagery and meteorological forcing data. It computes an instantaneous ET flux at the time of the Landsat overpass. The model calculates net surface radiation, soil heat flux, and sensible heat flux to the air. The residual energy flux, which is the latent heat flux driving evapotranspiration, is obtained by subtracting the soil and sensible heat fluxes from net radiation at the surface. For each pixel, the latent heat flux \( \lambda E \) is given by the surface energy balance equation:

\[
\lambda E = R_n - G - H
\]

Where:

- \( \lambda E \) is the latent heat flux (energy used for evapotranspiration),  
- \( R_n \) is net radiation at the surface,  
- \( G \) is soil heat flux, and  
- \( H \) is sensible heat flux.

The SEBAL-based evapotranspiration (SEBAL ET) serves as a proxy for the **actual water consumed by the crop**, a concept validated by <a href="https://doi.org/10.1029/2020WR028654" target="_blank">Bose et al. (2021)</a>. In this study, daily (24-hour) evapotranspiration was estimated assuming that instantaneous ET variations are negligible over the day <a href="https://doi.org/10.1061/(ASCE)0733-9437(2007)133:4(380)" target="_blank">(Allen et al., 2007)</a>. For weekly ET, the evapotranspiration estimated on the Landsat acquisition day was assumed steady for the following seven days until the next image was available.


## Percolation
To account for percolation, Sentinel-1 Synthetic Aperture Radar (SAR) data available on Google Earth Engine (GEE) was utilized. Sentinel-1 C-band (5 cm wavelength) data has been extensively used for soil moisture estimation and shows promising results up to 100 mm depth (<a href="https://doi.org/10.1109/TGRS.2018.2858004" target="_blank">Bauer-Marschallinger et al., 2019</a>; <a href="https://doi.org/10.1016/j.asr.2022.03.019" target="_blank">Bhogapurapu et al., 2022</a>). Ground sensors are ideal but impractical for large-scale deployment. Hence, Sentinel-1 was selected to maintain global scalability.

After estimating soil moisture from Sentinel-1, field capacity soil moisture was derived from the <a href="https://doi.org/10.5281/zenodo.2784001" target="_blank">Hengl & Gupta (2019)</a> dataset available on GEE. Percolation was then computed as:

\[
Percolation =
\begin{cases}
\text{Soil Moisture} - \text{Field Capacity}, & \text{if } \geq 0 \\
0, & \text{else}
\end{cases}
\]

## Precipitation
Accurate estimation of net water requirement depends on incorporating precipitation events. Precipitation for the operational week is derived from GPM IMERG early run data.

The cumulative precipitation over the week before the latest Landsat overpass is considered as **nowcast precipitation**. Additionally, the cumulative forecasted precipitation for the next week after the Landsat overpass is included.

## Net Water Requirement
The net water requirement is calculated as:

$$
\text{Net Water Requirement} = \text{SEBAL ET} - \text{PET} - \text{Per} + P_{nc} + P_{fc}
$$

Where:  

- **SEBAL ET**: Estimated actual evapotranspiration  
- **PET**: Potential evapotranspiration (Penman-Monteith method)  
- **Per**: Percolation losses  
- \(P_{nc}\): Nowcast precipitation   
- \(P_{fc}\): Forecasted precipitation   



**Interpretation:**

1. If net water requirement is **positive**, the water supply (including precipitation and percolation) meets or exceeds crop demand, so no additional irrigation is needed.

2. If net water requirement is **negative**, the water supply is insufficient, indicating additional irrigation is required for the upcoming week.

