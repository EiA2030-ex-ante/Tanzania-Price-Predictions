# Predictive mapping of wholesale grain prices for rural areas in Tanzania

This repository contains the data, data preparation steps, and code used for a study on spatio-temporal prediction of wholesale crop prices for rural areas in Tanzania. This project builds and evaluates Random Forest models to predict monthly crop prices at high spatial resolution. The whole workflow used is compiled in this document- https://eia2030-ex-ante.github.io/Tanzania-Price-Predictions/

# Background
Timely and location-specific data on wholesale crop prices in rural areas is often scarce in Tanzania. Available data is usually aggregated at the regional or national level (e.g., Dar es Salaam or zone-level averages), which fails to reflect the substantial price variability observed across rural markets. This lack of granular market information limits the ability of policymakers to respond effectively to local food price fluctuations, plan interventions, or invest in rural value chains.

To address this, the study develops and tests a spatially explicit machine learning framework to predict monthly prices for eight key staple crops across rural Tanzania. These crops include maize, rice, sorghum, bulrush millet, finger millet, wheat, beans, and potatoes. The models use Random Forest algorithms trained on observed prices from 44 markets, leveraging environmental, temporal, and market-access covariates. Key predictors include proximity to cities and ports, population density, travel time, and bioclimatic variables.

# Repository Structure
# 1 Data
This folder contains the cleaned market price data used for model training and validation. The raw data was sourced from the Tanzania Ministry of Industry and Trade (MIT), and includes monthly wholesale prices for key staple crops reported at various markets across the country. The cleaned dataset—https://raw.githubusercontent.com/EiA2030-ex-ante/Tanzania-Price-Predictions/main/Data/Tanzania_Price_Data_AllCrops_with_Coordinates4.csv — covers 44 markets and includes these crops: maize, rice, sorghum, bulrush millet, finger millet, wheat, beans, and potatoes.

This folder also includes the spatio-temporal raster outputs for all the crops listed above, generated from the prediction models. These outputs provide monthly price surface maps that estimate wholesale crop prices at high spatial resolution across rural Tanzania.

# 2. Code
This folder contains the R Markdown (.Rmd) file that document the full workflow—from data cleaning and preprocessing to model fitting, validation, and prediction. The code includes Random Forest model training, spatial cross-validation routines (e.g., leave-N-markets-out), and the generation of final monthly price prediction maps.

The entire modeling process is fully reproducible, with results also compiled and published here: https://eia2030-ex-ante.github.io/Tanzania-Price-Predictions/.

# License
This project is licensed under the Creative Commons Attribution 4.0 International License. See the LICENSE file for details.

# Acknowledgments
This work was made possible through the OneCGIAR Initiative on Excellence in Agronomy (INV-005431), as well as through the Guiding Acid Soil Management Investments in Africa (GAIA) project (Grant no: INV-029117), supported by the Bill & Melinda Gates Foundation (BMGF). We would like to thank all funders supporting research through contributions to the CGIAR Trust Fund: https://www.cgiar.org/funders/.


