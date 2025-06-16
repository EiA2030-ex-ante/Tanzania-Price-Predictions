# Predictive mapping of wholesale grain prices for rural areas in Tanzania

This repository contains the data, data preparation steps, and code used for a study on spatio-temporal prediction of wholesale crop prices for rural areas in Tanzania. This project builds and evaluates Random Forest models to predict monthly crop prices at high spatial resolution. The whole workflow used is compiled in this document- https://eia2030-ex-ante.github.io/Tanzania-Price-Predictions/

# Background
Timely and location-specific data on wholesale crop prices in rural areas is often scarce in Tanzania. Available data is usually aggregated at the regional or national level (e.g., Dar es Salaam or zone-level averages), which fails to reflect the substantial price variability observed across rural markets. This lack of granular market information limits the ability of policymakers to respond effectively to local food price fluctuations, plan interventions, or invest in rural value chains.

To address this, the study develops and tests a spatially explicit machine learning framework to predict monthly prices for eight key staple crops across rural Tanzania. These crops include maize, rice, sorghum, bulrush millet, finger millet, wheat, beans, and potatoes. The models use Random Forest algorithms trained on observed prices from 44 markets, leveraging environmental, temporal, and market-access covariates. Key predictors include proximity to cities and ports, population density, travel time, and bioclimatic variables.

# Repository Structure
# 1 Data/
This folder contains the cleaned market price data used for model training and validation. The raw data was sourced from the Tanzania Ministry of Industry and Trade (MIT), and includes monthly wholesale prices for key staple crops reported at various markets across the country. The cleaned dataset covers a total of 44 markets and includes these crops; maize, rice, sorghum, bulrush millet, finger millet, wheat, beans, and potatoes.

This folder also includes the spatio-temporal raster outputs for all the crops (maize, rice, sorghum, bulrush millet, finger millet, wheat, beans, and potatoes) from the prediction models—price surfaces that estimate monthly crop prices across rural Tanzania at high resolution.

# 2. Code/
This folder contains the R Markdown (.Rmd) file that document the full workflow—from data cleaning and preprocessing to model fitting, validation, and prediction. The code includes Random Forest model training, spatial cross-validation routines (e.g., leave-N-markets-out), and the generation of final monthly price prediction maps.

The entire modeling process is fully reproducible, with results also compiled and published here: https://eia2030-ex-ante.github.io/Tanzania-Price-Predictions/.
