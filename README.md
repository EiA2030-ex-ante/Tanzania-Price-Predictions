# Predictive mapping of wholesale grain prices for rural areas in Tanzania

This repository contains the data, data preparation steps, and code used for a study on spatio-temporal prediction of wholesale crop prices for rural areas in Tanzania. This project builds and evaluates Random Forest models to predict monthly crop prices at high spatial resolution. The whole workflow used is compiled in this document- https://eia2030-ex-ante.github.io/Tanzania-Price-Predictions/

# Background
Timely and location-specific data on wholesale crop prices in rural areas is often scarce in Tanzania. Available data is usually aggregated at the regional or national level (e.g., Dar es Salaam or zone-level averages), which fails to reflect the substantial price variability observed across rural markets. This lack of granular market information limits the ability of policymakers to respond effectively to local food price fluctuations, plan interventions, or invest in rural value chains.

To address this, the study develops and tests a spatially explicit machine learning framework to predict monthly prices for eight key staple crops across rural Tanzania. These crops include maize, rice, sorghum, bulrush millet, finger millet, wheat, beans, and potatoes. The models use Random Forest algorithms trained on observed prices from 44 markets, leveraging environmental, temporal, and market-access covariates. Key predictors include proximity to cities and ports, population density, travel time, and bioclimatic variables.

# Repository Structure
