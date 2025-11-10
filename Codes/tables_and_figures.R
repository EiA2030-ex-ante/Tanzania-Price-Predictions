stats <- mypts %>%
  filter(pkg > 0) %>%  
  group_by(Crop) %>%
  summarise(
    Mean = mean(pkg, na.rm = TRUE),
    Median = median(pkg, na.rm = TRUE),
    Minimum = min(pkg, na.rm = TRUE),
    Maximum = max(pkg, na.rm = TRUE),
    Std_Dev = sd(pkg, na.rm = TRUE),
    IQR = IQR(pkg, na.rm = TRUE),
    Observations = n()
  ) %>%
  arrange(Crop)
print(stats)

library(officer)
library(flextable)

# Convert to flextable
ft <- flextable(stats)

# Format: bold headers, autofit column widths, gridlines
ft <- ft %>%
  bold(part = "header") %>%
  autofit() %>%
  theme_vanilla() %>%
  border_remove() %>%
  hline_top(part = "all", border = fp_border()) %>%
  hline_bottom(part = "all", border = fp_border()) 

# Save to Word
doc <- read_docx() %>%
  body_add_flextable(ft)

print(doc, target = "stats_table.docx")

#-----------------------------------------------------------------------
unique(mypts$Crop)
mypts <- as.data.frame(mypts)
mypts <- mypts %>%
  mutate(Crop = fct_recode(Crop,
                           "Maize"          = "Maize",
                           "Rice"           = "Rice",
                           "Sorghum"        = "Sorghum",
                           "B.Millet" = "B.Millet",
                           "F.Millet"  = "F.Millet",
                           "Wheat"          = "Wheat",
                           "Bean"          = "Beans",
                           "Potato"         = "Potato"
  ))

png("boxplot_pkg_distribution_per_crop.png", width = 800, height = 600, res = 125)
# Boxplot of pkg distribution per Crop
ggplot(mypts %>% filter(pkg > 0), aes(x = Crop, y = pkg, fill = Crop)) +
  geom_boxplot(outlier.colour = "black", outlier.alpha = 0.4) +
  labs(
    title = "",
    x = "",
    y = "Price (TZs/kg)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, size = 10),   # horizontal labels
    axis.text.y = element_text(size = 10),                        # standardised y-axis labels
    axis.title.y = element_text(size = 10),                       # standardised y-axis title
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"                                      # remove legend
  )
dev.off()

#----------------------------------------------------------------------
# Cumulative probability plots
# List of crops excluding maize
#setwd("H:/Tanzania Price data/Datasets/Pred-plots")

crops <- c("beans", "rice", "sorghum", "bmillet", "fmillet", "wheat", "potato")

# Store CDF data for plotting
cdf_list <- list()

for (crop in crops) {
  
  # Get predicted raster list for the crop
  pred_list <- get(paste0("predictions_", crop))  # e.g., predictions_beans
  
  # Prepare a data frame for all pixels and months
  df_crop <- data.frame()
  
  for (month in 1:12) {
    vals <- values(pred_list[[month]])
    vals <- vals[!is.na(vals)]
    
    # Compute cumulative probabilities
    vals_sorted <- sort(vals)
    probs <- seq_along(vals_sorted) / length(vals_sorted)
    
    df_month <- data.frame(
      price = vals_sorted,
      cum_prob = probs,
      month = month,
      crop = crop
    )
    
    df_crop <- rbind(df_crop, df_month)
  }
  
  cdf_list[[crop]] <- df_crop
}

# Combine all crops
cdf_data <- bind_rows(cdf_list)

# convert month to factor with labels
cdf_data$month <- factor(cdf_data$month, labels = month.name)

# Define nicer labels for crops
crop_labels <- c(
  beans = "Bean",
  bmillet = "Bulrush Millet",
  fmillet = "Finger Millet",
  potato = "Potato",
  rice = "Rice",
  sorghum = "Sorghum",
  wheat = "Wheat"
)

# png("cumulative_probability_plots.png", width = 2400, height = 2400, res = 300)
# 12 colors 
month_colors <- viridis::mako(12)

ggplot(cdf_data, aes(x = price, y = cum_prob, color = month)) +
  geom_line(size = 1) +
  facet_wrap(~crop, scales = "free", ncol = 2, 
             labeller = labeller(crop = crop_labels)) +   # Use full names
  scale_color_manual(values = month_colors) +
  labs(x = "Predicted price (TZS/kg)", y = "Cumulative probability", color = "") +
  theme_minimal(base_size = 10) +
  theme(
    legend.justification = c("right", "top"),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 0, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.key.width = unit(1.5, "lines")
  )
dev.off()

#----------------------------------------------------------------------
# Variable importance plots
library(dplyr)
library(ggplot2)

# Extract variable importance
imp <- importance(rf) %>% as.data.frame()
imp$Variable <- rownames(imp)

# Choose importance column
importance_col <- if ("IncNodePurity" %in% names(imp)) {
  "IncNodePurity"
} else if ("%IncNodePurity" %in% names(imp)) {
  "%IncNodePurity"
} else {
  "IncMSE"
}

# Create nicer labels for publication
pretty_labels <- c(
  maize       = "Maize",
  rice        = "Rice",
  sorghum     = "Sorghum",
  bmillet     = "Bulrush millet",
  fmillet     = "Finger millet",
  wheat       = "Wheat",
  beans       = "Bean",
  potato      = "Potato",
  Month       = "Month",
  Year        = "Year",
  ttport_1    = "Travel time to port",
  ttcity_u5   = "Travel time to cities",
  popdens     = "Population density",
  bio_3       = "bio3",
  bio_6       = "bio6",
  bio_9       = "bio9",
  bio_12      = "bio12",
  bio_18      = "bio18",
  rain.sum.lag = "Lagged rainfall"
)

# Arrange and keep all
imp_plot <- imp %>%
  arrange(desc(.data[[importance_col]])) %>%
  mutate(
    PrettyVar = recode(Variable, !!!pretty_labels),
    PrettyVar = factor(PrettyVar, levels = PrettyVar)
  ) %>%
  slice_head(n = 19)

# png("variable_importance_plot.png", width = 800, height = 600, res = 125)
# Plot
ggplot(imp_plot, aes(x = PrettyVar, y = .data[[importance_col]])) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    x = "Predictor",
    y = "Importance",
    title = ""
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
dev.off()

#----------------------------------------------------------------------
# Validation and robustness checks plots
library(tidyverse)
library(patchwork)

# -----------------------------
# 1. Read CSVs
# -----------------------------
comparison_df_wide3 <- read.csv("Train-Test-Validation-results.csv")
comparison_df_wide  <- read.csv("Temporal-split-results.csv")
comparison_df_wide_ar <- read.csv("Temporal-split-results_AR.csv")

# -----------------------------
# 2. Add Validation column & rename crops
# -----------------------------
df1 <- comparison_df_wide3 %>%
  mutate(
    Validation = "Train–Test Split",
    Crop = recode(Crop, "Beans" = "Bean")
  )

df2 <- comparison_df_wide %>%
  mutate(
    Validation = "Temporal Split (Standard RF)",
    Crop = recode(Crop, "Beans" = "Bean")
  )

df3 <- comparison_df_wide_ar %>%
  mutate(
    Validation = "Temporal Split (Autoregressive RF)",
    Crop = recode(Crop, "Beans" = "Bean")
  )

# -----------------------------
# 3. Combine all datasets
# -----------------------------
all_df <- bind_rows(df1, df2, df3)

# Optional: trim whitespace from column names
names(all_df) <- trimws(names(all_df))
all_df$Crop <- trimws(all_df$Crop)

write.csv(all_df, "all_validation_results.csv", row.names = FALSE)

# -----------------------------
# 4. Reshape RMSE for plotting
# -----------------------------
rmse_long <- all_df %>%
  dplyr::select(Crop, Validation, RMSE_Pooled, RMSE_Crop_Specific) %>%
  pivot_longer(
    cols = starts_with("RMSE"),
    names_to = "Model",
    values_to = "RMSE"
  ) %>%
  mutate(Model = recode(Model,
                        "RMSE_Pooled" = "Pooled",
                        "RMSE_Crop_Specific" = "Crop-Specific"))

# -----------------------------
# 5. Reshape R² for plotting
# -----------------------------
r2_long <- all_df %>%
  dplyr::select(Crop, Validation, R_squared_Pooled, R_squared_Crop_Specific) %>%
  pivot_longer(
    cols = starts_with("R_squared"),
    names_to = "Model",
    values_to = "R2"
  ) %>%
  mutate(Model = recode(Model,
                        "R_squared_Pooled" = "Pooled",
                        "R_squared_Crop_Specific" = "Crop-Specific"))

# Define the desired order
crop_order <- c("Maize", "Rice", "Sorghum", "B.Millet", "F.Millet", "Wheat", "Bean", "Potato")

# Apply to both long datasets
rmse_long <- rmse_long %>%
  mutate(Crop = factor(Crop, levels = crop_order))

r2_long <- r2_long %>%
  mutate(Crop = factor(Crop, levels = crop_order))

# -----------------------------
# 6. Plot RMSE (points + lines instead of bars)
# -----------------------------
p1 <- ggplot(rmse_long, aes(x = Crop, y = RMSE, color = Model, group = Model)) +
  geom_point(size = 1.5) +
  geom_line(size = 0.7) +
  facet_wrap(~Validation, ncol = 1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(y = "RMSE (TZS/kg)", x = "") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    plot.margin = margin(5, 5, 5, 5),
    clip = "on"
  )


# -----------------------------
# 7. Plot R²
# -----------------------------
p2 <- ggplot(r2_long, aes(x = Crop, y = R2, color = Model, group = Model)) +
  geom_point(size = 1.5) +
  geom_line(size = 0.7) +
  facet_wrap(~Validation, ncol =1, scales = "free_y") +
  theme_minimal(base_size = 14) +
  labs(y = expression(R^2), x = "") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2),
    plot.margin = margin(5, 5, 5, 5),
    clip = "on"
  )

# -----------------------------
# 8. Combine & save
# -----------------------------
png("validation_plots.png", width = 2400, height = 2400, res = 300)
combined_plot <- p1 | p2
print(combined_plot)
dev.off()


#----------------------------------------------------------------------
# Load libraries
library(ggplot2)

# data
market_holdout_df <- read.csv("market_holdout_df.csv")

lm_fit <- lm(Avg_R2 ~ Markets_Heldout, data = market_holdout_df)

# Fit segmented (piecewise) model with breakpoint estimation
seg_fit <- segmented(lm_fit, seg.Z = ~Markets_Heldout)

# Add predicted values to data
market_holdout_df$seg_fit <- predict(seg_fit)

# Plot
# png("market_holdout_plot.png", width = 800, height = 600, res = 125)
ggplot(market_holdout_df, aes(x = Markets_Heldout, y = Avg_R2)) +
  geom_ribbon(aes(ymin = Avg_R2 - SD_R2, ymax = Avg_R2 + SD_R2), fill = "grey", alpha = 0.4) +
  geom_line(aes(y = seg_fit), color = "black", size = 1) +
  geom_point(color = "black") +
  labs(title = "",
       x = "Number of Markets Held Out",
       y = "Mean R²") +
  theme_minimal() +
  theme(
    #axis.line = element_line(color = "black", size = 0.6),  
    #axis.line.x = element_line(color = "black", size = 0.6),  
    #axis.line.y = element_line(color = "black", size = 0.6), 
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 10),  
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
  )
dev.off()

#----------------------------------------------------------------------
#-- Boxplot of Regional Avg vs Predicted----------------------------------
#We need the box plot for all months, or maybe violin plots?
# --- Expand analysis across multiple months instead of only August ---

# loop across months of 2024, extract comparison_df for each month
all_months <- 1:12
comparison_list <- list()

for (m in all_months) {
  # Load monthly prediction raster 
  pred_file <- sprintf("H:/Tanzania Price data/Datasets/Pred-plots/Maize/maize_price_rf_pred_%02d.tif", m)
  if (!file.exists(pred_file)) next
  
  pred_maize <- rast(pred_file) |> project(crs(tza1))
  names(pred_maize) <- "Predicted"
  
  # Build comparison_df for this month
  pred_df <- as.data.frame(pred_maize, xy = TRUE, na.rm = TRUE)
  nat_df  <- as.data.frame(rnatavg,   xy = TRUE, na.rm = TRUE)
  reg_df  <- as.data.frame(regional_raster, xy = TRUE, na.rm = TRUE)
  
  tmp_df <- pred_df %>%
    left_join(nat_df %>% rename(National_Avg = national_avg), by = c("x", "y")) %>%
    left_join(reg_df %>% rename(Regional_Avg = Reg_max_july_pkg), by = c("x", "y")) %>%
    mutate(
      Pred_vs_Nat = National_Avg - Predicted,
      Pred_vs_Reg = Regional_Avg - Predicted,
      Month = m,
      Year = 2024
    )
  
  # Attach region names
  pts <- vect(tmp_df, geom = c("x", "y"), crs = crs(tza1))
  tmp_df$Region <- terra::extract(tza1["NAME_1"], pts)[,2]
  
  comparison_list[[as.character(m)]] <- tmp_df
}

# Combine all months into one dataframe
comparison_all <- bind_rows(comparison_list) %>%
  filter(!is.na(Region))
head(comparison_all, 5)

# --- Faceted boxplot/violin plot ---

# Get unique months from comparison_all
unique_months <- sort(unique(comparison_all$Month))

# Loop over months and save plots
for (m in unique_months) {
  df_month <- subset(comparison_all, Month == m)
  
  # Keep only regions with at least one non-missing value
  df_month <- df_month %>%
    group_by(Region) %>%
    filter(!all(is.na(Pred_vs_Reg))) %>%
    ungroup()
  
  png(paste0("predicted_vs_regional_prices_", m, ".png"), width = 800, height = 600, res = 125)
  
  print(
    ggplot(df_month, aes(x = Region, y = Pred_vs_Reg)) +
      geom_boxplot(fill = "#A6CEE3", alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
      labs(
        x = "",
        y = "Regional Price - Predicted Price (TZS/kg)",
        title = ""
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text = element_text(size = 10),  
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2))
  )
  
  dev.off()
}

#-------------------------------------------------------------------------------------
# Convert months to words
comparison_all$Month <- factor(
  comparison_all$Month,
  levels = 1:12,
  labels = month.name
)

png("Regional_box_plots.png", width = 2400, height = 2400, res = 200)
ggplot(comparison_all, aes(x = Region, y = Pred_vs_Reg)) +
  geom_boxplot(fill = "#A6CEE3", alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    x = "",
    y = "Regional Price - Predicted Price (TZS/kg)",
    title = ""
  ) +
  theme_minimal(base_size = 10) +  # bump up base size
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 10),  # facet labels
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  ) +
  facet_wrap(~ Month, ncol = 3)   # 3x4 grid

dev.off()

#-------------------------------------------------------------------------------------
# Plotting difference maps for all months national vs predicted
library(modelbased)
library(parameters)
library(ggplot2)
library(RColorBrewer)
library(effsize)
library(tidyverse)
library(gghalves)
library(viridis)
library(terra)

# Tanzania extent and base raster
tza_extent <- ext(tza1) |> floor()
base_rast <- crop(rast(res=1/12), tza_extent)
base_rast <- project(base_rast, crs(tza1))

# Convert prices to spatial points
mypts <- vect(prices.monthly.long, geom = c("Longitude", "Latitude"), crs = crs(tza0), keepgeom = TRUE)
district_info <- terra::extract(tza2, mypts)
prices.monthly.long$Region_GADM <- district_info$NAME_1
prices.monthly.long$District_GADM <- district_info$NAME_2

# Fix Songwe mismatch
prices.monthly.long[Region_GADM == "Mbeya" & District_GADM == "Mbeya Rural", 
                    `:=`(Region_GADM = "Songwe", District_GADM = "Songwe")]

# Loop through all months
# Loop through all months
for (m in 1:12) {
  
  message("Processing month: ", m)
  
  # Filter points for this month and crop
  month_pts <- as.data.frame(mypts) %>%
    filter(Crop == "Maize", Year == 2024, Month == m)
  
  if(nrow(month_pts) == 0) next
  
  # --- National maximum (Dar es Salaam) ---
  nat_max <- month_pts %>%
    filter(Region == "Dar es Salaam") %>%
    summarise(Max_pkg = max(pkg, na.rm = TRUE)) %>%
    pull(Max_pkg)
  
  if(is.na(nat_max)) nat_max <- 0
  
  rnat <- base_rast
  values(rnat) <- nat_max
  rnatavg <- tryCatch(terra::mask(rnat, tza1), error = function(e) rnat)
  names(rnatavg) <- "national_max"
  
  # --- Regional maximum ---
  reg_max <- month_pts %>%
    group_by(Region_GADM) %>%
    summarise(Reg_max_pkg = max(pkg, na.rm = TRUE)) %>%
    ungroup()
  
  missing_regions <- setdiff(tza1$NAME_1, reg_max$Region_GADM)
  if(length(missing_regions) > 0){
    reg_max <- bind_rows(reg_max, tibble(Region_GADM = missing_regions, Reg_max_pkg = NA))
  }
  
  if(!"Njombe" %in% reg_max$Region_GADM){
    iringa_price <- reg_max %>% filter(Region_GADM == "Iringa")
    njombe_price <- iringa_price
    njombe_price$Region_GADM <- "Njombe"
    reg_max <- bind_rows(reg_max, njombe_price)
  }
  
  tza1_reg <- merge(tza1, reg_max, by.x = "NAME_1", by.y = "Region_GADM", all.x = TRUE)
  regional_raster <- rasterize(tza1_reg, base_rast, field = "Reg_max_pkg")
  if(all(is.na(values(regional_raster)))) {
    regional_raster <- base_rast
    values(regional_raster) <- NA
  }
  regional_raster <- tryCatch(resample(regional_raster, base_rast), error = function(e) regional_raster)
  names(regional_raster) <- "regional_max"
  
  # --- Predicted raster ---
  pred_file <- sprintf("H:/Tanzania Price data/Datasets/Pred-plots/Maize/maize_price_rf_pred_%02d.tif", m)
  if(!file.exists(pred_file)) next
  pred_maize <- rast(pred_file) |> project(crs(tza1))
  names(pred_maize) <- "pred_maize"
  
  # --- Extract values safely ---
  pred_values <- values(pred_maize)
  nat_values <- values(rnatavg)
  reg_values <- values(regional_raster)
  
  valid_idx <- which(!is.na(pred_values) & !is.na(nat_values) & !is.na(reg_values))
  if(length(valid_idx) == 0) next
  
  comparison_df <- data.frame(
    Predicted = pred_values[valid_idx],
    National_Avg = nat_values[valid_idx],
    Regional_Avg = reg_values[valid_idx]
  ) %>%
    mutate(
      Pred_vs_Nat = National_Avg - Predicted,
      Pred_vs_Reg = Regional_Avg - Predicted,
      Month = m
    )
  
  print(head(comparison_df))
  message("Month ", m, ": % Predicted > National = ", mean(comparison_df$Pred_vs_Nat > 0))
  message("Month ", m, ": % Predicted > Regional = ", mean(comparison_df$Pred_vs_Reg > 0))
}

#-------------------------------------------------------------------------------------
# Using Tmap for raster plotting
library(tmap)

# Store difference rasters
diff_rasters <- list()

for (m in 1:12) {
  
  # Load predicted raster
  pred_file <- sprintf("H:/Tanzania Price data/Datasets/Pred-plots/Maize/maize_price_rf_pred_%02d.tif", m)
  if(!file.exists(pred_file)) next
  pred_maize <- rast(pred_file) |> project(crs(tza1))
  names(pred_maize) <- "pred_maize"
  
  # Get Dar es Salaam max price for this month
  month_pts <- as.data.frame(mypts) %>% filter(Crop=="Maize", Year==2024, Month==m)
  nat_max <- month_pts %>% filter(Region=="Dar es Salaam") %>% summarise(Max_pkg = max(pkg, na.rm=TRUE)) %>% pull(Max_pkg)
  if(is.na(nat_max)) nat_max <- 0
  
  # Create national max raster
  rnat <- base_rast
  values(rnat) <- nat_max
  rnatavg <- tryCatch(terra::mask(rnat, tza1), error = function(e) rnat)
  
  # Difference raster
  diff_nat_pred <- rnatavg - pred_maize
  names(diff_nat_pred) <- paste0("Month_", m)
  
  diff_rasters[[m]] <- diff_nat_pred
}

# Stack all months
diff_stack <- rast(diff_rasters)

# Plot using tmap
tm_shape(diff_stack) +
  tm_raster(palette = "-RdYlBu", n = 7, title = "National - Predicted (TZS/kg)") +
  tm_facets(ncol = 4) +  # 4 columns, 3 rows
  tm_layout(main.title = "",
            legend.outside = TRUE)

#-------------------------------------------------------------------------------------
# Using ggplot for raster plotting
# Create a list to store rasters for each month
diff_rasters <- list()

for (m in 1:12) {
  # Load predicted raster
  pred_file <- sprintf("H:/Tanzania Price data/Datasets/Pred-plots/Maize/maize_price_rf_pred_%02d.tif", m)
  if(!file.exists(pred_file)) next
  pred_maize <- rast(pred_file) |> project(crs(tza1))
  
  # National max raster (Dar es Salaam)
  month_pts <- as.data.frame(mypts) %>% filter(Crop == "Maize", Year == 2024, Month == m)
  nat_max <- month_pts %>% filter(Region == "Dar es Salaam") %>% summarise(Max_pkg = max(pkg, na.rm=TRUE)) %>% pull(Max_pkg)
  if(is.na(nat_max)) nat_max <- 0
  rnat <- base_rast
  values(rnat) <- nat_max
  rnat <- tryCatch(terra::mask(rnat, tza1), error = function(e) rnat)
  
  # Difference raster
  diff_rasters[[m]] <- rnat - pred_maize
}

# Convert rasters to data frames for ggplot
plot_df <- do.call(rbind, lapply(1:12, function(m){
  r <- diff_rasters[[m]]
  df <- as.data.frame(r, xy=TRUE)
  df$Month <- m
  colnames(df)[3] <- "Diff"
  df
}))

# Convert Tanzania boundaries to sf
tza1_sf <- st_as_sf(tza1)

# Convert sf to data.frame suitable for ggplot
library(dplyr)
library(sf)

tza1_df <- tza1_sf %>%
  st_cast("MULTIPOLYGON") %>%  # ensure multipolygons are explicit
  st_cast("POLYGON") %>%       # break into individual polygons
  st_cast("LINESTRING") %>%    # convert boundaries to lines
  mutate(poly_id = row_number()) %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(X = X, Y = Y, L1 = L1) %>%
  group_by(L1) %>%
  mutate(id = L1)


# Replace numeric months with month names
plot_df$Month <- factor(plot_df$Month, levels = 1:12,
                        labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

png("National_vs_Predicted_diff.png", width = 2400, height = 1800, res = 200)
# Then plot as before
ggplot(plot_df) +
  geom_raster(aes(x = x, y = y, fill = Diff)) +
  geom_polygon(data = tza1_df, aes(x = X, y = Y, group = id),
               fill = NA, color = "gray30", size = 0.5) +
  scale_fill_gradient2(
    low = "darkred", mid = "white", high = "darkgreen", midpoint = 0,
    name = ""
  ) +
  facet_wrap(~Month, ncol = 4) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.2)
  )
dev.off()
#-------------------------------------------------------------------------------------


