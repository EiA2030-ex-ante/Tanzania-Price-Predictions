library(randomForest)
library(tidyverse)

set.seed(1983)
# Unique markets
unique_markets <- unique(mypts$Market)
total_markets <- length(unique_markets)

# Store results
market_holdout_results <- list()

# Loop over number of markets to hold out
for (num_holdout in 1:total_markets) {
  r2_vals <- c()
  
  for (rep in 1:5) {  # the held out markets are repeated five times to  average out randomness in which markets are held out.
    heldout_markets <- sample(unique_markets, num_holdout)
    
    train_subset <- mypts %>% filter(!Market %in% heldout_markets)
    test_subset  <- mypts %>% filter(Market %in% heldout_markets)
    
    model_rf <- randomForest(pkg ~ maize + rice + sorghum + bmillet + fmillet + wheat + beans + 
                               Month + Year + 
                               latitude + longitude +
                               ttport_1 + ttcity_u5 + popdens + 
                               bio_3 + bio_6 + bio_9 + bio_12 + bio_18 +
                               rain.sum.lag,
                             data = train_subset,
                             ntree = 500,
                             importance = TRUE,
                             na.action = na.omit)
    
    preds <- predict(model_rf, newdata = test_subset)
    actual_vals <- test_subset$pkg
    r2_score <- summary(lm(actual_vals ~ preds))$r.squared
    r2_vals <- c(r2_vals, r2_score)
  }
  
  market_holdout_results[[num_holdout]] <- data.frame(
    Markets_Heldout = num_holdout,
    Avg_R2 = mean(r2_vals),
    SD_R2 = sd(r2_vals)
  )
}

# Combine into a data frame
market_holdout_df <- bind_rows(market_holdout_results)
print(market_holdout_df)

# plot
ggplot(market_holdout_df, aes(x = Markets_Heldout, y = Avg_R2)) +
  #geom_line(color = "blue", size = 1) +
  geom_smooth(method = "loess", se = FALSE, color = "blue", linetype = "dashed") +
  geom_point(color = "blue") +
  #geom_errorbar(aes(ymin = Avg_R2 - SD_R2, ymax = Avg_R2 + SD_R2), 
  #              width = 0.3, color = "gray40") +
  labs(title = "Impact of Market Holdout on Predictive R²",
       x = "Number of Markets Held Out",
       y = "Mean R²") +
  theme_minimal()

