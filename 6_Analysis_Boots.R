library(dplyr)
library(ggplot2)
n <- 30000; setup <- 5
# source('0_regularizedParameterSetup.R')
# source('1_regularizedTrueACE.R')
# df_true <- rbind(data.frame(Time = 1:6, TrueRisk = true_risk0, treatment = 0),
#                  data.frame(Time = 1:6, TrueRisk = true_risk1, treatment = 1))
# write.table(df_true, file = sprintf('Boots_Results/n_%d_setup%d_true.txt', n, setup),
#             row.names = FALSE, col.names = TRUE,
#             quote = FALSE, sep = '\t')
rst <- read.table(sprintf('Boots_Results/n_%d_setup%d.txt', n, setup), header = TRUE)
rst_true <- read.table(sprintf('Boots_Results/n_%d_setup%d_true.txt', n, setup), header = TRUE)

rst <- rst %>%
  left_join(rst_true, by = c('Time', 'treatment'))

# rst %>% select(treatment, iterative, matched, Computation_Fit, Computation_Data) %>%
#   distinct() %>%
#   group_by(treatment, iterative, matched) %>%
#   mutate(Computation_Data = ifelse(is.na(Computation_Data), 0 ,Computation_Data)) %>%
#   summarise(mean(Computation_Fit + Computation_Data, na.rm=TRUE)) %>% View()

### 
### 1) For each dataID, mean/sd of 100 bootstrap sample (Boots_mean, Boots_sd);
### 2) mean(Boots_mean), sd(Boots_mean), mean(Boots_sd);
### 3) For each dataID, confidence interval (percentiles/HPD/Aysmp) based on 100 boots samples; 
### coverage of true
## 90% confidence interval of every risk
summary_each_simu <- rst %>% group_by(treatment, iterative, matched, Time, dataID) %>%
  summarise(Boots_mean = mean(Risk), 
            Boots_sd = sd(Risk),
            Quantile_lower = quantile(Risk, 0.05),
            Quantile_upper = quantile(Risk, 0.95),
            Quantile_cover = I(Quantile_upper > TrueRisk[1] & Quantile_lower < TrueRisk[1]),
            HPD_lower = coda::HPDinterval(coda::mcmc(Risk), prob = 0.9)[1],
            HPD_upper = coda::HPDinterval(coda::mcmc(Risk), prob = 0.9)[2],
            HPD_cover = I(HPD_upper > TrueRisk[1] & HPD_lower < TrueRisk[1]),
            Aysm_lower = Boots_mean - qnorm(0.95) * Boots_sd,
            Aysm_upper = Boots_mean + qnorm(0.95) * Boots_sd, 
            Aysm_cover = I(Aysm_upper > TrueRisk[1] & Aysm_lower < TrueRisk[1])) %>%
  ungroup()

### Summarize the 100 statistics
summary_all <- summary_each_simu %>% 
  group_by(treatment, iterative, matched, Time) %>%
  summarise(mean_Boots_mean = mean(Boots_mean),
            mean_Boots_sd = mean(Boots_sd),
            sd_Boots_mean = sd(Boots_mean),
            Quantile_cover_rate = mean(Quantile_cover),
            HPD_cover_rate = mean(HPD_cover),
            Aysm_cover_rate = mean(Aysm_cover))
  
  
summary_each_simu %>% 
  ggplot(mapping = aes(x = factor(Time), y = Boots_mean, fill = factor(treatment))) + 
  geom_boxplot(position = position_dodge(width = 0.8)) +
  facet_grid(
    rows = vars(iterative),
    cols = vars(matched),
    labeller = labeller(
      iterative = c(`0` = 'Non-iterative', `1` = 'Iterative'),
      matched = c(`0` = 'Complete', `1` = 'Match')
    )
  ) + 
  labs(x = 'Time', y = 'Boots_mean', fill = 'Treatment') +
  theme_minimal() + 
  theme(legend.position = 'bottom')

  