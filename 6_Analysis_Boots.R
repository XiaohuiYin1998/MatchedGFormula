library(dplyr)
library(ggplot2)

n     <- as.integer(Sys.getenv('GFORM_N',     unset = '30000'))
setup <- as.integer(Sys.getenv('GFORM_SETUP', unset = '5'))

results_dir <- 'Boots_Results'
rst      <- read.table(sprintf('%s/n_%d_setup%d.txt',      results_dir, n, setup),
                       header = TRUE)
rst_true <- read.table(sprintf('%s/n_%d_setup%d_true.txt', results_dir, n, setup),
                       header = TRUE)

rst <- rst %>%
  left_join(rst_true, by = c('Time', 'treatment'))

# Per-simulation summary: bootstrap mean / sd / CI coverage for each estimator.
summary_each_simu <- rst %>%
  group_by(treatment, iterative, matched, Time, dataID) %>%
  summarise(Boots_mean     = mean(Risk),
            Boots_sd       = sd(Risk),
            Quantile_lower = quantile(Risk, 0.05),
            Quantile_upper = quantile(Risk, 0.95),
            Quantile_cover = I(Quantile_upper > TrueRisk[1] & Quantile_lower < TrueRisk[1]),
            HPD_lower      = coda::HPDinterval(coda::mcmc(Risk), prob = 0.9)[1],
            HPD_upper      = coda::HPDinterval(coda::mcmc(Risk), prob = 0.9)[2],
            HPD_cover      = I(HPD_upper > TrueRisk[1] & HPD_lower < TrueRisk[1]),
            Aysm_lower     = Boots_mean - qnorm(0.95) * Boots_sd,
            Aysm_upper     = Boots_mean + qnorm(0.95) * Boots_sd,
            Aysm_cover     = I(Aysm_upper > TrueRisk[1] & Aysm_lower < TrueRisk[1]),
            .groups        = 'drop')

# Aggregate across the 100 simulated datasets.
summary_all <- summary_each_simu %>%
  group_by(treatment, iterative, matched, Time) %>%
  summarise(mean_Boots_mean     = mean(Boots_mean),
            mean_Boots_sd       = mean(Boots_sd),
            sd_Boots_mean       = sd(Boots_mean),
            Quantile_cover_rate = mean(Quantile_cover),
            HPD_cover_rate      = mean(HPD_cover),
            Aysm_cover_rate     = mean(Aysm_cover),
            .groups             = 'drop')

print(summary_all)

p1 <- summary_each_simu %>%
  ggplot(aes(x = factor(Time), y = Boots_mean, fill = factor(treatment))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(rst_true,
             mapping = aes(x = factor(Time), y = TrueRisk),
             position = position_dodge(width = 0.8),
             color = 'red', shape = 8, size = 1) +
  facet_grid(
    rows = vars(iterative),
    cols = vars(matched),
    labeller = labeller(
      iterative = c(`1` = 'Iterative', `0` = 'Non-iterative'),
      matched   = c(`1` = 'Matched',   `0` = 'Complete')
    )
  ) +
  labs(x = 'Time', y = 'Boots_mean', fill = 'Treatment') +
  theme_minimal() +
  theme(legend.position = 'bottom')

print(p1)

p2 <- summary_each_simu %>%
  ggplot(aes(x = factor(Time), y = Boots_mean, fill = factor(treatment))) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_point(
    data = rst_true,
    mapping = aes(x = factor(Time), y = TrueRisk),
    position = position_dodge(width = 0.8),
    color = 'red', shape = 8, size = 1
  ) +
  facet_grid(
    rows = vars(treatment),
    cols = vars(iterative, matched),
    scales = 'free_y',
    labeller = labeller(
      iterative = c(`1` = 'Iterative', `0` = 'Non-iterative'),
      matched   = c(`1` = '(matched)', `0` = '(complete)')
    )
  ) +
  labs(x = 'Time', y = 'Boots_mean', fill = 'Treatment') +
  theme_minimal() +
  theme(legend.position = 'bottom')

print(p2)
