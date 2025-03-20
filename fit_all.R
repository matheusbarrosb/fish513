# Load packages ----------------------------------------------------------------
package_list = c("dplyr", "cmdstanr", "ggplot2", "here", "tidyr", "ggpubr")
lapply(package_list, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

# Load functions ---------------------------------------------------------------
function_directory = file.path(here::here(), "R", "Functions/")
function_files     = list.files(function_directory)
for (i in 1:length(function_files)) source(paste0(function_directory, function_files[i]))

# Data processing --------------------------------------------------------------

# 1. UNEXPLOITED SPECIES (basically Myctobase)
unexp_data = read.csv(here::here("Data/Myctobase/myctobase_clean.csv"))
#View(unexp_data)

# Sum across all species for single trend
constant   = 1e-10 # avoids NaNs
unexp_data = unexp_data %>%
  group_by(year) %>%
  summarise(obs = sum(-log(n_m3+constant), na.rm = TRUE),
            n)

# fix random NA
unexp_data[length(unexp_data$year),1] = 2018
  
# unexp_data %>%  
# ggplot(aes(x = year, y = obs/mean(obs))) +
#   geom_line() +
#   geom_point() +
#   custom_theme() +
#   ylim(0, 4)


# 2. EXPLOITED SPECIES (RAM + PELGAS)
taxonomy = read.csv(here::here("Data/RAM/taxonomy.csv"))
ts       = read.csv(here::here("Data/RAM/ts.csv"))
stock    = read.csv(here::here("Data/RAM/stock.csv"))

stock2 = data.frame(stock$stockid, stock$scientificname, stock$region)
names(stock2) = c("stockid", "scientificname", "region")

# Get list of forage fish from taxonomy database 
taxonomy = taxonomy %>% filter(FisheryType == "Forage Fish") # 44 forage fish stocks

# Filter time series data
exp_data = ts %>%
  left_join(stock2, by = "stockid") %>%
  filter(scientificname %in% taxonomy$scientificname)

# # join PELGAS data
# pelgas = read.csv(here::here("Data/PELGAS/PELGAS_clean.csv"))
# pelgas = pelgas %>% 
#   select(year, sciName, wbiom) %>%
#   rename(tsyear = year,
#          scientificname = sciName) 
# 
# exp_data = full_join(exp_data, pelgas, by = c("tsyear", "scientificname"))

exp_data = exp_data %>%
  filter(tsyear %in% c(1960:2015)) %>% # only including years
  na.exclude() %>%
  group_by(tsyear) %>%
  summarise(obs  = median(tsvalue, na.rm = TRUE),
            sd   = sd(tsvalue, na.rm = TRUE),
            se   = sd/sqrt(n())) 

exp_data %>%
  ggplot(aes(x = tsyear, y = obs)) +
  geom_line() +
 # geom_errorbar(aes(ymin = obs - se, ymax = obs + se), width = 0.2) +
  custom_theme()

names(exp_data) = c("year", "obs", "sd", "se")

# Preparing data for Stan ------------------------------------------------------

# 1. UNEXPLOITED SPECIES (basically Myctobase)
unexp_mean = mean(unexp_data$obs, na.rm = TRUE)
data_list_unexp = list(
  N            = length(unique(unexp_data$year)),           # 28 years of data
  S            = 1,                                         # 1 single time series  
  M            = 1,                                         # = S, could simplify model specification?
  y            = unexp_data$obs,                            # observations
  states       = 1:1,
  n_obsvar     = 1,
  proVariances = c(1,0),
  obsVariances = 1,
  trends       = c(0,0),
  est_trend    = 0,
  est_nu       = 0,                                         # = 1 turns on overdispersion
  family       = 1,                                         # = 1 for gaussian
  n_provar     = 1,
  n_trends     = 1,
  n_pos        = dim(unexp_data)[1],                        # = N if S = 1
  row_indx_pos = rep(1, dim(unexp_data)[1]),
  col_indx_pos = as.numeric(as.factor(unexp_data$year)),
  est_A        = c(0,0),
  n_A          = 0
)

# 2. EXPLOITED SPECIES (RAM + PELGAS)
exp_mean = mean(exp_data$obs, na.rm = TRUE)

data_list_exp = list(
  N            = length(unique(exp_data$year)),             # 56 years of data
  S            = 1,                                         # 1 single time series  
  M            = 1,                                        
  y            = exp_data$obs/exp_mean,                            # observations
  states       = 1:1,
  n_obsvar     = 1,
  proVariances = c(1,0),
  obsVariances = 1,
  trends       = c(0,0),
  est_trend    = 0,
  est_nu       = 1,                                         # = 1 turns on overdispersion
  family       = 1,                                         # = 1 for gaussian
  n_provar     = 1,
  n_trends     = 1,
  n_pos        = dim(exp_data)[1],                        # = N if S = 1
  row_indx_pos = rep(1, dim(exp_data)[1]),
  col_indx_pos = as.numeric(as.factor(exp_data$year)),
  est_A        = c(0,0),
  n_A          = 0
)

# Model fitting ----------------------------------------------------------------

# compile model
model_directory = file.path(here::here(), "MARSS", "MARSS.stan")
model_file      = cmdstan_model(model_directory)

# mcmc settings 
mcmc_list = list(n_mcmc      = 1000,
                 n_burn      = 200,
                 n_chain     = 3,
                 n_thin      = 1,
                 step_size   = 0.4,
                 adapt_delta = 0.9)

# 1. UNEXPLOITED SPECIES (basically Myctobase)
fit_unexp = model_file$sample(
  data            = data_list_unexp,
  seed            = 2025,
  chains          = mcmc_list$n_chain,
  parallel_chains = mcmc_list$n_chain,
  iter_warmup     = mcmc_list$n_burn,
  iter_sampling   = mcmc_list$n_mcmc,
  adapt_delta     = 0.97,
  step_size       = 0.05,
  refresh         = 100,
  max_treedepth   = 20
)

# 2. EXPLOITED SPECIES (RAM + PELGAS)
fit_exp = model_file$sample(
  data            = data_list_exp,
  seed            = 2025,
  chains          = mcmc_list$n_chain,
  parallel_chains = mcmc_list$n_chain,
  iter_warmup     = mcmc_list$n_burn,
  iter_sampling   = mcmc_list$n_mcmc,
  adapt_delta     = 0.97,
  step_size       = 0.05,
  refresh         = 100,
  max_treedepth   = 20
)

# Extract and plot predictions -------------------------------------------------
plot_directory = file.path(here::here(), "Plots", "Final")

# 1. UNEXPLOITED SPECIES (basically Myctobase)
preds_unexp = fit_unexp$summary(variables = "pred",  ~ quantile(.x, probs = c(0.5, 0.1, 0.9)))

preds_unexp$year = unexp_data$year
preds_unexp$obs  = data_list_unexp$y

names(preds_unexp) = c("par", "mean", "low", "upp", "year", "obs")

preds_unexp %>%
  ggplot() +
  geom_line(aes(year, y = mean, color = "Estimated trend"), size = 1) +
  geom_ribbon(aes(x = year, ymin = low, ymax = upp), alpha = 0.15) +
  geom_point(aes(x = year, y = obs)) +
  custom_theme() +
  theme(strip.text.x = element_text(face = "italic"),
        legend.position = "top",
        legend.title = element_blank()) +
  ylab(bquote("Count/m"^3)) +
  xlab("Year") +
 # ylim(0, 3) +
  ggtitle("Unexploited pelagics")

ggsave("non_exploited_trend.pdf", path = plot_directory, width = 5, height = 3)

plot_sigmas(fit = fit_unexp, legend_coord = c(0.85, 0.85))
ggsave("non_exploited_sigmas.pdf", path = plot_directory, width = 7.5, height = 2.6)


# 2. EXPLOITED SPECIES (RAM + PELGAS)
preds_exp = fit_exp$summary(variables = "pred", ~ quantile(.x, probs = c(0.5, 0.1, 0.9)))

preds_exp$year = exp_data$year
preds_exp$obs  = exp_data$obs/exp_mean

names(preds_exp) = c("par", "mean", "low", "upp", "year", "obs")

preds_exp %>%
  ggplot() +
  geom_line(aes(year, y = mean, color = "Estimated trend"), size = 1) +
  geom_ribbon(aes(x = year, ymin = low, ymax = upp), alpha = 0.15) +
  geom_point(aes(x = year, y = obs)) +
  custom_theme() +
  theme(strip.text.x = element_text(face = "italic"),
        legend.position = "top",
        legend.title = element_blank()) +
  ylim(0,2.5) +
  ggtitle("Exploited pelagics") +
  xlab("Year") +
  ylab("Scaled biomass")

ggsave("exploited_trend.pdf", path = plot_directory, width = 5, height = 3)

plot_sigmas(fit = fit_exp, legend_coord = c(0.25, 0.85))
ggsave("exploited_sigmas.pdf", path = plot_directory, width = 7.5, height = 2.6)


# Export model predictions as .csv files ---------------------------------------
export_path = file.path(here::here(), "Data", "MARSS_outputs/")
write.csv(preds_unexp, file = paste0(export_path, "predictions_unexploited.csv"))
write.csv(preds_exp,   file = paste0(export_path, "predictions_exploited.csv"))








