library(cmdstanr)
library(tidyverse)

species_i_pre_df = readRDS("species_i_pre_df_EXAMPLE.RDS")
#=== Massage data
species_i_df <- species_i_pre_df %>% ungroup() %>% 
  dplyr::select(ID,Year,Site, log.spawner) %>%  # get just the columns that I need
  # drop_na(log.spawner) %>%
  group_by(ID,Site) %>% 
  complete(Year = min(Year):max(Year),
           fill = list(log.spawner = NA)) %>%
  ungroup() %>%
  pivot_wider(names_from = c("Site","ID"), values_from = "log.spawner") %>% 
  arrange(Year) %>% 
  column_to_rownames(var = "Year") %>% # make the years rownames
  # filter(if_any(everything(), ~ !is.na(.)))%>%
  as.matrix() %>% # turn into a matrix with year down the rows
  t() # make time across the columns

#== Set up for Stan
site_vec <- species_i_pre_df %>%  distinct(Site)%>% pull(Site) 
id_vec <- species_i_pre_df %>% ungroup() %>% distinct(ID,Indiv.pop.map) %>% pull(Indiv.pop.map) 

n_states <- max(as.numeric(forcats::fct_inorder(id_vec))) # site_vec
r_unequal <- seq(1, nrow(species_i_df))
r_equal <- rep(1, nrow(species_i_df))
uq_unequal <- seq(1, n_states)
uq_equal <- rep(1, n_states)

#= MARSS specifications
marss_list = list(states = as.numeric(fct_inorder(id_vec)),# as.numeric(forcats::fct_inorder(id_vec)), #site_vec
                  obsVariances = r_equal, #as.numeric(factor(id_vec)), 
                  proVariances = uq_equal,
                  trends = uq_equal,
                  stocks = uq_unequal) 

# set up data_list
source("R/marss_stan_functions.R")
data_list_tmp <- setup_data(y = species_i_df,
                            est_nu = FALSE,
                            est_trend = FALSE,
                            family = "gaussian", 
                            # mcmc_list = mcmc_list,
                            marss = marss_list)

mcmc_list  = list(n_mcmc = 3000, n_burn = 500, n_chain = 3, n_thin = 1,step_size=0.4,adapt_delta=0.9)

#=== Run model
file <- file.path(cmdstan_path(), "pinnipeds", "marss_MA1_cmd.stan")
mod <- cmdstan_model(file)
MA1_fit <- mod$sample(
  data = data_list_tmp$data,
  seed = 124,
  chains = mcmc_list$n_chain,
  parallel_chains = mcmc_list$n_chain,
  iter_warmup = mcmc_list$n_burn,
  iter_sampling = mcmc_list$n_mcmc,
  adapt_delta=0.97,
  step_size=0.05,
  refresh = 500 # print update every 500 iters
)

preds = MA1_fit$summary(variables = "pred")
