# read-in data -----------------------------------------------------------------

raw_data = read.csv(here::here("Data/PELGAS/raw_data.csv"),sep = ";") 

# extract years ----------------------------------------------------------------

raw_data$year = as.numeric(substr(raw_data$cruise,7,10))

# keep species with 10 + data points only ---------------------------------------

raw_data = raw_data %>% 
  group_by(spname) %>% 
  filter(n() >= 10)

# export data ------------------------------------------------------------------
write.csv(raw_data, here::here("Data/PELGAS/pelgas_clean.csv"), row.names = FALSE)