suppressMessages(source("~/projects/mscthesis/scripts/final/source_fct_pkg.R"))

settings <- get_settings()
settings$rpmodel_accl$method_optim <- "numerical"
settings$rpmodel_accl$energy_balance <- "on"
settings$rpmodel_accl$method_eb      <- "tealeaves"
settings$rpmodel_accl$method_jmaxlim <- "farquhar89"
settings$rpmodel_accl$kphio_calib <- settings$rpmodel_accl$kphio_calib*4

load("~/data/rsofun_benchmarking/df_drivers_topt_v3.4.Rdata")
df_drivers_k19     <- df_drivers_topt
df_evaluation_k19  <- readRDS("~/data/mscthesis/final/df_obs_opt_postQC.rds")

df_background <- run_accl_to_sim(settings, df_drivers_k19, df_evaluation_k19)
