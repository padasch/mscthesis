# START OF SCRIPT ----

# SOURCE ALL ----
suppressMessages(source("~/projects/mscthesis/scripts/final/source_fct_pkg.R"))

# ADD SCRIPT ----

## LOAD NEEDED DATA ----
## Driver and Data
df_drivers_p21    <- readRDS("~/data/mscthesis/final/df_drivers_p21.rds")
df_evaluation_p21 <- readRDS("~/data/mscthesis/raw/leaf_traits_peng2021/df_P21_clean_with_cl.rds") %>%
    dplyr::select(-author) %>%
    rename_with(tolower) %>%
    dplyr::filter(!(is.na(vcmax) & is.na(jmax)),
                  !is.na(date),
                  !is.na(lat),
                  !is.na(lon)) %>%
    group_by(site, date) %>%
    nest() %>%
    ungroup() %>%
    rename(sitename = site,
           data_raw = data)

# ......................................................................... ####
settings <- get_settings()
settings$rpmodel_accl$method_jmaxlim <- "farquhar89"
t_start <- Sys.time()
df_dummy <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
t_ana <- difftime(Sys.time(), t_start, units = "min") %>% round(2)


## Numerical, no EB
settings$rpmodel_accl$method_optim <- "numerical"
t_start <- Sys.time()
df_noeb <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
t_noeb <- difftime(Sys.time(), t_start, units = "min") %>% round(2)
p_noeb  <- plot_two_long_df(
    df_x = make_long_df(df_in = df_noeb, dataset = "no_energy_balance"),
    df_x_dataset = "no_energy_balance",
    df_y = make_df_evaluation_p21_long(df_in = df_noeb),
    df_y_dataset = "peng21",
    kphio_vcmax_corr = F
)



## EB via plantecophys
settings$rpmodel_accl$method_eb <- "plantecophys"
t_start <- Sys.time()
df_eb_pl <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
t_eb_pl <- difftime(Sys.time(), t_start, units = "min") %>% round(2)
p_eb_pl  <- plot_two_long_df(
    df_x = make_long_df(df_in = df_eb_pl, dataset = "eb_plantecophys"),
    df_x_dataset = "eb_plantecophys",
    df_y = make_df_evaluation_p21_long(df_in = df_eb_pl),
    df_y_dataset = "peng21",
    kphio_vcmax_corr = F
)


## EB via tealeaves
t_start <- Sys.time()
settings$rpmodel_accl$method_eb <- "tealeaves"
t_start <- Sys.time()
df_eb_te <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
t_eb_te <- difftime(Sys.time(), t_start, units = "min") %>% round(2)
p_eb_te  <- plot_two_long_df(
    df_x = make_long_df(df_in = df_eb_te, dataset = "eb_tealeaves"),
    df_x_dataset = "eb_tealeaves",
    df_y = make_df_evaluation_p21_long(df_in = df_eb_te),
    df_y_dataset = "peng21",
    kphio_vcmax_corr = F
)


## Computing time (run in console directly is quicker than in chunk!)
cat(" Computation time: \n",
    "_______________________________ \n",
    "Analytical Model:      ", paste0(t_ana), "min. \n", #ifelse(t_ana > 1, paste0(t_ana, " min"), paste0(t_ana*60, " sec")), "\n",
    "Numerical Model:       ", paste0(t_noeb), "min. \n", #ifelse(t_noeb > 1, paste0(t_noeb, " min"), paste0(t_noeb*60, " sec")), "\n",
    "Num. + Planteco Model: ", paste0(t_eb_pl), "min. \n", #ifelse(t_eb_pl > 1, paste0(t_eb_pl, " min"), paste0(t_eb_pl*60, " sec")), "\n",
    "Num. + tealeaves Model:", paste0(t_eb_te), "min. \n") #ifelse(t_eb_te > 1, paste0(t_eb_te, " min"), paste0(t_eb_te*60, " sec")), "\n")

# END OF SCRIPT ----
beep(3)
