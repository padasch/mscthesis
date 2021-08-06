get_settings <- function(){
    settings <- list(
        
        ## Global options
        verbose = T,
        
        ## Acclimated P-Model setup
        
        rpmodel_accl = list(method_optim      = "analytical",     # Options: analytical or numerical
                            method_jmaxlim    = "smith37",        # Options: smith37 (i.e. wang17) or farquhar89 (i.e. smith19)
                            method_ftemp      = "kumarathunge19", # Options: kattge07 or kumarathunge2019
                            energy_balance    = "off",             # Options: on or off, only applies if optimal_method = "numeric"
                            method_eb         = "plantecophys",   # Options: plantecophys or tealeaves package (latter takes significantly longer!)
                            kphio_calib       = 0.09423773,       # Options: numeric [0, 1], calibrated via rsofun v3.3
                            apar_soilm_calib  = 0.33349283,       # Options: numeric, calibrated via rsofun v3.3
                            bpar_soilm_calib  = 1.45602286,       # Options: numeric, calibrated via rsofun v3.3
                            tau = list(tc_air = 30,               # Options: numeric in days
                                       ppfd   = 30,               # Options: numeric in days
                                       vpd    = 30,               # Options: numeric in days
                                       patm   = 30,               # Options: numeric in days
                                       co2    = 30)               # Options: numeric in days
                            ),
        
        ## Instantaneous P-Model setup
        rpmodel_inst = list(method_ftemp      = "kumarathunge19", # Options: kattge07 or kumarathunge19
                            method_rd25       = "atkin15",        # Options: atkin15 or kumarathunge19
                            method_rd         = "heskel16",       # Options: q10, arrhenius or heskel2016
                            method_jmaxlim    = "smith37",        # Options: smith37 or farquhar89
                            method_vcmax25    = "rpmodel_accl",   # Options: prescribed or rpmodel_accl
                            method_ci         = "prentice14",     # Options: prescribed numeric or prentice14
                            q10               = 2),               # Options: numeric
        
        ## Simulation P-Model setup
        rpmodel_exp = list(ppfd               = "metainfo")       # Options: "metainfo" or numeric in mol/m2/s
        )
    
    return(settings)
}

# .............................................................................


run_rpmodel_accl <- function(settings = NA,     # Options: Setting to NA takes default settings from get_settings()
                             df_drivers,        # Options: Input of drivers for sites, needs: sitename (chr) and forcing (nested tibble, see rsofun v3.3 setup)
                             df_evaluation,     # Options: Input of evaluation sitename, date and target value: sitename | date | target_var
                             target_var){       # Options: Character specifying target_var

    
    ## Getting default settings if none are given:
    if (length(settings) == 1) {
        settings <- get_settings()
    }
    
    ## Dampen drivers to get acclimated values:
    df_drivers <- df_drivers %>% 
        mutate(forcing = purrr::map(forcing, ~mutate(., tc_growth_air = dampen_vec(temp, settings$rpmodel_accl$tau$tc_air))),
               forcing = purrr::map(forcing, ~mutate(., ppfd_growth   = dampen_vec(ppfd, settings$rpmodel_accl$tau$ppfd))),
               forcing = purrr::map(forcing, ~mutate(., patm_growth   = dampen_vec(patm, settings$rpmodel_accl$tau$patm))),
               forcing = purrr::map(forcing, ~mutate(., co2_growth    = dampen_vec(co2, settings$rpmodel_accl$tau$co2))),
               forcing = purrr::map(forcing, ~mutate(., vpd_growth    = dampen_vec(vpd, settings$rpmodel_accl$tau$vpd))))
    
    ## Reducing dataframe to sites where acclimated variables should be calculated:
    df_rpmodel_accl <- df_drivers %>%
        unnest(c("forcing", "siteinfo")) %>%
        mutate(
            chi = NA,
            ci = NA,
            xi = NA,
            gs = NA,
            vcmax = NA,
            vcmax25 = NA,
            jmax = NA,
            jmax25 = NA,
            kphio = NA) %>% 
        right_join(df_evaluation) %>% 
        drop_na(ppfd) # To make sure, no rows with NA's enter loop below
    
    ## Calculating acclimated variables:
    for (row in 1:nrow(df_rpmodel_accl)) {
        message('\014')
        message("rpmodel_accl: ", round(row / nrow(df_rpmodel_accl) * 100), " %")
        
        df_loop <- rpmodel(## Acclimated inputs
                           tc_growth_air = df_rpmodel_accl$tc_growth_air[row],
                           tc_home       = df_rpmodel_accl$tc_home[row],
                           vpd           = df_rpmodel_accl$vpd_growth[row],
                           co2           = df_rpmodel_accl$co2_growth[row],
                           ppfd          = df_rpmodel_accl$ppfd_growth[row],
                           patm          = df_rpmodel_accl$patm_growth[row],
                           
                           ## Local parameters
                           kphio         = settings$rpmodel_accl$kphio_calib,
                           apar_soilm    = settings$rpmodel_accl$apar_soilm_calib,
                           bpar_soilm    = settings$rpmodel_accl$bpar_soilm_calib,
                           
                           ## Calculation methods
                           method_optim   = settings$rpmodel_accl$method_optim, 
                           method_jmaxlim = settings$rpmodel_accl$method_jmaxlim,
                           method_ftemp   = settings$rpmodel_accl$method_ftemp,
                           energy_balance = settings$rpmodel_accl$energy_balance,
                           method_eb      = settings$rpmodel_accl$method_eb)
        
        df_rpmodel_accl$chi[row]     <- df_loop$chi
        df_rpmodel_accl$ci[row]      <- df_loop$ci
        df_rpmodel_accl$xi[row]      <- df_loop$xi
        df_rpmodel_accl$gs[row]      <- df_loop$gs
        df_rpmodel_accl$vcmax[row]   <- df_loop$vcmax
        df_rpmodel_accl$vcmax25[row] <- df_loop$vcmax25
        df_rpmodel_accl$jmax[row]    <- df_loop$jmax
        df_rpmodel_accl$jmax25[row]  <- df_loop$jmax25
        df_rpmodel_accl$kphio[row]   <- df_loop$kphio
    }
    
    ## Return final dataframe, ready for instantaneous P-Model
    out <- df_rpmodel_accl %>%
        ## Get nested acclimated data
        dplyr::select(sitename, date, chi, ci, xi, gs, vcmax, vcmax25, jmax, jmax25, kphio) %>% 
        group_by(sitename, date) %>% 
        nest() %>% 
        ungroup() %>% 
        rename(rpm_accl = data) %>% 
        
        ## Get nested forcing data
        left_join(df_drivers %>%
                      dplyr::select(sitename, forcing, siteinfo) %>% 
                      unnest(c(forcing, siteinfo)) %>% 
                      dplyr::select(c(df_drivers$forcing[[1]] %>% names()), "tc_home", "sitename")) %>%
        
        group_by(sitename, date, rpm_accl) %>% 
        nest() %>% 
        ungroup() %>% 
        rename(forcing = data) %>% 
        
        ## Get nested evaluation data
        left_join(df_evaluation)
    
    return(out)
}

# .............................................................................
#---------------------------------#
# Function to get anet ~ tc_leaf  #
#---------------------------------#
sim_rpmodel <- function(df_in, settings){
    
    ## Define empty arrays
    tc_leaf_array    <- seq(0, 40, 0.1 )
    agross_array    <- rep(NA, length(tc_leaf_array))
    anet_array	    <- rep(NA, length(tc_leaf_array))
    vcmax_array     <- rep(NA, length(tc_leaf_array))
    jmax_array      <- rep(NA, length(tc_leaf_array))
    rd_array        <- rep(NA, length(tc_leaf_array))
    ci_array        <- rep(NA, length(tc_leaf_array))
    gs_array        <- rep(NA, length(tc_leaf_array))
    ac_array        <- rep(NA, length(tc_leaf_array))
    aj_array        <- rep(NA, length(tc_leaf_array))
    min_a_array     <- rep(NA, length(tc_leaf_array))
    dummy_array     <- rep(NA, length(tc_leaf_array))
    
    ## Define input for instantaneous P-Model
    df_out <- df_in %>% 
        mutate(rpm_sim = list(tibble()))
    
    df_sim <- df_out %>%
        unnest(c(rpm_accl, forcing, metainfo)) 
    
    if (settings$rpmodel_exp$ppfd == "metainfo") {
        df_sim$ppfd <- df_sim$ppfd_measurement
    }
    
    if (settings$rpmodel_inst$method_vcmax25 == "prescribed") {
        df_sim$vcmax25 <- df_sim$vcmax25_pft
        
        # Acclimated jmax_to_vcmax taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment:
        df_sim$jmax25  <- df_sim$vcmax25 * (2.56 - (0.0375 * df_sim$tc_home) + (-0.0202 * (df_sim$tc_growth_air - df_sim$tc_home)))
    }
    
    ## Get starting time of computation
    start.time <- Sys.time()
    
    ## Start loop over sites
    for (site in 1:nrow(df_sim)) {
        
        if (settings$verbose){
            message('\014')
            message("rpmodel_sim: ", round(site / nrow(df_sim) * 100), "%")
        }
        
        sitename_temp <- df_sim$sitename[site]
        date_temp     <- df_sim$date[site]
        
        
        for(i in 1:length(tc_leaf_array)){
            
            # if (settings$verbose){
            #     message('\014')
            #     message("rpmodel_sim: Slice: ", round(site / nrow(df_sim) * 100), "% | Run: ", round(i / length(tc_leaf_array) * 100), "% | \n")
            # }
            
            out_pmodel_inst_opt = rpmodel_inst(vcmax25    = df_sim$vcmax25[site],
                                               jmax25     = df_sim$jmax25[site],
                                               xi         = df_sim$xi[site],
                                               gs         = df_sim$gs_accl[site],
                                               tc_leaf    = tc_leaf_array[i],
                                               vpd        = df_sim$vpd[site],
                                               co2        = df_sim$co2[site],
                                               fapar      = df_sim$fapar[site],
                                               ppfd       = df_sim$ppfd[site],
                                               patm       = df_sim$patm[site],
                                               kphio      = df_sim$kphio[site],
                                               tc_growth  = df_sim$tc_growth_air[site],
                                               tc_home    = df_sim$tc_home[site],
                                               
                                               # settings
                                               settings    = settings
            )
            
            anet_array[i]      = out_pmodel_inst_opt$anet
            agross_array[i]    = out_pmodel_inst_opt$agross
            vcmax_array[i]     = out_pmodel_inst_opt$vcmax
            jmax_array[i]      = out_pmodel_inst_opt$jmax
            rd_array[i]        = out_pmodel_inst_opt$rd
            ci_array[i]        = out_pmodel_inst_opt$ci
            gs_array[i]        = out_pmodel_inst_opt$gs
            ac_array[i]        = out_pmodel_inst_opt$ac
            aj_array[i]        = out_pmodel_inst_opt$aj
            min_a_array[i]     = out_pmodel_inst_opt$min_a
            dummy_array[i]     = out_pmodel_inst_opt$dummy
        }
        
        # Return final tibble with simulated values
        df_out$rpm_sim[site] <- list(tibble(tc_leaf_sim = tc_leaf_array,
                                            anet_sim    = anet_array,
                                            agross_sim  = agross_array,
                                            vcmax_sim   = vcmax_array,
                                            jmax_sim    = jmax_array,
                                            rd_sim      = rd_array,
                                            ci_sim      = ci_array,
                                            gs_sim      = gs_array,
                                            ac_sim      = ac_array,
                                            aj_sim      = aj_array,
                                            min_a_sim   = min_a_array,
                                            dummy_sim   = dummy_array))
    }
    
    return(df_out)
}

extract_tc_opt_sim <- function(df_in){
    df_temp <- df_in %>%
        mutate(tc_opt_sim    = purrr::map_dbl(rpm_sim, ~.$tc_leaf_sim[which.max(.$anet_sim)]),
               tc_opt_ac_sim = purrr::map_dbl(rpm_sim, ~.$tc_leaf_sim[which.max(.$ac_sim)]),
               tc_opt_aj_sim = purrr::map_dbl(rpm_sim, ~.$tc_leaf_sim[which.max(.$aj_sim)]),
               min_a_sim     = purrr::map_chr(rpm_sim, ~.$min_a_sim[which.max(.$agross_sim)]),
               min_a_sim     = as.factor(min_a_sim),
               tc_growth_air = purrr::map_dbl(forcing, ~dplyr::select(., tc_growth_air) %>% pull)) %>% 
        
        # Get a_growth from predictions
        mutate(agrowth_sim  = purrr::map_dbl(rpm_sim, ~.$anet_sim[ which(round(.$tc_leaf_sim, 1) == round(tc_growth_air, 1)) ] )) %>%
        
        # Get a_growth from data
        unnest(fit_opt) %>% 
        mutate(agrowth_dat  = aopt - ( b * ( tc_growth_air - tc_opt_obs) ^ 2)) %>%
        
        # Clean up
        dplyr::select(sitename, date, tc_opt_sim, tc_opt_ac_sim, tc_opt_aj_sim, min_a_sim, agrowth_sim, agrowth_dat) %>% 
        group_by(sitename, date) %>% 
        nest() %>% 
        ungroup() %>% 
        rename(fit_sim = data)
    
    df_out <- left_join(df_in, df_temp)
}


#---------------------------------#
# Function for instant rpmodel    #
#---------------------------------#
rpmodel_inst <- function(vcmax25, jmax25, xi, gs, tc_leaf, vpd, co2, fapar, ppfd, patm, kphio, tc_growth, tc_home, settings){
    
    ## Prepare output:
    out <- list(
        agross = NA,
        anet   = NA,
        vcmax  = NA,
        jmax   = NA,
        rd     = NA,
        ci     = NA,
        gs     = NA,
        aj     = NA,
        ac     = NA,
        min_a  = NA,
        dummy  = NA)
    
    dummy_var <- NA
    
    ## Get instantaneous conditions:
    # Instantaneous partial pressure of CO2
    ca = co2_to_ca(co2, patm)

    # Instantaneous gammastar:
    gammastar = calc_gammastar(tc_leaf, patm)
    
    # Instantaneous K:
    kmm = calc_kmm(tc_leaf, patm)
    
    # Instantaneous viscosity:
    ns_star <- calc_viscosity_h2o( tc_leaf, patm ) / calc_viscosity_h2o( tc = 25.0, p = 101325.0 )
    
    # Instantaneous vcmax and jmax:
    vcmax  <- vcmax25 * calc_ftemp_inst_vcmax(tc_leaf, tc_growth,          tcref = 25.0, method_ftemp = settings$rpmodel_inst$method_ftemp)
    jmax   <- jmax25  * calc_ftemp_inst_jmax( tc_leaf, tc_growth, tc_home, tcref = 25.0, method_ftemp = settings$rpmodel_inst$method_ftemp)
    
    
    ## Call analytical or numerical model
    
    if (T) { #(settings$rpmodel_accl$method_optim == "analytical") { # TODO: ALWAYS USING ANALYTICAL VERSION BECAUSE UNSURE WHAT GS-PREDICTION I SHOULD TAKE
        
        # Instantaneous leaf-internal CO2: ci
        ci <- calc_ci(ca, gammastar, xi, vpd, patm, settings$rpmodel_inst$method_ci)
        
        # Electron Transport Rate: Aj
        aj <- calc_aj(kphio, ppfd, jmax, gammastar, ci, ca, fapar, j_method = settings$rpmodel_inst$method_jmaxlim, model = "analytical")$aj
        
        # Carboxylation Rate: Ac
        ac <- calc_ac(ci, gammastar, kmm, vcmax, model = "analytical")$ac
        
    } else if (settings$rpmodel_accl$method_optim == "numerical") {
        
        # Electron Transport Rate: Aj
        aj_out <- calc_aj(kphio, ppfd, jmax, gammastar, ci, ca, fapar, j_method = settings$rpmodel_inst$method_jmaxlim, model = settings$rpmodel_accl$method_optim)
        aj     <- aj_out$aj
        ci_j   <- aj_out$ci
        
        # Carboxylation Rate: Ac
        ac_out <- calc_ac(ci, gammastar, kmm, model = settings$rpmodel_accl$method_optim)
        ac     <- ac_out$ac
        ci_c   <- ac_out$ci
        
        # Instantaneous leaf-internal CO2: ci
        ci <- max(ci_c, ci_j)
    }
    
    # Dark Respiration: Rd
    rd <- calc_rd(tc_leaf, vcmax25, tc_growth, q10 = settings$rpmodel_inst$q10, method_rd25 = settings$rpmodel_inst$method_rd25, method_rd_scale = settings$rpmodel_inst$method_rd)
    
    # Gross assimilation rate: A_gross
    agross <- min(aj, ac)
    
    if (aj == ac) {
        min_a <- "colimit"
    } else if (agross == ac) {
        min_a <- "ac"
    } else if (agross == aj) {
        min_a <- "aj"
    }
    
    # Net assimilation rate: A_net
    anet <- agross - rd
    
    # Stomatal conductance: gs
    gs <- aj/(ca - ci)
    
    # Definition of output:
    out <- list(
        agross = agross,
        anet = anet,
        vcmax = vcmax,
        jmax = jmax,
        rd = rd,
        ci = ci,
        gs = gs,
        aj = aj,
        ac = ac,
        min_a = min_a,
        dummy = dummy_var)
        
    return(out)
}

# .............................................................................
run_accl_to_sim <- function(settings,
                            df_drivers,
                            df_evaluation,
                            target_var) {
    
    df_out <- run_rpmodel_accl(settings = settings,
                     df_drivers = df_drivers,
                     df_evaluation = df_evaluation,
                     target_var = target_var) %>% 
        
        sim_rpmodel(df_in = .,
                    settings = settings) %>% 
        
        extract_tc_opt_sim(df_in = .)
    
    return(df_out)
}

# .............................................................................

plot_obs_vs_pred_mina <- function(df_in) {
    
    # df_temp <- df_in %>% unnest(fit_sim)
    
    fit     <- lm(tc_opt_obs ~ tc_opt_sim, data = df_temp)
    sry     <- summary(fit)
    r2      <- sry$adj.r.squared %>% round(2)
    rmse    <- sqrt( ( c(crossprod(sry$residuals)) / length(sry$residuals) ) ) %>% round(2)
    slope   <- sry$coefficients[2, 1] %>% round(2)
    b0_sign <- sry$coefficients[1, 4] <= 0.05 # T = Intercept is signficantly different from 0
    b1_sign <- !(between(1, confint(fit, 'tc_opt_sim', level=0.95)[1] %>% round(2), confint(fit, 'tc_opt_sim', level=0.95)[2]) %>% round(2)) # T = Slope is significantly different from 1
    bias    <- mean(df_temp$tc_opt_sim - df_temp$tc_opt_obs) %>% round(2)
    
    n_ac <- df_temp %>% dplyr::select(sitename, date, tc_opt_sim, min_a_sim) %>% distinct() %>% dplyr::filter(min_a_sim == "ac") %>% nrow()
    n_aj <- df_temp %>% dplyr::select(sitename, date, tc_opt_sim, min_a_sim) %>% distinct() %>% dplyr::filter(min_a_sim == "aj") %>% nrow()
    
    df_metrics <- tibble(r2 = r2,
                         rmse = rmse,
                         slope = slope,
                         b0_sign_not_0   = b0_sign,
                         b1_sign_not_one = b1_sign,
                         bias = bias,
                         n_minac = n_ac,
                         n_minaj = n_aj)
    
    ## Check
    
    plot <- df_temp %>%
        # TODO Add non-sense samples to have at least one sample of each limitation in data frame to make plotting work. Does not affect plot nor fit.
        add_row(min_a_sim = "aj", tc_opt_obs = -1, tc_opt_sim = -1, .before = 1) %>%
        add_row(min_a_sim = "ac", tc_opt_obs = -1, tc_opt_sim = -1, .before = 1) %>%
        ggplot() +
        geom_abline(linetype = "dotted") +
        # geom_point(aes(y = tc_opt_obs, x = tc_opt_sim, color = dataset_yr, shape = min_a_sim)) +
        # geom_point(aes(y = tc_opt_obs, x = tc_opt_sim, shape = min_a_sim, color = min_a_sim), alpha = 0.5) +
        
        # Separating for Ac and Aj differentiation
        geom_point(data = . %>% dplyr::filter(min_a_sim == "ac"), aes(y = tc_opt_obs, x = tc_opt_sim, color = "Ac"))+
        geom_point(data = . %>% dplyr::filter(min_a_sim == "aj"), aes(y = tc_opt_obs, x = tc_opt_sim, color = "Aj")) +
        
        geom_smooth(aes(y = tc_opt_obs, x = tc_opt_sim), method = "lm", size = 0.5, fullrange = T, color = "black") +
        xlim(0, 40) +
        ylim(0, 40) +
        xlab("Predicted T_opt [°C]") +
        ylab("Observed T_opt [°C]") +
        scale_shape_manual(name = "shape", values = c(4, 16)) +
        scale_color_manual(name = "min(Ac, Aj)", values = c(Ac = "grey50", Aj = "tomato")) +
        # labs(color = "Dataset",
        labs(color = "Limit.",
             shape = "Limit.",
             title = "to be added",
             subtitle = bquote(R^2  ~  " = "  ~ .(r2)  ~  " | bias = "  ~ .(bias)  ~  " | slope = "  ~ .(slope)),
             caption = paste0("Signficant different b0 ~ 0: ", ifelse(b0_sign, "N", "Y"), " | b1 ~ 1: ", ifelse(b1_sign, "N", "Y"),
                              " | Lim: Ac = ", n_ac, " Aj = ", n_aj,
                              " | RMSE = ", rmse)
        )
    
    out <- list(plot = plot,
                df_metrics   = df_metrics)
    
    return(out)
}
# .............................................................................

plot_obs_vs_pred <- function(df_in, x, y) {
    
    df_temp <- df_in
    
    df_temp$y <- df_temp[[y]]
    df_temp$x <- df_temp[[x]]
    
    fit     <- lm(x ~ y, data = df_temp)
    sry     <- summary(fit)
    r2      <- sry$adj.r.squared %>% round(2)
    rmse    <- sqrt( ( c(crossprod(sry$residuals)) / length(sry$residuals) ) ) %>% round(2)
    slope   <- sry$coefficients[2, 1] %>% round(2)
    b0_sign <- sry$coefficients[1, 4] <= 0.05 # T = Intercept is signficantly different from 0
    b1_sign <- !(between(1, confint(fit, x, level=0.95)[1] %>% round(2), confint(fit, x, level=0.95)[2]) %>% round(2)) # T = Slope is significantly different from 1
    bias    <- mean(df_temp$y - df_temp$x) %>% round(2)
    
    df_metrics <- tibble(r2 = r2,
                         rmse = rmse,
                         slope = slope,
                         b0_sign_not_0   = b0_sign,
                         b1_sign_not_one = b1_sign,
                         bias = bias)
    
    plot <- df_temp %>%
        ggplot() +
        aes(y = y, x = x) +
        geom_abline(linetype = "dotted") +
        geom_smooth(method = "lm", size = 0.5, fullrange = T, color = "black") +
        geom_point() +
        ylab(paste(y)) +
        xlab(paste(x)) +
        labs(subtitle = bquote(R^2  ~  " = "  ~ .(r2)  ~  " | bias = "  ~ .(bias)  ~  " | slope = "  ~ .(slope) ~  " | RMSE = "  ~ .(rmse)),
             caption = paste0("Signficant different b0 ~ 0: ", ifelse(b0_sign, "N", "Y"), " | b1 ~ 1: ", ifelse(b1_sign, "N", "Y"))
        ) +
        theme_bw()
    
    out <- list(plot = plot,
                df_metrics   = df_metrics)
    
    return(out)
}


# TEST ZONE ####

if (F) {
    
    load("~/data/rsofun_benchmarking/df_drivers_topt_v3.4.Rdata")
    df_evaluation <- readRDS("~/data/mscthesis/final/df_obs_opt_postQC.rds")
    df_drivers     <- df_drivers_topt
    target_var     <- "tc_opt_obs"
   
    ## Individual variables
    tc_growth_air = 20
    tc_home = 20
    vpd = 2000
    co2 = 400
    fapar = 1
    ppfd = 0.0015
    patm = 10^5 
    elv = 10
    kphio = 0.1
    beta = 146.0 
    meanalpha = 1.0
    apar_soilm = 0.0
    bpar_soilm = 0.73300
    c4 = FALSE
    ## Calculation methods
    method_optim   = "analytical"
    method_optci   = "prentice14"
    method_jmaxlim = "smith37"
    method_ftemp   = "kumarathunge19"
    method_eb      = "plantecophys"
    energy_balance = "on"
    ## Other settings
    do_ftemp_kphio = TRUE
    do_soilmstress = FALSE
    returnvar = NULL
    verbose = FALSE
    
     
    rpmodel(tc_growth_air = 20,
    tc_home = 20,
    vpd = 2000,
    co2 = 400,
    fapar = 1,
    ppfd = 500,
    patm = 10^5 ,
    elv = 10,
    kphio = 0.1,
    beta = 146.0 ,
    meanalpha = 1.0,
    apar_soilm = 0.0,
    bpar_soilm = 0.73300,
    c4 = FALSE,
    ## Calculation methods,
    method_optim   = "analytical",
    method_optci   = "prentice14",
    method_jmaxlim = "smith37",
    method_ftemp   = "kumarathunge19",
    ## Other settings,
    do_ftemp_kphio = TRUE,
    do_soilmstress = FALSE,
    returnvar = NULL,
    verbose = FALSE)
    
    
    
    test <- tibble(site = c("a", "b"),
                   vars = list(tibble(tc_growth_air = 20,
                         tc_home = 20,
                         vpd = 2000,
                         co2 = 400,
                         fapar = 1,
                         ppfd = 500,
                         patm = 10^5 ,
                         elv = 10,
                         kphio = 0.1,
                         beta = 146.0 ,
                         meanalpha = 1.0,
                         apar_soilm = 0.0,
                         bpar_soilm = 0.73300,
                         c4 = FALSE,
                         ## Calculation methods,
                         method_optim   = "analytical",
                         method_optci   = "prentice14",
                         method_jmaxlim = "smith37",
                         method_ftemp   = "kumarathunge19",
                         ## Other settings,
                         do_ftemp_kphio = TRUE,
                         do_soilmstress = FALSE,
                         returnvar = NULL,
                         verbose = FALSE),
                         
                    tibble(tc_growth_air = 20,
                         tc_home = 20,
                         vpd = 2000,
                         co2 = 400,
                         fapar = 1,
                         ppfd = 500,
                         patm = 10^5 ,
                         elv = 10,
                         kphio = 0.1,
                         beta = 146.0 ,
                         meanalpha = 1.0,
                         apar_soilm = 0.0,
                         bpar_soilm = 0.73300,
                         c4 = FALSE,
                         ## Calculation methods,
                         method_optim   = "analytical",
                         method_optci   = "prentice14",
                         method_jmaxlim = "smith37",
                         method_ftemp   = "kumarathunge19",
                         ## Other settings,
                         do_ftemp_kphio = TRUE,
                         do_soilmstress = FALSE,
                         returnvar = NULL,
                         verbose = FALSE)))
    
    
    test2 <- test %>% mutate(out = purrr::map(vars, ~rpmodel(tc_growth_air = tc_growth_air,
                                                   tc_home = tc_home,
                                                   vpd = vpd,
                                                   co2 = co2,
                                                   fapar = fapar,
                                                   ppfd = ppfd,
                                                   patm = 10^5 ,
                                                   elv = elv,
                                                   kphio = settings$rpmodel_accl$kphio_calib,
                                                   beta = beta ,
                                                   meanalpha = meanalpha,
                                                   apar_soilm = apar_soilm,
                                                   bpar_soilm = bpar_soilm,
                                                   c4 = c4,
                                                   ## Calculation##methods,
                                                   method_optim   = method_optim,
                                                   method_optci   = method_optci,
                                                   method_jmaxlim = method_jmaxlim,
                                                   method_ftemp   = method_ftemp,
                                                   ## Other##settings,
                                                   do_ftemp_kphio = do_ftemp_kphio,
                                                   do_soilmstress = do_soilmstress,
                                                   returnvar = returnvar,
                                                   verbose = verbose)))
    
}





if (F) {
    
    # predicted vs. observed
    p1 <- ggplot(df_exp) +
        geom_abline(linetype = "dotted") +
        geom_point(aes(y = tc_opt_obs, x = tc_opt_rpm, color = sitename)) +
        geom_smooth(aes(y = tc_opt_obs, x = tc_opt_rpm), method = "lm", size = 0.5, fullrange = T, color = "black") +
        xlim(0, 40) +
        ylim(0, 40)
    
    # predicted vs. tc_growth
    p2 <- ggplot(df_exp) +
        geom_abline(linetype = "dotted") +
        geom_point(aes(y = tc_opt_rpm, x = tc_growth, color = sitename)) +
        geom_smooth(aes(y = tc_opt_rpm, x = tc_growth), method = "lm", size = 0.5, fullrange = T, color = "black") +
        xlim(0, 40) +
        ylim(0, 40)
    
    # observed vs. tc_growth
    p3 <- ggplot(df_exp) +
        geom_abline(linetype = "dotted") +
        geom_point(aes(y = tc_opt_obs, x = tc_growth, color = sitename)) +
        geom_smooth(aes(y = tc_opt_obs, x = tc_growth), method = "lm", size = 0.5, fullrange = T, color = "black") +
        xlim(0, 40) +
        ylim(0, 40)
    
    # patchwork
    p_multi <- (p1 / p2 / p3 ) +
        plot_layout(guides='collect') &
        theme(legend.position = "bottom",
              legend.key.size = unit(0.25, "cm")) &
        plot_annotation(tag_levels = 'A',
                        title = paste0(settings$NAME),
                        subtitle = paste0("Pred-Obs: adj. R^2 = ", round(fit_exp$adj.r.squared, 3)))
}
