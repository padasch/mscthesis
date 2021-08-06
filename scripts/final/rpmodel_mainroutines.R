# MAIN RPMODEL ####
rpmodel <- function(
    ## Forcings
    tc_growth_air,
    tc_home,
    vpd,
    co2,
    fapar = 1,
    ppfd,
    patm = NA,
    elv = NA,
    
    ## Local parameters
    # kphio = ifelse(do_ftemp_kphio, ifelse(do_soilmstress, 0.087182, 0.081785), 0.049977),
    kphio,
    beta = 146.0,
    soilm = stopifnot(!do_soilmstress),
    meanalpha = 1.0,
    apar_soilm = 0.0,
    bpar_soilm = 0.73300,
    c4 = FALSE,
    soilmstress = 1,
    
    ## Calculation methods
    method_optim   = "analytical",
    method_optci   = "prentice14",
    method_jmaxlim = "wang17",
    method_ftemp   = "kumarathunge19",
    energy_balance = "off",
    method_eb      = "plantecophys",
    
    ## Other settings
    do_ftemp_kphio = TRUE,
    do_soilmstress = FALSE,
    returnvar = NULL,
    verbose = FALSE 
){
  
  ## DEBUG ZONE ####
  ## Overwriting kphio input to check for its effect:
  if (F) {
    if (method_jmaxlim == "wang17" | method_jmaxlim == "smith19") {kphio <- 0.085}
    if (method_jmaxlim == "smith19" | method_jmaxlim == "farquhar89") {kphio <- 0.275}
  }
  
  
  ## Check arguments ####
  if (identical(NA, elv) && identical(NA, patm)){
      stop("Aborted. Provide either elevation (arugment elv) or atmospheric pressure (argument patm).")
  } else if (!identical(NA, elv) && identical(NA, patm)){
      if (verbose) warning("Atmospheric pressure (patm) not provided. Calculating it as a function of elevation (elv), assuming standard atmosphere (101325 Pa at sea level).")
      patm <- patm(elv)
  }

  
  ## Get constants ####
  c_molmass <- 12.0107  # molecular mass of carbon (g)
  kPo <- 101325.0       # standard atmosphere, Pa (Allen, 1973)
  kTo <- 25.0           # base temperature, deg C (Prentice, unpublished)
  rd_to_vcmax <- 0.015  # Ratio of Rdark to Vcmax25, number from Atkin et al., 2015 for C3 herbaceous
  
  ## Get energy balance ####
  if (method_optim == "numerical" && energy_balance == "on") {
      
      if (method_jmaxlim == "wang17") {method_jmaxlim <- "smith37"}
      if (method_jmaxlim == "smith19") {method_jmaxlim <- "farquhar89"}
      
      tc_growth_leaf <- calc_tc_leaf_final(tc_air    = tc_growth_air,
                                           tc_home   = tc_home,
                                           ppfd      = ppfd * 3600 * 24, # Input has to be in mol/m2/d
                                           fapar     = fapar,
                                           co2       = co2,
                                           patm      = patm,
                                           vpd       = vpd,
                                           kphio     = kphio,
                                           method_jmaxlim_inst = method_jmaxlim,
                                           method_eb = method_eb)
    } else {
        tc_growth_leaf <- tc_growth_air
    }
    
    ## Get model-independent variables ####
    if (do_ftemp_kphio){ ftemp_kphio <- calc_ftemp_kphio( tc_growth_leaf, c4 ) }
    if (do_soilmstress) { soilmstress <- soilmstress( soilm, meanalpha, apar_soilm, bpar_soilm ) }
    ca         <- co2_to_ca( co2, patm )
    gammastar  <- calc_gammastar( tc_growth_leaf, patm )
    kmm        <- calc_kmm( tc_growth_leaf, patm )
    ns_star    <- calc_viscosity_h2o( tc_growth_leaf, patm ) / calc_viscosity_h2o( kTo, kPo )  # (unitless)
    xi         <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
    kphio_accl <- kphio * ftemp_kphio
    
        
    ## Analytical model ####
    if (method_optim == "analytical") {
        

        ### Get chi ####
        ## The heart of the P-model: calculate ci:ca ratio (chi) and additional terms

        if (c4){
            out_optchi <- calc_chi_c4()
            
        } else if (method_optci == "prentice14"){
            out_optchi <- calc_optimal_chi( kmm, gammastar, ns_star, ca, vpd, beta )
            
        } else {
            stop("rpmodel(): argument method_optci not idetified.")
            
        }
        
        
        ## Get vcmax and jmax ####
        if (c4){
            
            out_lue_vcmax <- calc_lue_vcmax_c4(
                kphio,
                ftemp_kphio,
                c_molmass,
                soilmstress
            )
            
        } else if (method_jmaxlim == "wang17" | method_jmaxlim == "smith37"){
            
            ## apply correction by Jmax limitation
            out_lue_vcmax <- calc_lue_vcmax_wang17(
                out_optchi,
                kphio,
                ftemp_kphio,
                c_molmass,
                soilmstress
            )
            
            
            
        } else if (method_jmaxlim == "smith19" | method_jmaxlim == "farquhar89"){
            
            out_lue_vcmax <- calc_lue_vcmax_smith19(
                out_optchi,
                kphio,
                ftemp_kphio,
                c_molmass,
                soilmstress
            )
            
        } else if (method_jmaxlim=="none"){
            
            out_lue_vcmax <- calc_lue_vcmax_none(
                out_optchi,
                kphio,
                ftemp_kphio,
                c_molmass,
                soilmstress
            )
            
        } else {
            
            stop("rpmodel(): argument method_jmaxlim not idetified.")
            
        }
        
        ## Corrolaries ####
        chi        <- out_optchi$chi
        xi         <- out_optchi$xi
        ci         <- chi * ca # leaf-internal CO2 partial pressure (Pa)
        iwue       <- ( ca - ci ) / 1.6  # intrinsic water use efficiency (in Pa)
        
        ## Scalings with Iabs:
        iabs <- fapar * ppfd
        
        # Gross Primary Productivity
        gpp <- iabs * out_lue_vcmax$lue   # in g C m-2 s-1
        
        # Carboxylation rate
        ftemp_vcmax       <- calc_ftemp_inst_vcmax(tcleaf = tc_growth_leaf, tcgrowth = tc_growth_air, method_ftemp = method_ftemp)
        vcmax             <- iabs * out_lue_vcmax$vcmax_unitiabs
        # vcmax             <- kphio_accl * ppfd * out_optchi$mjoc * calc_mprime( out_optchi$mj ) / out_optchi$mj ## Testing formulation from rsofun directly, works better than via lue()-function
        vcmax25           <- vcmax / ftemp_vcmax
        ac                <- calc_ac(ci, gammastar, kmm, vcmax, model = method_optim)$ac
        
        ## Electron transport rate
        jmax <- calc_jmax(kphio_accl, iabs, ci, gammastar, method = settings$rpmodel_accl$method_jmaxlim)
        
        ## TODO: OLD CODE COMMENTED OUT
        ## ..............................  
        # fact_jmaxlim      <- vcmax * (ci + 2.0 * gammastar) / (kphio_accl * iabs * (ci + kmm))
        # jmax              <- 4.0 * kphio_accl * iabs / sqrt( (1.0/fact_jmaxlim)^2 - 1.0 )
        ## ..............................  
        
        ftemp_jmax        <- calc_ftemp_inst_jmax(tcleaf = tc_growth_leaf, tcgrowth = tc_growth_air, tchome = tc_home, method_ftemp = method_ftemp)
        jmax25            <- jmax / ftemp_jmax
        aj                <- calc_aj(kphio_accl, ppfd, jmax, gammastar, ci, ca, fapar, j_method = method_jmaxlim, model = method_optim)$aj
        
        ## Dark respiration
        ftemp_inst_rd      <- calc_ftemp_inst_rd( tc_growth_leaf )
        rd_unitiabs        <- rd_to_vcmax * (ftemp_inst_rd / ftemp_vcmax) * out_lue_vcmax$vcmax_unitiabs
        rd                 <- iabs * rd_unitiabs
        
        ## Gross assimilation
        assim <- ifelse(aj < ac , aj, ac)
        
        ## Average stomatal conductance
        gs <- assim / (ca - ci)
        
    } else if (method_optim == "numerical") {
        
        ## Numerical Model ####
        ## Get optimal vcmax, jmax and gs
      
        if (method_jmaxlim == "wang17") {method_jmaxlim <- "smith37"}
        if (method_jmaxlim == "smith19") {method_jmaxlim <- "farquhar89"}
        
        out_num <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_growth_leaf,
                                                  patm = patm,
                                                  co2 = co2,
                                                  vpd = vpd,
                                                  ppfd = ppfd * 3600 * 24,
                                                  fapar = fapar,
                                                  kphio = kphio,
                                                  method_jmaxlim_inst = method_jmaxlim)
        
        chi     <- out_num$chi_mine
        ci      <- out_num$ci_mine
        gs      <- out_num$gs_mine    / 3600 / 24 # Numeric output is in mol/m2/d, get it into mol/m2/s
        vcmax   <- out_num$vcmax_mine / 3600 / 24 # Numeric output is in mol/m2/d, get it into mol/m2/s
        jmax    <- out_num$jmax_mine  / 3600 / 24 # Numeric output is in mol/m2/d, get it into mol/m2/s
        vcmax25 <- vcmax /  calc_ftemp_inst_vcmax(tcleaf = tc_growth_leaf, tcgrowth = tc_growth_air, method_ftemp = method_ftemp)
        jmax25  <- jmax  /  calc_ftemp_inst_jmax(tcleaf = tc_growth_leaf, tcgrowth = tc_growth_air, tchome = tc_home, method_ftemp = method_ftemp)
    }
    
    ## Output definition ####
    out <- tibble(
        # gpp             = gpp,   # remove this again later
        # ca              = ca,
        # gammastar       = gammastar,
        # kmm             = kmm,
        # ns_star         = ns_star,
        chi               = chi,
        xi                = xi,
        # mj              = out_optchi$mj,
        # mc              = out_optchi$mc,
        ci                = ci,
        # iwue              = iwue,
        gs                = gs,
        vcmax             = vcmax,
        vcmax25           = vcmax25,
        jmax              = jmax,
        jmax25            = jmax25,
        kphio             = kphio_accl,
        tc_growth_leaf    = tc_growth_leaf
        # rd              = rd
    )
    
    return( out )
}

# ANALYTICAL RPMODEL FUNCTIONS ####
calc_optimal_chi <- function(kmm, gammastar, ns_star, ca, vpd, beta ){
  
  # Input:    - float, 'kmm' : Pa, Michaelis-Menten coeff.
  #           - float, 'ns_star'  : (unitless) viscosity correction factor for water
  #           - float, 'vpd' : Pa, vapor pressure deficit
  # Output:   float, ratio of ci/ca (chi)
  # Features: Returns an estimate of leaf internal to ambient CO2
  #           partial pressure following the "simple formulation".
  # Depends:  - kc
  #           - ns
  #           - vpd
  
  ## Avoid negative VPD (dew conditions), resolves issue #2 (https://github.com/stineb/rpmodel/issues/2)
  vpd <- ifelse(vpd < 0, 0, vpd)
  
  ## leaf-internal-to-ambient CO2 partial pressure (ci/ca) ratio
  xi  <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )
  chi <- gammastar / ca + ( 1.0 - gammastar / ca ) * xi / ( xi + sqrt(vpd) )
  
  ## more sensible to use chi for calculating mj - equivalent to code below
  # # Define variable substitutes:
  # vdcg <- ca - gammastar
  # vacg <- ca + 2.0 * gammastar
  # vbkg <- beta * (kmm + gammastar)
  #
  # # Check for negatives, vectorized
  # mj <- ifelse(ns_star>0 & vpd>0 & vbkg>0,
  #              mj(ns_star, vpd, vacg, vbkg, vdcg, gammastar),
  #              rep(NA, max(length(vpd), length(ca)))
  #              )
  
  ## alternative variables
  gamma <- gammastar / ca
  kappa <- kmm / ca
  
  ## use chi for calculating mj
  mj <- (chi - gamma) / (chi + 2 * gamma)
  
  ## mc
  mc <- (chi - gamma) / (chi + kappa)
  
  ## mj:mv
  mjoc <- (chi + kappa) / (chi + 2 * gamma)
  
  # format output list
  out <- list(
    xi = xi,
    chi = chi,
    mc = mc,
    mj = mj,
    mjoc = mjoc
  )
  return(out)
}


calc_mprime <- function( mc ){
  # Input:  mc   (unitless): factor determining LUE
  # Output: mpi (unitless): modified m accounting for the co-limitation
  #                         hypothesis after Prentice et al. (2014)
  
  kc <- 0.41          # Jmax cost coefficient
  
  mpi <- mc^2 - kc^(2.0/3.0) * (mc^(4.0/3.0))
  
  # Check for negatives:
  mpi <- ifelse(mpi>0, sqrt(mpi), NA)
  
  return(mpi)
}


calc_lue_vcmax_wang17 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
  
  ## Include effect of Jmax limitation
  len <- length(out_optchi[[1]])
  mprime <- calc_mprime( out_optchi$mj )
  
  out <- list(
    
    mprime = mprime,
    
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * mprime * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * mprime / out_optchi$mj * soilmstress,
    
    ## complement for non-smith19
    omega      = rep(NA, len),
    omega_star = rep(NA, len)
    
  )
  
  return(out)
}



calc_lue_vcmax_smith19 <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
  
  len <- length(out_optchi[[1]])
  
  # Adopted from Nick Smith's code:
  # Calculate omega, see Smith et al., 2019 Ecology Letters
  omega <- function( theta, c_cost, m ){
    
    cm <- 4 * c_cost / m                        # simplification term for omega calculation
    v  <- 1/(cm * (1 - theta * cm)) - 4 * theta # simplification term for omega calculation
    
    # account for non-linearities at low m values
    capP <- (((1/1.4) - 0.7)^2 / (1-theta)) + 3.4
    aquad <- -1
    bquad <- capP
    cquad <- -(capP * theta)
    m_star <- (4 * c_cost) / polyroot(c(aquad, bquad, cquad))
    
    omega <- ifelse(  m < Re(m_star[1]),
                      -( 1 - (2 * theta) ) - sqrt( (1 - theta) * v),
                      -( 1 - (2 * theta))  + sqrt( (1 - theta) * v)
    )
    return(omega)
  }
  
  ## constants
  theta <- 0.85    # should be calibratable?
  c_cost <- 0.05336251
  
  
  ## factors derived as in Smith et al., 2019
  omega <- omega( theta = theta, c_cost = c_cost, m = out_optchi$mj )          # Eq. S4
  omega_star <- 1.0 + omega - sqrt( (1.0 + omega)^2 - (4.0 * theta * omega) )       # Eq. 18
  
  ## Effect of Jmax limitation
  mprime <- out_optchi$mj * omega_star / (8.0 * theta)
  
  ## Light use efficiency (gpp per unit absorbed light)
  lue <- kphio * ftemp_kphio * mprime * c_molmass * soilmstress
  
  # calculate Vcmax per unit aborbed light
  vcmax_unitiabs  <- kphio * ftemp_kphio * out_optchi$mjoc * omega_star / (8.0 * theta) * soilmstress   # Eq. 19
  
  out <- list(
    lue            = lue,
    vcmax_unitiabs = vcmax_unitiabs,
    omega          = omega,
    omega_star     = omega_star
  )
  
  return(out)
}



calc_lue_vcmax_none <- function(out_optchi, kphio, ftemp_kphio, c_molmass, soilmstress){
  ## Do not include effect of Jmax limitation
  len <- length(out_optchi[[1]])
  
  out <- list(
    
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * out_optchi$mj * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * out_optchi$mjoc * soilmstress,
    
    ## complement for non-smith19
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
  )
  
  return(out)
}



calc_lue_vcmax_c4 <- function( kphio, ftemp_kphio, c_molmass, soilmstress ){
  
  len <- length(kphio)
  out <- list(
    ## Light use efficiency (gpp per unit absorbed light)
    lue = kphio * ftemp_kphio * c_molmass * soilmstress,
    
    ## Vcmax normalised per unit absorbed PPFD (assuming iabs=1), with Jmax limitation
    vcmax_unitiabs = kphio * ftemp_kphio * soilmstress,
    
    ## complement for non-smith19
    omega               = rep(NA, len),
    omega_star          = rep(NA, len)
  )
  
  return(out)
}


calc_chi_c4 <- function(){
  
  # (Dummy-) ci:ca for C4 photosynthesis
  out <- list( chi=9999, mc=1, mj=1, mjoc=1 )
  return(out)
}


# INSTANT RPMODEL FUNCTIONS ####
run_rpmodel_accl <- function(settings = NA,     # Options: Setting to NA takes default settings from get_settings()
                             df_drivers,        # Options: Input of drivers for sites, needs: sitename (chr) and forcing (nested tibble, see rsofun v3.3 setup)
                             df_evaluation,     # Options: Input of evaluation sitename, date and target value: sitename | date | target_var
                             target_var = "toremve"){       # Options: Character specifying target_var
  
  
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
      kphio = NA,
      tc_growth_leaf = NA) %>% 
    right_join(df_evaluation) %>% 
    drop_na(ppfd) # To make sure, no rows with NA's enter loop below
  
  ## Calculating acclimated variables:
  for (row in 1:nrow(df_rpmodel_accl)) {
    message('\014')
    message("rpmodel_accl: ", round(row / nrow(df_rpmodel_accl) * 100), " %")
    
    df_loop <- rpmodel(
      
      ## Acclimated inputs
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
      method_eb      = settings$rpmodel_accl$method_eb
      )
    
    ## Definition of rpmodel output
    df_rpmodel_accl$chi[row]     <- df_loop$chi
    df_rpmodel_accl$ci[row]      <- df_loop$ci
    df_rpmodel_accl$xi[row]      <- df_loop$xi
    df_rpmodel_accl$gs[row]      <- df_loop$gs
    df_rpmodel_accl$vcmax[row]   <- df_loop$vcmax
    df_rpmodel_accl$vcmax25[row] <- df_loop$vcmax25
    df_rpmodel_accl$jmax[row]    <- df_loop$jmax
    df_rpmodel_accl$jmax25[row]  <- df_loop$jmax25
    df_rpmodel_accl$kphio[row]   <- df_loop$kphio
    df_rpmodel_accl$tc_growth_leaf[row]   <- df_loop$tc_growth_leaf
  }
  
  ## Return final dataframe, ready for instantaneous P-Model
  out <- df_rpmodel_accl %>%
    ## Get nested acclimated data
    dplyr::select(sitename, date, chi, ci, xi, gs, vcmax, vcmax25, jmax, jmax25, kphio, tc_growth_leaf) %>% 
    nest(rpm_accl = c(-sitename, -date, -tc_growth_leaf)) %>% 
    
    ## Get nested forcing data (note: tc_home is saved in siteinfo)
    left_join(df_drivers %>%
                dplyr::select(sitename, forcing, siteinfo) %>% 
                unnest(c(forcing, siteinfo)) %>% 
                dplyr::select(any_of(df_drivers$forcing[[1]] %>% names()), "tc_home", "sitename")) %>%
    
    ## Create nested dataframes
    nest(forcing = c(df_drivers$forcing[[1]] %>% dplyr::select(-date) %>% names(), "tc_home", "tc_growth_leaf")) %>% 
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
  } else if (is.numeric(settings$rpmodel_exp$ppfd)) {
    df_sim$ppfd <- settings$rpmodel_exp$ppfd
  }
  
  if (settings$rpmodel_inst$method_vcmax25 == "prescribed") {
    df_sim$vcmax25 <- df_sim$vcmax25_pft
    
    # Acclimated jmax_to_vcmax taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment:
    df_sim$jmax25  <- df_sim$vcmax25 * (2.56 - (0.0375 * df_sim$tc_home) + (-0.0202 * (df_sim$tc_growth_air - df_sim$tc_home)))
  }
  
  ## Get starting time of computation
  # start.time <- Sys.time()
  
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
  gs <- anet/(ca - ci)
  
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

# HELPER FUNCTIONS ####
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
    rpmodel_exp = list(ppfd               = "metainfo")       # Options: "ambient", metainfo" or numeric in mol/m2/s
  )
  
  return(settings)
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
