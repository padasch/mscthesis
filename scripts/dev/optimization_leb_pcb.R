# Optimization of leaf energy balance and plant cost balance ####
# Optimal gs, vcmax, jmax ####

calc_optimal_gs_vcmax_jmax <- function(par,
                                       tc_leaf,
                                       tc_air,
                                       ppfd,
                                       fapar = 1,
                                       co2,
                                       patm,
                                       vpd,
                                       kphio,
                                       beta = 146.0,
                                       method_jmaxlim_inst,
                                       maximize = FALSE,
                                       return_all = FALSE,
                                       units_out  = "per-s") {
        
    
    ## 1: Parameters to be optimized:
    vcmax <- par[1]
    jmax  <- par[2]
    gs    <- par[3]
    
    ## 2: Get photosynthetic variables based on environmental conditions:
    kmm       <- calc_kmm(tc_leaf, patm)
    gammastar <- calc_gammastar(tc_leaf, patm)
    ns_star   <- calc_viscosity_h2o(tc_leaf, patm) / calc_viscosity_h2o(25, 101325)
    ca        <- co2_to_ca(co2, patm)
    kphio     <- kphio * calc_ftemp_kphio( tc_leaf, c4 = F)
    iabs      <- ppfd * fapar
    
    
    ## 3: Calculate assimilation rates with to-be-optimized jmax, vcmax and gs:
    
    ## 3.1: Electron transport is limiting
    ## Solve quadratic equation system using: A(Fick's Law) = A(Jmax Limitation)
    ## This leads to a quadratic equation:
    ## A * ci^2 + B * ci + C  = 0
    ## 0 = a + b*x + c*x^2
    
    ## Jmax Limitation following Smith (1937):
    if (method_jmaxlim_inst == "smith37") {
        ## A = gs * (ca - ci)
        ## A = kphio * iabs (ci-gammastar)/ci+2*gammastar) * L
        ## L = 1 / sqrt(1 + ((4 * kphio * iabs)/jmax)^2)
        
        ## with
        L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
        A <- -gs
        B <- gs * ca - 2 * gammastar * gs - L * kphio * iabs
        C <- 2 * gammastar * gs * ca + L * kphio * iabs * gammastar
        
        ci_j <- QUADM(A, B, C)
        a_j  <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L  
        
        c_cost <- 0.103 # As estimated by Wang et al. (2017)
    }
    
    ## Jmax Limitation following Farquhar (1989):
    if (method_jmaxlim_inst == "farquhar89") {
        ## A = gs * (ca - ci)
        ## A = j/4 * (ci-gammastar)/ci+2*gammastar)
        ## j = (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4 * kphio * theta * iabs * jmax))) / (2*theta)
        
        ## with
        theta <- 0.85
        j <- (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4 * kphio * theta * iabs * jmax))) / (2 * theta)
        A <- -gs
        B <- gs * ca - 2 * gammastar * gs - j/4
        C <- 2 * gammastar * gs * ca + gammastar * j/4
        
        ci_j <- ci_j <- QUADM(A, B, C)
        a_j <- j/4 * (ci_j - gammastar)/(ci_j + 2 * gammastar)
        
        c_cost <- 0.053 # As estimated by Smith et al. (2019)
    }
    
    ## 4: Rubisco is limiting
    ## Solve Eq. system
    ## A = gs (ca- ci)
    ## A = Vcmax * (ci - gammastar)/(ci + Kmm)
    
    ## This leads to a quadratic equation:
    ## A * ci^2 + B * ci + C  = 0
    ## 0 = a + b*x + c*x^2
    
    ## with
    A <- -1.0 * gs
    B <- gs * ca - gs * kmm - vcmax
    C <- gs * ca * kmm + vcmax * gammastar
    
    ci_c <- QUADM(A, B, C)
    a_c <- vcmax * (ci_c - gammastar) / (ci_c + kmm)
    
    ## 5. Take minimum of the two assimilation rates and maximum of the two ci
    ci <- max(ci_c, ci_j)
    a_gross <- min( a_j, a_c ) # Original approach using min()
    # a_gross <- -QUADP(A = 1 - 1E-07, B = a_c + a_j, C = a_c*a_j) # Altenrative approach using hyperbolic minumum to avoid discontinuity (see Duursma et al (2015), Eq. (5))
    
    ## 6. Calculate costs as defined by the criteria: (1.6*ns_star*gs*vpd + beta*vcmax)/Ac + jmax*c/Aj = min! # TODO -> fix formulation
    cost_transp <- 1.6 * ns_star * gs * vpd
    cost_vcmax  <- beta * vcmax
    cost_jmax   <- c_cost * jmax
    
    net_assim <- -(cost_transp + cost_vcmax + cost_jmax) / a_gross # Original
    
    # net_assim <- -(cost_transp/ns_star + cost_vcmax + cost_jmax) / a_gross
    # net_assim <- -((cost_transp/ns_star + cost_vcmax) / a_c  + cost_jmax/a_j) # TODO -> try out!
    
    if (maximize) net_assim <- -net_assim
    
    if (return_all){
        
        ## Turn per-day units back into per-second
        if (units_out == "per-s") {
            vcmax <- vcmax         / (3600 * 24) # Final unit: [mol/m2/s]
            jmax <- jmax           / (3600 * 24) # Final unit: [mol/m2/s]
            gs <- gs               / (3600 * 24) # Final unit: [mol/m2/s/Pa]
            a_c <- a_c             / (3600 * 24) # Final unit: [mol/m2/s]
            a_j <- a_j             / (3600 * 24) # Final unit: [mol/m2/s]
            a_gross <- a_gross     / (3600 * 24) # Final unit: [mol/m2/s]
            net_assim <- net_assim / (3600 * 24) # Final unit: [-]
        }
        
        ## Output
        return(
            tibble(
                vcmax_mine = vcmax,
                jmax_mine = jmax,
                gs_mine = gs,
                ci_mine = ci,
                chi_mine = ci / ca,
                a_c_mine = a_c,
                a_j_mine = a_j,
                a_gross = a_gross,
                ci_c_mine = ci_c,
                ci_j_mine = ci_j,
                cost_transp = cost_transp,
                cost_vcmax = cost_vcmax,
                cost_jmax = cost_jmax,
                net_assim = net_assim,
                method_jmaxlim_inst = method_jmaxlim_inst
            )
        )
    } else {
        return( net_assim )
    }
}
# .............................................................................
minimize_costs_gs_vcmax_jmax <- function(tc_leaf,
                                         patm,
                                         co2,
                                         vpd,
                                         ppfd,
                                         kphio,
                                         method_jmaxlim_inst,
                                         vcmax_start = NA,
                                         jmax_start  = NA,
                                         gs_start    = NA
                                         ) {
    
    ## TODO: Input for optimization has to be in per-day to work properly:
    ppfd        <- ppfd * 3600 * 24
    vcmax_start <- 20
    jmax_start  <- 20
    gs_start    <- 0.3
    
    
    ## Get optimal parameters (TODO: Output order of magnitude depends on lower/upper boundaries OoM)
    out_optim <- optimr::optimr(

        ## Optimization inputs:
        par        = c( vcmax_start,       jmax_start      , gs_start),
        lower      = c( vcmax_start*10^-4, jmax_start*10^-4, gs_start*10^-4 ),
        upper      = c( vcmax_start*10^4,  jmax_start*10^4 , gs_start*10^4  ),
        fn         = calc_optimal_gs_vcmax_jmax,
        method     = "L-BFGS-B",
        control    = list(maxit = 1000),

        ## Function inputs:
        tc_leaf    = tc_leaf,
        patm       = patm,
        co2        = co2,
        vpd        = vpd,
        ppfd       = ppfd,
        kphio      = kphio,
        method_jmaxlim_inst = method_jmaxlim_inst,
        maximize   = TRUE)

    
    optimized_par <- calc_optimal_gs_vcmax_jmax(par        = out_optim$par,
                                                tc_leaf    = tc_leaf,
                                                patm       = patm,
                                                co2        = co2,
                                                vpd        = vpd,
                                                ppfd       = ppfd,
                                                kphio      = kphio,
                                                method_jmaxlim_inst = method_jmaxlim_inst,
                                                units_out  = "per-s",
                                                return_all = T)
    
    return(optimized_par)
    
}

# .............................................................................
# Optimal leaf temperature ####
# .............................................................................
LeafEnergyBalance <- function(Tleaf, # This parameter is optimized
                              Tair,
                              gs,
                              PPFD,
                              VPD,
                              Patm,
                              Wind = 2,
                              Wleaf = 0.02,
                              StomatalRatio = 1, # 2 for amphistomatous
                              LeafAbs = 0.5, # in shortwave range, much less than PAR
                              returnwhat = c("balance", "fluxes")
) {
    
    ## Ouput: Difference between an assumed tc_leaf and the tc_leaf calculated from tc_air and gs
    
    ## Function needs the following units:
    # Tleaf         degC
    # Tair          degC
    # gs            mol/m2/s, stomatal conductance to H2O!
    # PPFD          umol/m2/s
    # VPD           kPa
    # Patm          kPa
    # Wind          m/s
    # Wleaf         m
    # StomatalRatio unitless (1 or 2)
    # LeafAbs       [0, 1]
    
    ## Input is in SI units, adjustments needed:
    gs   <- gs * (esat(Tair, patm / 1000) - VPD)  # mol H20/m2/s/Pa to mol H2O/m2/s by multiplying with water vapor pressure
    PPFD <- PPFD * 10^6  # mol/m2/s to umol/m2/s
    Patm <- Patm / 1000  # Pa to kPa
    VPD  <- VPD  / 1000  # Pa to kPa
    
    # Output
    returnwhat <- match.arg(returnwhat)
    
    # Constants
    Boltz <- 5.67 * 10^-8     # w M-2 K-4
    Emissivity <- 0.95        # -
    LatEvap <- 2.54           # MJ kg-1
    CPAIR <- 1010.0           # J kg-1 K-1
    
    H2OLV0 <- 2.501e6         # J kg-1
    H2OMW <- 18e-3            # J kg-1
    AIRMA <- 29.e-3           # mol mass air (kg/mol)
    AIRDENS <- 1.204          # kg m-3
    UMOLPERJ <- 4.57
    DHEAT <- 21.5e-6          # molecular diffusivity for heat
    
    # Density of dry air
    AIRDENS <- Patm*1000/(287.058 * Tk(Tair))
    
    # Latent heat of water vapour at air temperature (J mol-1)
    LHV <- (H2OLV0 - 2.365E3 * Tair) * H2OMW
    
    # Const s in Penman-Monteith equation  (Pa K-1)
    SLOPE <- (esat(Tair + 0.1) - esat(Tair)) / 0.1
    
    # Radiation conductance (mol m-2 s-1)
    Gradiation <- 4.*Boltz*Tk(Tair)^3 * Emissivity / (CPAIR * AIRMA)
    
    # See Leuning et al (1995) PC&E 18:1183-1200 Appendix E
    # Boundary layer conductance for heat - single sided, forced convection
    CMOLAR <- Patm*1000 / (8.314 * Tk(Tair))   # .Rgas() in package...
    Gbhforced <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR
    
    # Free convection
    GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
    Gbhfree <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR
    
    # Total conductance to heat (both leaf sides)
    Gbh <- 2*(Gbhfree + Gbhforced)
    
    # Heat and radiative conductance
    Gbhr <- Gbh + 2*Gradiation
    
    # Boundary layer conductance for water (mol m-2 s-1)
    Gbw <- StomatalRatio * 1.075 * Gbh  # Leuning 1995
    gw <- gs*Gbw/(gs + Gbw)
    
    # Longwave radiation
    # (positive flux is heat loss from leaf)
    Rlongup <- Emissivity*Boltz*Tk(Tleaf)^4
    
    # Rnet
    Rsol <- 2*PPFD/UMOLPERJ   # W m-2
    Rnet <- LeafAbs*Rsol - Rlongup   # full
    
    # Isothermal net radiation (Leuning et al. 1995, Appendix)
    ea <- esat(Tair) - 1000*VPD
    ema <- 0.642*(ea/Tk(Tair))^(1/7)
    Rnetiso <- LeafAbs*Rsol - (1 - ema)*Boltz*Tk(Tair)^4 # isothermal net radiation
    
    # Isothermal version of the Penmon-Monteith equation
    GAMMA <- CPAIR*AIRMA*Patm*1000/LHV
    ET <- (1/LHV) * (SLOPE * Rnetiso + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbhr/gw)
    
    # Latent heat loss
    lambdaET <- LHV * ET
    
    # Heat flux calculated using Gradiation (Leuning 1995, Eq. 11)
    Y <- 1/(1 + Gradiation/Gbh)
    H2 <- Y*(Rnetiso - lambdaET)
    
    # Heat flux calculated from leaf-air T difference.
    # (positive flux is heat loss from leaf)
    H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf)
    
    # Leaf-air temperature difference recalculated from energy balance.
    # (same equation as above!)
    Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR))
    
    # Difference between input Tleaf and calculated, this will be minimized.
    # EnergyBal <- (Tleaf - Tleaf2)           # OLD, needed to work with uniroot()
    EnergyBal <- (Tleaf - Tleaf2)^2         # NEW, needed to work with optimr()
    # EnergyBal <- abs(Tleaf - Tleaf2)        # NEW, needs more iterations than ()^2
    
    if(returnwhat == "balance"){
        
        return(EnergyBal)                      # OLD

    }
    
    if(returnwhat == "fluxes"){
        
        l <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)
        return(l)
    }
}

# .............................................................................
diff_tcleaf_in_and_tcleaf_eb <- function(tc_leaf, # This parameter is optimized
                                         tc_air,
                                         ppfd,
                                         co2,
                                         patm,
                                         vpd,
                                         kphio,
                                         method_jmaxlim_inst,
                                         method_eb = "plantecophys") {
    

                                           
    ## Output: difference in tc_leaf assumed for gs and tc_leaf from energy balance
    
    ## 1: Get optimal gs, vcmax and jmax at given tc_leaf
    optim_vars <- minimize_costs_gs_vcmax_jmax(tc_leaf,
                                               patm,
                                               co2,
                                               vpd,
                                               ppfd,
                                               kphio,
                                               method_jmaxlim_inst)
    
    ## 2. Get optimal gs for water for calculating tc_leaf via energy balance:
    gs_water <- optim_vars$gs_mine / 1.6
    
    if (method_eb == "plantecophys") {
        ## 2.1: Via plantecophys energy balance
        
        tc_leaf_leb <- optimr::optimr(
                                       # Parameter boundaries to optimize within:
                                       par       = 15,
                                       lower     = max(1, 15 - tc_air), # OLD: 0
                                       upper     = 15 + tc_air, # OLD: 40
                                       # Function to optimize and its inputs:
                                       fn        = LeafEnergyBalance,
                                       Tair      = tc_air,
                                       gs        = gs_water,
                                       VPD       = vpd,
                                       Patm      = patm,
                                       PPFD      = ppfd,
                                       # Optimr settings:
                                       method    = "L-BFGS-B",
                                       control   = list(maxit = 10000, maximize = TRUE))$par
        
    } else if (method_eb == "tealeaves") {
        ## 2.2: Via tealeaves energy balance
        
        # Get relative humidity from vpd:
        RH <- (100 - ((vpd) / esat(tc_air, patm/1000))) / 100
        
        # Get leaf parameters:
        leaf_par <- make_leafpar(
            replace = list(
                g_sw = set_units(gs_water, "mol/m^2/s/Pa")))
        
        # Get environmental parameters:
        enviro_par <- make_enviropar(
            replace = list(
                T_air = set_units(tc_air + 273.15, "K"),
                RH    = as_units(RH),
                P     = set_units(patm, "Pa")))
        
        # Get physical constants:
        constants  <- make_constants()
        
        # Get tc_leaf:
        tc_leaf_leb <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)$T_leaf %>% 
            set_units("degree_Celsius") %>% drop_units()
    }
    
    # 3: Get squared difference between from assumed tc_leaf and actual tc_leaf
    eps <- (tc_leaf_leb - tc_leaf)^2
    
    return(eps)
}
# .............................................................................
calc_tc_leaf_from_tc_air <- function(tc_air,
                                     ppfd,
                                     patm,
                                     co2,
                                     vpd,
                                     kphio,
                                     method_jmaxlim_inst,
                                     method_eb) {
      
    
    ## IN DEV:
    sol_optimize <- optimize(diff_tcleaf_in_and_tcleaf_eb,
                             interval  = c(max(tc_air - 15, 1), tc_air + 15),
                             tc_air    = tc_air,
                             ppfd      = ppfd,
                             co2       = co2,
                             patm      = patm,
                             vpd       = vpd,
                             kphio		 = kphio,
                             method_jmaxlim_inst = method_jmaxlim_inst,
                             method_eb = method_eb)
    
    return(sol_optimize$minimum)
    
  ## TODO: OPTIM() AND OPTIMR() BELOW CRASH FOR SOME REASON...
    
    # out_optim <- optimr::optimr(
    #     
    #     ## Optimization inputs:
    #     par        = tc_air,
    #     lower      = 1,
    #     upper      = 40,
    #     fn         = diff_tcleaf_in_and_tcleaf_eb,
    #     method     = "L-BFGS-B",
    #     control    = list(maxit = 1000),
    #     
    #     ## Function inputs:
    #     tc_air = tc_air,
    #     ppfd = ppfd,
    #     patm = patm,
    #     co2 = co2,
    #     vpd = vpd,
    #     kphio = kphio,
    #     method_jmaxlim_inst = method_jmaxlim_inst,
    #     method_eb = method_eb)
    #     
    # out_optim <- optim(
    # 
    #     ## Optimization inputs:
    #     par        = tc_air,
    #     lower      = 1,
    #     upper      = 40,
    #     fn         = diff_tcleaf_in_and_tcleaf_eb,
    #     method     = "L-BFGS-B",
    #     control    = list(maxit = 1000),
    # 
    #     ## Function inputs:
    #     tc_air = tc_air,
    #     ppfd = ppfd,
    #     patm = patm,
    #     co2 = co2,
    #     vpd = vpd,
    #     kphio = kphio,
    #     method_jmaxlim_inst = method_jmaxlim_inst,
    #     method_eb = method_eb)
    
    
    
    
    }
    


# .............................................................................
# Tests ####
if (F) {
    
    ## Get dummy variables:
    tc_air <- 25 # degC
    tc_leaf <- 25 # degC
    patm <- 101325 # Pa
    co2 <- 400 # ppm
    vpd <- 2000 # Pa
    ppfd <- 500e-6 # mol/m2/s
    vcmax <- 50e-6 # mol/m2/s
    jmax <- 100E-6 # mol/m2/s
    gs <- 30e-4 # mol CO2 /m2/s/Pa
    fapar <- 1 # -
    kphio <- 0.25 # -
    beta <- 146 # -
    c_cost <- 0.103 # -
    method_jmaxlim_inst <- "smith37"
    method_eb <- "plantecophys"
    method_eb <- "tealeaves"
    
    ## Get actual data:
    df_test <- readRDS("~/data/mscthesis/final/df_drivers_p21.rds") %>% 
      leftjoin()
    dplyr::select(forcing) %>% unnest(forcing) %>% mutate(vcmax = NA, jmax = NA, gs = NA)
    
    ## Checking optimal gs, vcmax, jmax
    for (r in 1:nrow(df_test)) {
      opt <- minimize_costs_gs_vcmax_jmax(tc_leaf = temp,
                                              patm,
                                              co2,
                                              vpd,
                                              ppfd,
                                              kphio,
                                              method_jmaxlim_inst)
      )
      
    }
    
    
    
    ## Function calls
    calc_tcleaf_from_eb(tc_air, gs*patm, vpd, patm, ppfd)
    
    calc_optimal_gs_vcmax_jmax(
        par        = c(vcmax, jmax, gs),
        tc_leaf    = tc_leaf,
        patm       = patm,
        co2        = co2,
        vpd        = vpd,
        ppfd       = ppfd,
        kphio      = kphio,
        method_jmaxlim_inst = method_jmaxlim_inst,
        maximize   = TRUE,
        return_all = T
    )
    
    minimize_costs_gs_vcmax_jmax(tc_leaf,
                                 patm,
                                 co2,
                                 vpd,
                                 ppfd,
                                 kphio,
                                 method_jmaxlim_inst)$gs_mine
    
    
    diff_tcleaf_in_and_tcleaf_eb(tc_air,
                                 ppfd,
                                 patm,
                                 co2,
                                 vpd,
                                 kphio,
                                 method_jmaxlim_inst,
                                 method_eb = "plantecophys")
    
    calc_tc_leaf_from_tc_air(tc_air = 30,
                             ppfd,
                             patm,
                             co2,
                             vpd,
                             kphio,
                             method_jmaxlim_inst,
                             method_eb = "tealeaves")
    
    ## Plots
    df_test <- readRDS("~/data/mscthesis/final/df_drivers_p21.rds") %>%
      dplyr::select(forcing) %>%
      unnest(forcing) %>%
      mutate(tc_leaf = purrr::map_dbl(., ~calc_tc_leaf_from_tc_air(tc_air,
                                                ppfd,
                                                patm,
                                                co2,
                                                vpd,
                                                kphio,
                                                method_jmaxlim_inst,
                                                method_eb = "tealeaves")))
 
}
