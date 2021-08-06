# ....................................
# SCRIPT FOR LEAF ENERGY BALANCE MODEL
# ....................................

# LeafEnergyBalance() taken from plantecophys package
LeafEnergyBalance <- function(Tleaf = 21.5, Tair = 20,
                              gs = 0.15,
                              PPFD = 1500, VPD = 2, Patm = 101,
                              Wind = 2, Wleaf = 0.02,
                              StomatalRatio = 1,   # 2 for amphistomatous
                              LeafAbs = 0.5,   # in shortwave range, much less than PAR
                              returnwhat=c("balance","fluxes")){
    
    # EDIT: Functions needed in here and not available directly from Plantecophys
    Tk <- function(x){x+273.15}
    
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
        
        # out <- list(tc_leaf       = Tleaf,     # NEW
        # 						tc_leaf_star  = Tleaf2,    # NEW
        # 						eps           = EnergyBal) # NEW
        # return(out)                            # NEW
    }
    
    if(returnwhat == "fluxes"){
        
        l <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)
        return(l)
    }
}

# .............................................................................

calc_stomatal_conductance <- function(tc_leaf   = 25,
                                      tc_growth = 25,
                                      tc_home   = 25,
                                      ppfd      = 1500 / 10^-6,
                                      co2       = 400,
                                      patm      = 101325,
                                      vpd       = 5000,
                                      xi        = 50,
                                      vcmax25   = 5e-5,
                                      jmax25    = 10e-5,
                                      kphio     = 0.10,
                                      toggles   = NA
){
    
    
    if (length(toggles) == 1){
        toggles <- list(XI_INST  = F,
                        FIXED_CI = F,
                        KUMA_PAR = T,
                        TRIAL_AJ = F)
    }
    
    
    # Beta (can also be piped for adjustable input)
    beta <- 146.0
    
    # Ambient CO2
    ca = co2_to_ca(co2, patm)
    
    # Get gammastar
    gammastar = calc_gammastar(tc_leaf, patm)
    
    # Get michaelis-menten constant
    kmm = calc_kmm(tc_leaf, patm)
    
    # Get temperature responses of vcmax and jmax
    vcmax  = vcmax25 * calc_ftemp_inst_vcmax(tc_leaf, tc_growth,          toggles$KUMA_PAR )
    jmax   = jmax25  * calc_ftemp_inst_jmax( tc_leaf, tc_growth, tc_home, toggles$KUMA_PAR )
    
    # Get viscosity
    ns      <- calc_viscosity_h2o( tc_leaf, patm )            # Pa s
    ns25    <- calc_viscosity_h2o( tc = 25.0, p = 101325.0 )  # Pa s
    ns_star <- ns / ns25                                      # (unitless)
    
    # TODO: Do we use instantaneous or acclimated xi?
    # Get xi
    if (toggles$XI_INST) {xi <- sqrt( (beta * ( kmm + gammastar ) ) / ( 1.6 * ns_star ) )}
    
    # Instant ci
    kv		= (ca - gammastar) / (1 + xi / sqrt(vpd))
    ci_j  = ca - kv
    ci_c  = ca - kv
    
    if (toggles$FIXED_CI) {
        ci_j  = 27.5
        ci_c  = 27.5
    }
    
    # Electron rate
    L     = 1.0 / sqrt(1.0 + ((4.0 * kphio * ppfd) / jmax)^2)
    a_j 	= L * kphio * ppfd * (ci_j - gammastar)/(ci_j + 2 * gammastar)
    gs_j	= a_j / kv
    
    if (toggles$TRIAL_AJ) {
        # Testing hyperbola formulation
        theta   <- 0.85
        j       <- (kphio * ppfd + jmax - sqrt(( kphio * ppfd + jmax)^2 - (4*kphio*theta*ppfd*jmax))) / (2*theta)
        m_j     <- ci_j / (ci_j + 2*gammastar)
        a_j     <- j/4 * m_j * (1.0 - gammastar/ci_j)
        gs_j	  <- a_j / kv
    }
    
    # Rubisco rate
    a_c = vcmax * (ci_c - gammastar)/(ci_c + kmm) * (1.0 - gammastar/ci_j)
    gs_c = a_c / kv
    
    # TODO: Does this make sense?
    # Get limiting gs
    gs <- min(gs_c, gs_j)
    
    return(gs)
}


# .............................................................................
calc_tc_leaf_from_tc_air <- function(tc_air   = 25,
                                     gs       = 0.15 / 101325,
                                     vpd      = 2 * 10^3,
                                     patm     = 101325,
                                     ppfd     = 1500 / 10^6,
                                     wind     = 2,
                                     wleaf    = 0.02,
                                     stoma_r  = 1,
                                     leaf_abs = 0.5
){
    
    # LeafEnergyBalance is equivalent to "maximize_this_tc_leaf" with its Tleaf getting optimized
    
    sol_optimr <-	optimr::optimr(
        # Parameter boundaries to optimize within:
        par       = 15,
        lower     = 15 - tc_air, # OLD: 0 
        upper     = 15 + tc_air, # OLD: 40 
        # Function to optimize and its inputs:
        fn        = LeafEnergyBalance,
        Tair      = tc_air,
        gs        = gs   * patm,        # input has to be in mol/m2/s
        VPD       = vpd  / 10^3,        # input has to be in kPa
        Patm      = patm / 10^3,        # input has to be in kPa
        PPFD      = ppfd / 10^6,        # input has to be in umol
        Wind      = wind,
        Wleaf     = wleaf,
        StomatalRatio = stoma_r,        # 2 for amphistomatous
        LeafAbs   = leaf_abs,
        # Optimr settings:
        method    = "L-BFGS-B",
        control   = list( maxit = 100, maximize = TRUE ))
    
    return(sol_optimr$par)
}

maximize_this_tc_leaf <- function(tc_leaf   = 25, # This gets optimized in optimr()
                                  tc_air    = 25,
                                  tc_growth = 25,
                                  tc_home   = 25,
                                  ppfd      = 1500 / 10^6,
                                  co2       = 400,
                                  patm      = 101325,
                                  vpd       = 5000,
                                  xi        = 50,
                                  vcmax25   = 5e-5,
                                  jmax25    = 10e-5,
                                  kphio     = 0.10,
                                  toggles   = NA,
                                  wind      = 2,
                                  wleaf     = 0.02,
                                  stoma_r   = 1,
                                  leaf_abs  = 0.5
){
    
    
    if (T) {
        #   T = get tc_leaf based on plantecophys LEB
        #   F = get tc_leaf based on tealeaves LEB
        tc_leaf_x <- calc_tc_leaf_from_tc_air(tc_air = tc_air,
                                              gs = calc_stomatal_conductance(tc_leaf   = tc_leaf, 
                                                                             tc_growth = tc_growth,
                                                                             tc_home   = tc_home,
                                                                             ppfd      = ppfd,
                                                                             co2       = co2,
                                                                             patm      = patm,
                                                                             vpd       = vpd,
                                                                             xi        = xi,
                                                                             vcmax25   = vcmax25,
                                                                             jmax25    = jmax25,
                                                                             kphio     = kphio,
                                                                             toggles   = toggles),
                                              wind     = wind,
                                              wleaf    = wleaf,
                                              stoma_r  = stoma_r,
                                              leaf_abs = leaf_abs)
    } else {
        gs <- calc_stomatal_conductance(tc_leaf   = tc_leaf, 
                                        tc_growth = tc_growth,
                                        tc_home   = tc_home,
                                        ppfd      = ppfd,
                                        co2       = co2,
                                        patm      = patm,
                                        vpd       = vpd,
                                        xi        = xi,
                                        vcmax25   = vcmax25,
                                        jmax25    = jmax25,
                                        kphio     = kphio,
                                        toggles   = toggles)
        
        # Get leaf parameters
        leaf_par <- make_leafpar(
            replace = list(
                g_sw = set_units(gs, "mol/m^2/s/Pa")))
        
        # Get environmental parameters Next, we'll change the air temperature to 25 degree C (= 298.15 K)
        enviro_par <- make_enviropar(
            replace = list(
                T_air = set_units(tc_air + 273.15, "K")))
        # Get physical constants
        constants  <- make_constants()
        
        # Get tc_leaf
        tc_leaf_x <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)$T_leaf %>% 
            set_units("degree_Celsius") %>% drop_units()
    }
    
    
    # Get difference between tc_leaf and tc_leaf_x
    eps <- (tc_leaf_x - tc_leaf)^2
    
    return(eps)
}

# .............................................................................
calc_tc_leaf_final <- function(tc_air    = 25,
                               tc_growth = 25,
                               tc_home   = 25,
                               ppfd      = 130,
                               fapar     = 1,
                               patm      = 101325,
                               co2       = 400,
                               vpd       = 1000,
                               kphio     = 0.05,
                               toggles   = NA,
                               wind      = 2,
                               wleaf     = 0.02,
                               stoma_r   = 1,
                               leaf_abs  = 0.5,
                               method_jmaxlim_inst = "smith37",
                               beta      = 146, 
                               c_cost    = 0.41
){
    
    # Get optimized tc_leaf
    
    # Call optimize()
    
    
    sol_optimize <- tryCatch(
        {
            sol_optimize <- optimize(maximize_this_tc_leaf,
                                     interval  = c(1, 40),
                                     tc_air    = tc_air,
                                     tc_growth = tc_growth,
                                     tc_home   = tc_home,
                                     ppfd      = ppfd,
                                     fapar     = fapar,
                                     co2       = co2,
                                     patm      = patm,
                                     vpd       = vpd,
                                     kphio		 = kphio,
                                     toggles   = toggles,
                                     wind      = wind,
                                     wleaf     = wleaf,
                                     stoma_r   = stoma_r,
                                     leaf_abs  = leaf_abs,
                                     method_jmaxlim_inst = method_jmaxlim_inst,
                                     beta      = beta,
                                     c_cost    = c_cost)
        },
        warning = function(cond){
            # print("There was a warning")
            return(NA)
        },
        error = function(cond){
            # print("This message will not be printed.")
            return(NA)
        },
        finally = {
            #pass
        })
    
    if (length(sol_optimize) == 1) {
        return (tc_air)
    }
    
    return(sol_optimize$minimum)
}
