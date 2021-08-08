# PHOTOSYNTHESIS FUNCTIONS ####
# .............................................................................
calc_ci <- function(ca, gammastar, xi, vpd, patm, ci_fixed_at_ppm = NA){
    
    if (is.numeric(ci_fixed_at_ppm)) {
        ci <- co2_to_ca(ci_fixed_at_ppm, patm)
        
    } else if (is.na(ci_fixed_at_ppm) | ci_fixed_at_ppm == "prentice14") {
        ci <- ca - (ca - gammastar) / (1 + xi / sqrt(vpd))
    }
    
    return(ci)
}

# .............................................................................
calc_aj <- function(kphio, ppfd, jmax, gammastar, ci, ca, fapar, theta = 0.85, j_method = "smith37", model = "analytical") {

    ## farquhar89: Jmax Limitation following Farquhar (1989):
    ## smith37:    Jmax Limitation following Smith (1937):
    
    ## Analytical model with given ci:
    if (model == "analytical") {
        
        if (j_method == "smith37") {
            j <- (4 * kphio * ppfd) / sqrt(1 + ((4 * kphio * ppfd)/(jmax))^2)
        }
        
        if (j_method == "farquhar89") {
            j <- (kphio * ppfd + jmax - sqrt(( kphio * ppfd + jmax)^2 - (4*kphio*theta*ppfd*jmax))) / (2*theta)
        }
        
        aj    <- j/4 * (1.0 - gammastar/ci) * ci / (ci + 2*gammastar)
        
        out = list(aj = aj,
                   j  =  j)
        
    ## Numerical model without given ci:
    } else if (model == "numerical") {
        
        iabs <- ppfd * fapar
        
        if (j_method == "smith37") {
            ## A = gs (ca - ci)
            ## A = kphio * iabs (ci-gammastar)/ci+2*gammastar) * L
            ## L = 1 / sqrt(1 + ((4 * kphio * iabs)/jmax)^2)
            
            ## with
            L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
            A <- -gs
            B <- gs * ca - 2 * gammastar * gs - L * kphio * iabs
            C <- 2 * gammastar * gs * ca + L * kphio * iabs * gammastar
            
            ci_j <- QUADM(A, B, C)
            j    <- 4 * kphio * iabs * gammastar
            aj  <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L  
        }
        
        if (j_method == "farquhar89") {
            ## A = gs (ca - ci)
            ## A = j/4 * (ci-gammastar)/ci+2*gammastar)
            ## j = (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4 * kphio * theta * iabs * jmax))) / (2*theta)
            
            ## with
            j <- (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4 * kphio * theta * iabs * jmax))) / (2 * theta)
            A <- -gs
            B <- gs * ca - 2 * gammastar * gs - j/4
            C <- 2 * gammastar * gs * ca + gammastar * j/4
            
            ci_j <- ci_j <- QUADM(A, B, C)
            aj <- j/4 * (ci_j - gammastar)/(ci_j + 2 * gammastar)
        }
        
        out = list(aj   = aj,
                   j    =  j,
                   ci_j = ci_j)
        
    }
    
    return(out)
}

# .............................................................................
calc_ac <- function(ci, gammastar, kmm, vcmax, model = "analytical") {
    
    if (model == "analytical") {
        ac <- vcmax * (ci - gammastar)/(ci + kmm)
        
        out <- list(ac = ac)
        
    } else if (model == "numerical") {
        ## Rubisco is limiting
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
        ac <- vcmax * (ci_c - gammastar) / (ci_c + kmm)
        
        out <- list(ac = ac,
                    ci = ci_c)
    }
    
    return(out)
}

# .............................................................................
calc_rd <- function(tc_leaf, vcmax25, tc_growth, q10 = 2, method_rd25 = "atkin15", method_rd_scale = "heskel16") {
    
    ## Get base rate of Rd at 25 degC:
    if (method_rd25 == "atkin15") {
        rd_to_vcmax <- 0.015 # Ratio of Rdark to Vcmax25, Atkin et al., 2015 for C3 herbaceous
    }
    
    if (method_rd25 == "kumarathunge19") {
        rd_to_vcmax <- 0.0360 - 0.0010 * tc_growth # Acclimated rd_to_vcmax taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
    }
    
    ## Get temperature scaling for Rd:
    if (method_rd_scale == "heskel16") {
        apar <- 0.1012
        bpar <- 0.0005
        f <- exp( apar * (tc_leaf - 25.0) - bpar * (tc_leaf^2 - 25.0^2) )
    }
    
    if (method_rd_scale == "arrhenius") {
        f <- calc_ftemp_arrh(tc_leaf + 273.15, dha = 20700) # dha: Activation energy taken from Kumarathunge et al. (2019), Table 1, Mature Natural Environment
    }
    
    if (method_rd_scale == "q10") {
        f <- q10^(tc_leaf - 25.0)
    }
    
    rd <- vcmax25 * rd_to_vcmax * f
    
    return(rd)
}
# .............................................................................
calc_jmax <- function(kphio, iabs, ci, gammastar, method, theta = NA, c = NA){ 
    
    if (method == "smith37" | method == "wang17") {
        c      <- 0.103 # Estimated by Wang et al. (2017)
        c_star <- 4*c
        jmax   <- 4*kphio*iabs* (sqrt(1 / (1 - (c_star*(ci + 2*gammastar)/(ci-gammastar)))^(2/3)) - 1)  
    }
    
    if (method == "farquhar89" | method == "smith19") {
        theta <- 0.85
        c     <- 0.053 # Estimated by Smith et al. (2019)
        m     <- (ci - gammastar) / (ci + 2 * gammastar)
        omega <- -(1 - 2 * theta) + sqrt((1-theta) * (1 / (4*c/m * (1-theta*4*c/m)) - 4*theta))
        jmax  <- kphio * iabs * omega
    }
    
    return(jmax)
}


# .............................................................................
calc_ftemp_arrh <- function(
    tk,
    dha,
    tkref = 298.15
){
    
    # Note that the following forms are equivalent:
    # ftemp = exp( dha * (tk - 298.15) / (298.15 * kR * tk) )
    # ftemp = exp( dha * (tc - 25.0)/(298.15 * kR * (tc + 273.15)) )
    # ftemp = exp( (dha/kR) * (1/298.15 - 1/tk) )
    
    kR   <- 8.3145     # Universal gas constant, J/mol/K
    ftemp <- exp( dha * (tk - tkref) / (tkref * kR * tk) )
    
    return(ftemp)
}

# .............................................................................
dampen_vec <- function( vec, tau ){
    
    if (length(vec)<365){
        
        rlang::abort("dampen_vec(): Aborting. Length of argument 'vec' must be at least 365.")
        
    } else {
        
        ## Add one year's data to the head of 'vec' to avoid "boundary effects"
        vec <- c(vec[1:365], vec)
        
        ## apply dampening
        vec_damped <- rep(NA, length(vec))
        vec_damped[1] <- vec[1] # this is the unwanted "boundary effect"
        for (idx in 2:length(vec)){
            dvar <- (1.0/tau) * (vec[idx] - vec_damped[idx - 1])
            vec_damped[idx] <- vec_damped[idx - 1] + dvar
        }
        
        ## remove the first year again, that was added artifically before
        vec_damped <- vec_damped[366:length(vec_damped)]
        
    }
    
    return( vec_damped )
}
# .............................................................................

calc_ftemp_kphio <- function( tc, c4 = FALSE ){
    
    if (c4){
        ftemp = -0.008 + 0.00375 * tc - 0.58e-4 * tc^2   # Based on calibrated values by Shirley
    } else {
        ftemp <- 0.352 + 0.022 * tc - 3.4e-4 * tc^2
    }
    
    ## avoid negative values
    ftemp <- ifelse(ftemp < 0.0, 0.0, ftemp)
    
    return(ftemp)
}

# .............................................................................

calc_ftemp_inst_vcmax <- function( tcleaf, tcgrowth, tcref = 25.0, method_ftemp = "kumarathunge19" ){
    
    # local parameters
    Hd     <- 200000 # deactivation energy (J/mol)
    Rgas   <- 8.3145 # universal gas constant (J/mol/K)
    tkref  <- tcref + 273.15  # to Kelvin
    tkleaf <- tcleaf + 273.15  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
    
    # Kattge2007 Parametrization
    Ha    <- 71513  # activation energy (J/mol)
    a_ent <- 668.39 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    b_ent <- 1.07   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    
    if (method_ftemp == "kumarathunge19"){
        # Kumarathunge2019 Implementation:
        # local parameters
        a_ent = 645.13 # offset of entropy vs. temperature relationship (J/mol/K)
        b_ent = 0.38   # slope  of entropy vs. temperature relationship (J/mol/K^2)
        
        # local variables
        Ha = 42600 + (1140 * tcgrowth) # Acclimation for vcmax
    }
    
    # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
    dent <- a_ent - (b_ent * tcgrowth)  # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
    fva <- calc_ftemp_arrh( tkleaf, Ha, tkref = tkref )
    fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
    fv  <- fva * fvb
    
    return( fv )
}

# .............................................................................

calc_ftemp_inst_jmax <- function( tcleaf, tcgrowth, tchome = NA, tcref = 25.0, method_ftemp = "kumarathunge19" ){
    
    # loal parameters
    Hd    <- 200000 # deactivation energy (J/mol)
    tkref <- tcref + 273.15  # to Kelvin
    tkleaf <- tcleaf + 273.15  # conversion of temperature to Kelvin, tcleaf is the instantaneous leaf temperature in degrees C.
    
    # Kattge2007 Parametrization
    Ha    <- 49884  # activation energy (J/mol)
    Rgas  <- 8.3145 # universal gas constant (J/mol/K)
    a_ent <- 659.70 # offset of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K)
    b_ent <- 0.75   # slope of entropy vs. temperature relationship from Kattge & Knorr (2007) (J/mol/K^2)
    
    # calculate entropy following Kattge & Knorr (2007), negative slope and y-axis intersect is when expressed as a function of temperature in degrees Celsius, not Kelvin !!!
    dent <- a_ent - b_ent * tcgrowth   # 'tcgrowth' corresponds to 'tmean' in Nicks, 'tc25' is 'to' in Nick's
    
    if(method_ftemp == "kumarathunge19"){
        # Kumarathunge2019 Implementation:
        # local parameters
        Ha    = 40710  # activation energy (J/mol)
        a_ent = 658.77 # offset of entropy vs. temperature relationship (J/mol/K)
        b_ent = 0.84   # slope of entropy vs. temperature relationship (J/mol/K^2)
        c_ent = 0.52   # 2nd slope of entropy vs. temperature (J/mol/K^2)
        
        # Entropy calculation, equations given in Celsius, not in Kelvin
        dent = a_ent - (b_ent * tchome) - c_ent * (tcgrowth - tchome)
    }
    
    fva <- calc_ftemp_arrh( tkleaf, Ha, tkref = tkref )
    fvb <- (1 + exp( (tkref * dent - Hd)/(Rgas * tkref) ) ) / (1 + exp( (tkleaf * dent - Hd)/(Rgas * tkleaf) ) )
    fv  <- fva * fvb
    
    return( fv )
}

# HELPER FUNCTIONS ####
# .............................................................................
co2_to_ca <- function( co2, patm ){
    # Input:    - float, annual atm. CO2, ppm (co2)
    #           - float, monthly atm. pressure, Pa (patm)
    # Output:   - ca in units of Pa
    # Features: Converts ca (ambient CO2) from ppm to Pa.
    ca   <- ( 1.0e-6 ) * co2 * patm         # Pa, atms. CO2
    return( ca )
}
# .............................................................................
QUADP <- function(A,B,C){
    
    if (any(is.na(c(A,B,C)))){
        return(NA)
    } else {
        if((B^2 - 4*A*C) < 0){
            warning("IMAGINARY ROOTS IN QUADRATIC")
            return(0)
        }
        
        if(identical(A,0)){
            if(identical(B,0)){
                return(0)
            } else {
                return(-C/B)
            }
        } else {
            return((- B + sqrt(B^2 - 4*A*C)) / (2*A))
        }
    }
    
}
# .............................................................................

QUADM <- function(A,B,C){
    
    if (any(is.na(c(A,B,C)))){
        return(NA)
    } else {
        if((B^2 - 4*A*C) < 0){
            warning("IMAGINARY ROOTS IN QUADRATIC")
            return(0)
        }
        
        if(identical(A,0)){
            if(identical(B,0)){
                return(0)
            } else {
                return(-C/B)
            }
        } else {
            return((- B - sqrt(B^2 - 4*A*C)) / (2*A))
        }
    }
    
}
# .............................................................................
Tk <- function(x){x+273.15}

# .............................................................................
# Calculation of vapor pressure in Pa
esat <- function(TdegC, Pa=101){  
    a <- 611.21
    b <- 17.502
    c <- 240.97
    f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
    esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
    return(esatval)
}

# .............................................................................
# OPTIMIZATION FUNCTIONS ####
# .........................................................

calc_optimal_tcleaf_vcmax_jmax <- function(tc_leaf = 25,
                                           patm = 101325,
                                           co2 = 400,
                                           vpd = 1000,
                                           ppfd = 130,
                                           fapar = 1,
                                           kphio = 0.05,
                                           beta = 146,
                                           c_cost = 0.41,
                                           vcmax_start = 20,
                                           gs_start = 0.5,
                                           jmax_start = 20,
                                           method_jmaxlim_inst = "smith37") {
    
    optimise_this_tcleaf_vcmax_jmax <- function( par, args, iabs, kphio, beta, c_cost, method_jmaxlim_inst, maximize=FALSE, return_all=FALSE ){
        
        ## Parameters to be optimized
        vcmax <- par[1]
        gs    <- par[2]
        jmax  <- par[3]
        
        ## Arguments to calculate variables
        tc_leaf <- args[1]
        patm    <- args[2]
        co2     <- args[3]
        vpd     <- args[4]
        
        ## Local variables based on arguments -> TODO: maybe move this to outside?
        kmm       <- calc_kmm(tc_leaf, patm)
        gammastar <- calc_gammastar(tc_leaf, patm)
        ns_star   <- calc_viscosity_h2o(tc_leaf, patm) / calc_viscosity_h2o(25, 101325)
        ca        <- co2_to_ca(co2, patm)
        vpd       <- vpd
        kphio     <- kphio * calc_ftemp_kphio( tc_leaf, c4 = F )
        
        ## Electron transport is limiting
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
        
        ## Rubisco is limiting
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
        
        ## Take minimum of the two assimilation rates and maximum of the two ci
        assim <- min( a_j, a_c )
        ci <- max(ci_c, ci_j)
        
        ## only cost ratio is defined. for this here we need absolute values. Set randomly
        cost_transp <- 1.6 * ns_star * gs * vpd
        cost_vcmax  <- beta * vcmax
        cost_jmax   <- c_cost * jmax
        
        ## Option B: This is equivalent to the P-model with its optimization of ci:ca.
        if (assim<=0) {
            net_assim <- -(999999999.9)
        } else {
            net_assim <- -(cost_transp + cost_vcmax + cost_jmax) / assim
        }
        
        if (maximize) net_assim <- -net_assim
        
        # print(par)
        # print(net_assim)
        
        if (return_all){
            return(tibble(vcmax_mine = vcmax, jmax_mine = jmax, gs_mine = gs, ci_mine = ci, chi_mine = ci / ca, a_c_mine = a_c, a_j_mine = a_j, assim = assim, ci_c_mine = ci_c, ci_j_mine = ci_j, cost_transp = cost_transp, cost_vcmax = cost_vcmax, cost_jmax = cost_jmax, net_assim = net_assim, method_jmaxlim_inst = method_jmaxlim_inst))
        } else {
            return( net_assim )
        }
    }
    
    out_optim <- optimr::optimr(
        par        = c( vcmax_start,       gs_start,       jmax_start ), # starting values
        lower      = c( vcmax_start*0.001, gs_start*0.001, jmax_start*0.001 ),
        upper      = c( vcmax_start*1000,  gs_start*1000,  jmax_start*1000 ),
        fn         = optimise_this_tcleaf_vcmax_jmax,
        args       = c(tc_leaf, patm, co2, vpd),
        iabs       = (ppfd * fapar),
        kphio      = kphio,
        beta       = beta,
        c_cost     = c_cost/4,
        method_jmaxlim_inst = method_jmaxlim_inst,
        method     = "L-BFGS-B",
        maximize   = TRUE,
        control    = list( maxit = 10000 )
    )
    
    varlist <- optimise_this_tcleaf_vcmax_jmax(
        par = out_optim$par,
        args = c(tc_leaf, patm, co2, vpd),
        iabs = (fapar * ppfd),
        kphio,
        beta,
        c_cost / 4,
        method_jmaxlim_inst,
        maximize = FALSE,
        return_all = TRUE
    )
    
    return(out_optim)
}

# .............................................................................
LeafEnergyBalance <- function(Tleaf = 21.5,  # Input in degC
                              Tair = 20,     # Input in degC
                              gs = 0.30,     # Input in mol/m2/d/Pa
                              PPFD = 130,    # Input in mol/m2/d
                              VPD = 1000,    # Input in Pa
                              Patm = 101325, # Input in Pa
                              Wind = 2,      # Input in m/s
                              Wleaf = 0.02,
                              StomatalRatio = 1, # 2 for amphistomatous
                              LeafAbs = 0.5, # in shortwave range, much less than PAR
                              returnwhat = c("balance", "fluxes")
) {
    
    
    ## Edits to make function runnable within rpmodel
    ## Correct input units to fit calculations below
    PPFD <- PPFD * 10^6 / (3600*24)               # mol/m2/d to umol/m2/s
    gs   <- gs   * 10^-6 * Patm / (3600*24)       # mol/m2/d/Pa to ppm mol/m2/s
    Patm <- Patm / 1000                           # Pa to kPa
    VPD  <- VPD  / 1000                           # Pa to kPa
    
    
    ## Original function:
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
        
        out <- list(tc_leaf       = Tleaf,     # NEW
                    tc_leaf_star  = Tleaf2,    # NEW
                    eps           = EnergyBal) # NEW
        return(out)                            # NEW
    }
    
    if(returnwhat == "fluxes"){
        
        l <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)
        return(l)
    }
}

calc_tc_leaf_from_tc_air <- function(tc_air   = 25,   # input in degC
                                     gs       = 0.30, # input in mol/m2/d/Pa
                                     vpd      = 1000, # input in Pa
                                     patm     = 101325, # input in Pa
                                     ppfd     = 130,  # input in mol/m2/d
                                     wind     = 2,    # input in m/s
                                     wleaf    = 0.02,
                                     stoma_r  = 1, 
                                     leaf_abs = 0.5) {
    
    
    # LeafEnergyBalance is equivalent to "maximize_this_tc_leaf" with its Tleaf getting optimized
    
    sol_optimr <-	optimr::optimr(
        # Parameter boundaries to optimize within:
        par       = 15,
        lower     = 15 - tc_air, # OLD: 0 
        upper     = 15 + tc_air, # OLD: 40 
        
        # Function to optimize and its inputs:
        fn        = LeafEnergyBalance,
        Tair      = tc_air,    # input in degC
        gs        = gs,        # input in mol/m2/d/Pa
        VPD       = vpd,       # input in Pa
        Patm      = patm,      # input in Pa
        PPFD      = ppfd,      # input in mol/m2/d
        Wind      = wind,      # input in m/s
        Wleaf     = wleaf,
        StomatalRatio = stoma_r,        # 2 for amphistomatous
        LeafAbs   = leaf_abs,
        
        # Optimr settings:
        method    = "L-BFGS-B",
        control   = list( maxit = 10000, maximize = TRUE )
    )
    
    return(sol_optimr$par)
}


maximize_this_tc_leaf <- function(tc_leaf   = 25, # This gets optimized in optimr()
                                  tc_air    = 25,
                                  tc_growth = 25,
                                  tc_home   = 25,
                                  ppfd      = 130, # mol/m2/d
                                  fapar     = 1,
                                  co2       = 400,
                                  patm      = 101325,
                                  vpd       = 1000,
                                  kphio     = 0.05,
                                  toggles   = NA,
                                  wind      = 2,
                                  wleaf     = 0.02,
                                  stoma_r   = 1,
                                  leaf_abs  = 0.5,
                                  method_jmaxlim_inst = "smith37",
                                  method_eb = "plantecophys",
                                  beta      = 146,
                                  c_cost    = 0.41) {
    
    varlist_opt_tcleaf <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_leaf,
                                                         patm = patm,
                                                         co2 = co2,
                                                         vpd = vpd,
                                                         ppfd = ppfd,
                                                         fapar = fapar,
                                                         kphio = kphio,
                                                         beta = beta,
                                                         c_cost = c_cost,
                                                         method_jmaxlim_inst = method_jmaxlim_inst,
                                                         vcmax_start = 20,
                                                         gs_start = 0.5,
                                                         jmax_start = 20)
    gs <- varlist_opt_tcleaf$gs_mine
    
    if (method_eb == "plantecophys") {
        ## Via plantecophys energy balance
        tc_leaf_x <- calc_tc_leaf_from_tc_air(tc_air = tc_air,
                                              # gs       = 0.15 / 101325,
                                              gs       = gs,
                                              wind     = wind,
                                              wleaf    = wleaf,
                                              stoma_r  = stoma_r,
                                              leaf_abs = leaf_abs)
        
    } else if (method_eb == "tealeaves") {
        ## Via tealeaves energy balance
        
        # Get relative humidity from vpd:
        RH <- (100 - ((vpd * 100) / esat(tc_air, patm/100))) / 100
        
        # Get leaf parameters
        leaf_par <- make_leafpar(
            replace = list(
                g_sw = set_units(gs, "mol/m^2/d/Pa")))
        
        # Get environmental parameters
        enviro_par <- make_enviropar(
            replace = list(
                T_air = set_units(tc_air + 273.15, "K"),
                RH    = as_units(RH),
                P     = set_units(patm, "Pa")))
        
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
                               method_eb = "plantecophys",
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
                                     method_eb = method_eb,
                                     beta      = beta,
                                     c_cost    = c_cost)
        },
        warning = function(cond){
            message("There was a warning")
            return(NA)
        },
        error = function(cond){
            message("This message will not be printed.")
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

# PLOTTING FUNCTIONS ####
# .............................................................................

plot_obs_vs_pred_mina <- function(df_in, unnest_var = "fit_sim") {
    
    df_temp <- df_in %>% unnest(unnest_var)
    
    fit     <- lm(tc_opt_obs ~ tc_opt_sim, data = df_temp)
    sry     <- summary(fit)
    r2      <- sry$adj.r.squared %>% round(2)
    rmse    <- sqrt( ( c(crossprod(sry$residuals)) / length(sry$residuals) ) ) %>% round(2)
    slope   <- sry$coefficients[2, 1] %>% round(2)
    b0_is_0 <- !(sry$coefficients[1, 4] <= 0.05) # T = Intercept is 0
    b1_is_1 <- (confint(fit, 'tc_opt_sim', level=0.95)[1] < 1) & (confint(fit, 'tc_opt_sim', level=0.95)[2] > 1) # T = Slope is 1
    bias    <- mean(df_temp$tc_opt_sim - df_temp$tc_opt_obs) %>% round(2)
    
    n_ac <- df_temp %>% dplyr::select(sitename, date, tc_opt_sim, min_a_sim) %>% distinct() %>% dplyr::filter(min_a_sim == "ac") %>% nrow()
    n_aj <- df_temp %>% dplyr::select(sitename, date, tc_opt_sim, min_a_sim) %>% distinct() %>% dplyr::filter(min_a_sim == "aj") %>% nrow()
    
    df_metrics <- tibble(r2 = r2,
                         rmse = rmse,
                         slope = slope,
                         b0_is_0   = b0_is_0,
                         b1_sign_not_one = b1_is_1,
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
        
        geom_smooth(data = df_temp, aes(y = tc_opt_obs, x = tc_opt_sim), method = "lm", size = 0.5, fullrange = T, color = "black") +
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
             subtitle = bquote(R^2  ~  " = "  ~ .(r2)  ~  " | bias = "  ~ .(bias)  ~  " | slope = "  ~ .(slope)  ~  " | rmse = "  ~ .(rmse)),
             caption = paste0("Intercept = 0: ", ifelse(b0_is_0, "Y", "N"), " | slope = 1: ", ifelse(b1_is_1, "Y", "N"),
                              " | #Lim: Ac = ", n_ac, ", Aj = ", n_aj)
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
# .............................................................................

plot_topt_vs_tgrowth <- function(df_in,
                                 unnest_var = c("fit_sim", "forcing"), 
                                 tc_opt = "tc_opt_sim") { # tc_opt: tc_opt_sim, tc_opt_obs
    
    ## Unnest data
    df_temp <- df_in %>% unnest(unnest_var)
    df_temp$tc_opt <- df_temp[[tc_opt]]
    
    ## Checks
    energy_balance_on <- ifelse(unique(df_temp$tc_growth_air == df_temp$tc_growth_leaf) != T, T, F) # T = energy balance model was on
    n_aj <- df_temp %>% dplyr::select(sitename, date, tc_opt_obs, min_a_sim) %>% distinct() %>% dplyr::filter(min_a_sim == "aj") %>% nrow()
    if (n_aj != 0) {message("Note: Aj limitations in tc_opt estimates but not visualized.")}
    
    
    ## Get intersection points
    ### For air temperatures
    sry     <- lm(tc_opt ~ tc_growth_air, data = df_temp) %>% summary()
    q       <- sry$coefficients[1, 1]
    m       <- sry$coefficients[2, 1]
    x       <- seq(0, 40, 0.1)
    y       <- round(x*m + q, 1)
    crs_air <- x[x == y] 
    if (length(crs_air) > 1) {crs_air <- crs_air[1]}
    
    ### For leaf temperatures
    sry     <- lm(tc_opt ~ tc_growth_leaf, data = df_temp) %>% summary()
    q       <- sry$coefficients[1, 1]
    m       <- sry$coefficients[2, 1]
    x       <- seq(0, 40, 0.1)
    y       <- round(x*m + q, 1)
    crs_leaf<- x[x == y] 
    if (length(crs_leaf) > 1) {crs_leaf <- crs_leaf[1]}
    
    
    ## Create plots
    if (energy_balance_on) {
        p_temp <- df_temp %>%
            pivot_longer(cols = c(tc_growth_leaf, tc_growth_air), names_to = "condition", values_to = "tc_growth") %>% 
            ggplot(aes(x = tc_growth, y = tc_opt, color = condition)) +
            scale_color_manual(values = c(tc_growth_air = "dodgerblue", tc_growth_leaf = "seagreen")) +
            labs(caption = paste("Intersection with one-to-one line: tc_growth_air at ", crs_air, "°C | tc_growth_leaf at ", crs_leaf, "°C"))
        
    } else {
        p_temp <- df_temp %>%
            ggplot(aes(x = tc_growth_air, tc_opt))
    }
    
    p_out <- p_temp +
        geom_abline(linetype = "dotted") +
        geom_point() +
        geom_smooth(method = "lm", fullrange = T)+
        xlim(0, 40) +
        ylim(0, 40) +
        xlab("T_growth [°C]") +
        ylab("T_opt simulated [°C]") +
        ggpmisc::stat_poly_eq(formula = y ~ x,
                              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
                              parse = TRUE) +
        theme_bw() +
        theme(legend.position = "bottom")

    out <- list(plot = p_out, df_metrics   = NA)
    
    return(out)
}

# .............................................................................
make_df_evaluation_p21_long <- function(df_in, ftemp_method = "kumarathunge19") {
    
    df_evaluation_p21_long <- df_in %>%
        mutate(tc_growth_air = purrr::map(forcing, ~ dplyr::select(., tc_growth_air)),
               tc_home = purrr::map(forcing, ~ dplyr::select(., tc_home))) %>% 
        unnest(c(tc_growth_air, tc_home)) %>% 
        mutate(tc_leaf   = purrr::map_dbl(data_raw, ~ mean(.$tleaf)),
               jmax      = purrr::map_dbl(data_raw, ~ mean(.$jmax))  * 10 ^-6,
               vcmax     = purrr::map_dbl(data_raw, ~ mean(.$vcmax)) * 10 ^-6,
               jmax25    = jmax  / calc_ftemp_inst_jmax(tc_leaf, tc_growth_air, tc_home, method_ftemp = ftemp_method),
               vcmax25   = vcmax / calc_ftemp_inst_vcmax(tc_leaf, tc_growth_air, method_ftemp = ftemp_method)) %>% 
        pivot_longer(cols = c(jmax, vcmax, jmax25, vcmax25), names_to = "variable", values_to = "peng21") %>% 
        dplyr::select(sitename, date, tc_growth_air, tc_home, tc_leaf, variable, peng21) 
    
    return(df_evaluation_p21_long)
}


make_long_df <- function(df_in,
                         v_unnest = c("rpm_accl", "forcing"),
                         v_vars   = c("vcmax", "jmax"),
                         dataset  = "dataset") {
    
    df_out <- df_in %>%
        unnest(v_unnest) %>%
        pivot_longer(cols = v_vars, names_to = "variable", values_to = "values") %>% 
        dplyr::select(sitename, date, variable, values)
    
    names(df_out)[names(df_out)=="values"] <- dataset
    
    return(df_out)
}

plot_two_long_df <- function(df_x,
                             df_x_dataset,
                             df_y,
                             df_y_dataset) {
    
    names(df_x)[names(df_x) == df_x_dataset] <- "x"
    names(df_y)[names(df_y) == df_y_dataset] <- "y"
    max        <- max(df_x$x, df_y$y)
    df_temp    <- left_join(df_x, df_y)
    
    p <- df_temp  %>% 
        ggplot() +
        aes(x = x, y = y) +
        geom_abline() +
        geom_point() +
        geom_smooth(method = "lm", fullrange = T) +
        ggpmisc::stat_poly_eq(formula = y ~ x,
                              aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                              parse = TRUE) +
        xlab(paste(df_x_dataset)) +
        ylab(paste(df_y_dataset)) +
        ylim(0, max) +
        xlim(0, max) +
        facet_wrap(~variable, scales = "free") 
    
    return(p)
}

get_instant_vcmax_jmax <- function(df_in, ftemp_method = "kumarathunge19") {
    
    df_out <- df_in %>%
        unnest(c(rpm_accl, forcing)) %>%
        mutate(tc_leaf_dat = purrr::map_dbl(data_raw, ~ mean(.$tleaf)),
               jmax    = jmax25  * calc_ftemp_inst_jmax(tc_leaf_dat, tc_growth_air, tc_home, method_ftemp = ftemp_method),
               vcmax   = vcmax25 * calc_ftemp_inst_vcmax(tc_leaf_dat, tc_growth_air, method_ftemp = ftemp_method)) %>% 
        nest(rpm_accl = c("chi", "ci", "xi", "gs", "vcmax", "vcmax25", "jmax", "jmax25", "kphio")) %>% 
        dplyr::select(sitename, date, rpm_accl) %>% 
        left_join(df_in %>% dplyr::select(-rpm_accl))
    
    return(df_out)
}

# TEST ZONE -------------------------------------------------------------------

if (F) { # F to avoid calculation when sourcing
    
    # Parameters:
    tc_air    <- 20
    ppfd      <- 130
    co2       <- 400
    patm      <- 101325
    vpd       <- 1000
    kphio		  <- 0.05
    fapar     <- 1
    method_jmaxlim_inst <- "smith37"
    method_eb <- "plantecophys"
    beta      <- 146
    c_cost    <- 0.41
    
    vcmax_start <- 20
    jmax_start  <- 20
    gs_start    <- 0.5
    
    tc_leaf <- calc_tc_leaf_final(tc_air = tc_air,
                                  ppfd = ppfd,
                                  fapar = fapar,
                                  co2 = co2,
                                  patm = patm,
                                  vpd = vpd,
                                  kphio = kphio,
                                  method_jmaxlim_inst = method_jmaxlim_inst,
                                  beta = beta,
                                  c_cost = c_cost,
                                  method_eb = method_eb)
    
    varlist_final <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_leaf,
                                                    patm = patm,
                                                    co2 = co2,
                                                    vpd = vpd,
                                                    ppfd = ppfd,
                                                    fapar = fapar,
                                                    kphio = kphio,
                                                    beta = beta,
                                                    c_cost = c_cost,
                                                    vcmax_start = vcmax_start,
                                                    gs_start = gs_start,
                                                    jmax_start = jmax_start,
                                                    method_jmaxlim_inst = method_jmaxlim_inst)
    
    
    tc_array  <- seq(0, 40, 1)
    var_array <- rep(NA, length(tc_array))
    
    for (t in 1:length(tc_array)) {
        var_array[t] <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_array[t],
                                                       patm = patm,
                                                       co2 = co2,
                                                       vpd = vpd,
                                                       ppfd = ppfd,
                                                       fapar = fapar,
                                                       kphio = kphio,
                                                       beta = beta,
                                                       c_cost = c_cost,
                                                       vcmax_start = vcmax_start,
                                                       gs_start = gs_start,
                                                       jmax_start = jmax_start,
                                                       method_jmaxlim_inst = method_jmaxlim_inst)$vcmax_mine
    }  
    
    plot(tc_array, var_array)
    
}
