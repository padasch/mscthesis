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
calc_aj <- function(kphio, ppfd, jmax, gammastar, ci = NA, ca, fapar, theta = 0.85, j_method = "smith37", model = "analytical", gs) {
    
    ## Get iabs ....................................................................................
    iabs <- ppfd * fapar
    
    ## Get J based on Jmax-Limitation formulation ..................................................
    ## smith37:    Jmax Limitation following Smith et al. (1937):
    if (j_method == "smith37") {j <- (4 * kphio * iabs) / sqrt(1 + ((4 * kphio * iabs)/(jmax))^2)}
    
    ## farquhar89: Jmax Limitation following Farquhar et al. (1989):
    if (j_method == "farquhar89") {j <- (kphio * iabs + jmax - sqrt(( kphio * iabs + jmax)^2 - (4*kphio*theta*iabs*jmax))) / (2*theta)}
    
    
    ## Find ci if not given
    if (model == "numerical" | is.na(ci)) {
        ## Solve Eq. system:
        ## A = gs (ca - ci)
        ## A = j/4 * (ci-gammastar)/ci+2*gammastar)
        ## This leads to a quadratic equation:
        ## A * ci^2 + B * ci + C  = 0
        ## 0 = a + b*x + c*x^2
        ## with
        A  <- -gs
        B  <- gs * ca - 2 * gammastar * gs - j/4
        C  <- 2 * gammastar * gs * ca + gammastar * j/4
        ci <- QUADM(A, B, C)
    }
    
    aj   <- j/4 * (ci - gammastar)/(ci + 2 * gammastar)
    
    out = list(aj   = aj,
               j    =  j,
               ci   = ci)
    
    return(out)
    
    # Test Case:
    calc_aj(kphio = 0.1, ppfd = 130, jmax=8, gs = 0.5,
            gammastar=4, ci = NA, ca=40, fapar=1, theta = 0.85, 
            j_method = "smith37", model = "analytical")
}

# .............................................................................
calc_ac <- function(ci = NA, ca, gammastar, kmm, vcmax, model = "analytical", gs) {
    
    ## Find ci if not given: 
    if (model == "numerical" | is.na(ci)) {
        ## Solve Eq. system:
        ## A = gs (ca- ci)
        ## A = Vcmax * (ci - gammastar)/(ci + Kmm)
        ## This leads to a quadratic equation:
        ## A * ci^2 + B * ci + C  = 0
        ## 0 = a + b*x + c*x^2
        ## with
        A  <- -1.0 * gs
        B  <- gs * ca - gs * kmm - vcmax
        C  <- gs * ca * kmm + vcmax * gammastar
        ci <- QUADM(A, B, C)
    }
    
    ac  <- vcmax * (ci - gammastar)/(ci + kmm)
    out <- list(ac = ac, ci = ci)
    
    return(out)
    
    
}

# .............................................................................
calc_ac <- function(ci, ca, gammastar, kmm, vcmax, model = "analytical", gs) {
    
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
        f <- (3.22 - 0.046 * tc_leaf)^(tc_leaf - 25.0)/10 # Q10 changes with tc_leaf acc. to Tjoelker et al. (2001)
    }
    
    rd <- vcmax25 * rd_to_vcmax * f
    
    return(rd)
}
# .............................................................................
calc_jmax <- function(kphio, iabs, ci, gammastar, method, theta = NA, c = NA){ 
    
    if (method == "smith37" | method == "wang17") {
        c      <- 0.103 # Estimated by Wang et al. (2017)
        c_star <- 4*c
        jmax   <- 4*kphio*iabs / (sqrt(1 / (1 - (c_star*(ci + 2*gammastar)/(ci-gammastar))^(2/3)) - 1))  
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
    ## Pa in kPa
    a <- 611.21
    b <- 17.502
    c <- 240.97
    f <- 1.0007 + 3.46 * 10^-8 * Pa * 1000
    esatval <- f * a * (exp(b * TdegC/(c + TdegC)))
    return(esatval)
}

VPDtoRH <- function(VPD, TdegC, Pa=101){
    ## VPD and Pa in kPa
    esatval <- esat(TdegC, Pa)
    e <- pmax(0, esatval - VPD*1000)
    RH <- 100 * e/esatval
    return(RH)
}

# .............................................................................
# ENERGY BALANCE + NUMERICAL####
# .........................................................
calc_optimal_gs_vcmax_jmax <- function(tc_air, # Input in: degC
                                       patm, # Input in: Pa
                                       co2, # Input in: Pa
                                       vpd, # Input in: Pa
                                       ppfd , # Input in: mol/m2/s
                                       fapar, # Input in: -
                                       kphio, # Input in: -
                                       beta = 146, # Input in: -
                                       vcmax_start = 15,  # Input in: -
                                       jmax_start  = 40,  # Input in: -
                                       gs_start    = 0.4, # Input in: -
                                       method_jmaxlim_inst = "smith37", 
                                       method_eb = "off",
                                       which_optimizer = "optimr" # optimr or gensa
                                       ) { 
    
    ## .................................................................................................
    ## Pick Optimizer and set settings based on energy balance model input
    if (which_optimizer == "optimr") {
        if (method_eb == "tealeaves") {
            ctrl <- list(maxit = 10)
        } else {
            ctrl <- list(maxit = 10000)
        }
        
        out_optim <- optimr::optimr(
            fn         = optimise_this_gs_vcmax_jmax,
            par        = c( vcmax_start,       jmax_start      , gs_start ), # starting values
            lower      = c( vcmax_start/1000,  jmax_start/1000,  gs_start/1000 ),
            upper      = c( vcmax_start*1000,  jmax_start*1000 , gs_start*1000 ),
            args       = c(tc_air, patm, co2, vpd),
            iabs       = (ppfd * fapar * 3600 * 24),
            kphio      = kphio,
            method_jmaxlim_inst = method_jmaxlim_inst,
            method_eb  = method_eb,
            method     = "L-BFGS-B",
            maximize   = TRUE,
            control    = ctrl
        )
    }
    
    if (which_optimizer == "gensa") {
        out_optim <- GenSA(
            fn         = optimise_this_gs_vcmax_jmax,
            par        = c( vcmax_start,       jmax_start      , gs_start ), # starting values
            lower      = c( vcmax_start/1000,  jmax_start/50,  gs_start/1000 ),
            upper      = c( vcmax_start*1000,  jmax_start*50 , gs_start*1000 ),
            args       = c(tc_air, patm, co2, vpd),
            iabs       = (ppfd * fapar * 3600 * 24),
            kphio      = kphio,
            method_jmaxlim_inst = method_jmaxlim_inst,
            method_eb  = method_eb,
            maximize   = TRUE,
        )
    }
    
    ## .................................................................................................
    ## Return optimized variables:
    varlist <- optimise_this_gs_vcmax_jmax(
        par = out_optim$par,
        args = c(tc_air, patm, co2, vpd),
        iabs = (ppfd * fapar * 3600 * 24),
        kphio = kphio,
        method_jmaxlim_inst = method_jmaxlim_inst,
        method_eb = method_eb,
        maximize = FALSE,
        return_all = TRUE
    )
    
    out <- list(out_optim = out_optim,
                varlist   = varlist)
    
    return(out)
    
    ## Test Case:
    calc_optimal_gs_vcmax_jmax(tc_air = 10,
                               patm = 101325,
                               co2 = 400,
                               vpd = 1000,
                               ppfd = 800e-6,
                               fapar = 1,
                               kphio = 0.05,
                               method_jmaxlim_inst = "smith37",
                               method_eb = "tealeaves",
                               which_optimizer = "optimr")
}

# .............................................................................
optimise_this_gs_vcmax_jmax <- function(par,
                                        args,
                                        iabs,
                                        kphio,
                                        beta = 146,
                                        c_cost = NA,
                                        method_jmaxlim_inst,
                                        method_eb = "off",
                                        maximize = FALSE,
                                        return_all = FALSE) {
    
    ## .................................................................................................
    ## Set inputs to local variables
    vcmax   <- par[1]
    jmax    <- par[2]
    gs      <- par[3]
    tc_air  <- args[1]
    patm    <- args[2]
    co2     <- args[3]
    vpd     <- args[4]
    
    ## .................................................................................................
    ## Call energy balance function
    if (method_eb == "off") {
        tc_leaf  <- tc_air
    } else {
        ppfd <- iabs # only true if fapar = 1!
        
        tc_leaf <- calc_tleaf_from_tair_with_fixed_gs(tc_air = tc_air, # degC
                                                      gs_c = gs   / 3600/24,       # mol/m2/s
                                                      ppfd = ppfd / 3600/24,     # mol/m2/s
                                                      vpd = vpd,       # Pa
                                                      patm = patm,     # Pa
                                                      method_eb = method_eb # model: plantecophys, tealeaves
        )
    }
    
    ## .................................................................................................
    ## Get variables from functions
    vpd_leaf  <- VPDairToLeaf(vpd/1000, tc_air, tc_leaf, patm/1000) * 1000
    kmm       <- calc_kmm(tc_leaf, patm)
    gammastar <- calc_gammastar(tc_leaf, patm)
    ns_star   <- calc_viscosity_h2o(tc_leaf, patm) / calc_viscosity_h2o(25, 101325)
    ca        <- co2_to_ca(co2, patm)
    kphio     <- kphio * calc_ftemp_kphio( tc_leaf, c4 = F )
    
    ## .................................................................................................
    ## Calcualte assimilation rates
    ## Aj following Smith (1937):
    if (method_jmaxlim_inst == "smith37") {
        aj_out <- calc_aj(kphio, ppfd = iabs, jmax, gammastar, ci = NA, ca, fapar=1, theta = 0.85, j_method = "smith37", model = "numerical", gs)
        a_j <- aj_out$aj
        ci_j <- aj_out$ci
        c_cost <- 0.103 # As estimated by Wang et al. (2017)
    }
    
    ## Aj following Smith (1937):
    if (method_jmaxlim_inst == "farquhar89") {
        aj_out <- calc_aj(kphio, ppfd = iabs, jmax, gammastar, ci = NA, ca, fapar=1, theta = 0.85, j_method = "farquhar89", model = "numerical", gs)
        a_j <- aj_out$aj
        ci_j <- aj_out$ci
        c_cost <- 0.053 # As estimated by Smith et al. (2019)
    }
    
    ## Ac following FvCB-Model: 
    ac_out <- calc_ac(ci=NA, ca, gammastar, kmm, vcmax, model = "numerical", gs)
    a_c  <- ac_out$ac
    ci_c <- ac_out$ci
    
    ## Take minimum of the two assimilation rates and maximum of the two ci
    assim <- min( a_j, a_c )
    ci    <- max(ci_c, ci_j)
    
    ## .................................................................................................
    ## Calculate individual and total costs
    cost_transp <- 1.6 * ns_star * gs * vpd_leaf
    cost_vcmax  <- beta * vcmax
    cost_jmax   <- c_cost * jmax
    
    if (assim<=0) {
        stop("Assimilation is negative!")
    } else {
        net_assim <- -(cost_transp + cost_vcmax + cost_jmax) / assim
    }
    
    if (maximize) net_assim <- -net_assim
    
    ## .................................................................................................
    ## Gather all variables
    if (return_all) {
        return(
            tibble(
                vcmax_mine = vcmax /(3600*24), # Turn output into per-seconds scale
                jmax_mine = jmax /(3600*24), # Turn output into per-seconds scale
                gs_mine = gs /(3600*24), # Turn output into per-seconds scale
                ci_mine = ci,
                chi_mine = ci / ca,
                a_c_mine = a_c /(3600*24), # Turn output into per-seconds scale
                a_j_mine = a_j /(3600*24), # Turn output into per-seconds scale
                assim = assim /(3600*24), # Turn output into per-seconds scale
                ci_c_mine = ci_c,
                ci_j_mine = ci_j,
                cost_transp = cost_transp,
                cost_vcmax = cost_vcmax,
                cost_jmax = cost_jmax,
                net_assim = net_assim,
                method_jmaxlim_inst = method_jmaxlim_inst,
                tc_leaf = tc_leaf,
                kphio = kphio
            )
        )
    } else {
        return( net_assim )
    }
}

# .............................................................................
LeafEnergyBalance <- function(Tleaf,  # Input in degC
                              Tair,     # Input in degC
                              gs,     # Input in mol/m2/s/Pa
                              PPFD,    # Input in mol/m2/s
                              VPD,    # Input in Pa
                              Patm, # Input in Pa
                              Wind = 2,      # Input in m/s
                              Wleaf = 0.02,
                              StomatalRatio = 2, # 1 = hypostomatous, 2 = amphistomatous
                              LeafAbs = 0.5, # in shortwave range, much less than PAR
                              returnwhat = c("diff")) { # sqrd, abs, diff, flux, verbose
    
    
    ## 0. Prerequisites
    ## 0.1 Scale units to fit input and calculations
    Wleaf <- 0.1                        # Assume same leaf size as tealeaves model
    PPFD  <- PPFD * 10^6                # mol/m2/s to umol/m2/s
    gs    <- gs * Patm                  # mol/m2/s/Pa to mol/m2/s
    Patm  <- Patm / 1000                # Pa to kPa
    VPD   <- VPD  / 1000                # Pa to kPa
    
    ## 0.2 Definition of constants
    Boltz <- 5.67 * 10^-8     # w M-2 K-4
    Emissivity <- 0.97        # - Original value in plantecophys: 0.95, 0.97 to fit tealeaves conditions
    LatEvap <- 2.54           # MJ kg-1
    CPAIR <- 1010.0           # J kg-1 K-1
    
    H2OLV0 <- 2.501e6         # J kg-1
    H2OMW <- 18e-3            # J kg-1
    AIRMA <- 29.e-3           # mol mass air (kg/mol)
    AIRDENS <- 1.204          # kg m-3
    UMOLPERJ <- 2.04          # umol/J from Meek et al. (1984) as ingestr scales, Bonan (2016) suggests 4.6
    DHEAT <- 21.5e-6          # molecular diffusivity for heat
    
    
    ## 1. Total Radiative Uptake
    Rsol    <- 2*PPFD/UMOLPERJ               # Short-wave uptake, W m-2
    Rlongup <- Emissivity*Boltz*Tk(Tleaf)^4  # Long-wave uptake, W m-2, positive is heat loss from leaf
    Rnet    <- LeafAbs*Rsol - Rlongup        # Net radiation uptake (- because Rlongup flux direction)
    
    
    ## 2. Net Radiation
    ## 2.1 Simplified terms
    AIRDENS <- Patm*1000/(287.058 * Tk(Tair))     # Density of dry air
    LHV     <- (H2OLV0 - 2.365E3 * Tair) * H2OMW  # Latent heat of water vapour at air temperature (J mol-1)
    SLOPE   <- (esat(Tair + 0.1) - esat(Tair)) / 0.1 # Const s in Penman-Monteith equation  (Pa K-1)
    Gradiation <- 4.*Boltz*Tk(Tair)^3 * Emissivity / (CPAIR * AIRMA) # Radiation conductance, both leaf sides (mol m-2 s-1)
    CMOLAR  <- Patm*1000 / (8.314 * Tk(Tair))      # .Rgas() in package...
    GRASHOF <- 1.6E8 * abs(Tleaf-Tair) * (Wleaf^3) # Grashof number
    ea      <- esat(Tair) - 1000*VPD               # Vapor pressure
    
    ## 2.2 Boundary layer conductance for heat and for water, see Leuning et al. (1995) PC&E 18:1183-1200 Appendix E
    Gbhforced  <- 0.003 * sqrt(Wind/Wleaf) * CMOLAR             # Boundary conductance under forced convection
    Gbhfree    <- 0.5 * DHEAT * (GRASHOF^0.25) / Wleaf * CMOLAR # Boundary conductance under free convection
    Gbh        <- 2*(Gbhfree + Gbhforced)                       # Total conductance (both leaf sides)
    Gradiation <- Gradiation * 2                                # Radiation conductance on both leaf sides
    Gbhr       <- Gbh + Gradiation                              # Total boundary layer conductance to heat
    Gbw        <- StomatalRatio * 1.075 * Gbh                   # Boundary layer conductance for water (mol m-2 s-1)

    ## 2.3 Isothermal net radiation, see Leuning et al. (1995) PC&E 18:1183-1200 Appendix D
    ema     <- 0.642*(ea/Tk(Tair))^(1/7)                 # Simplified term
    gw      <- gs*Gbw/(gs + Gbw)                         # Simplified term
    Rnetiso <- LeafAbs*Rsol - (1 - ema)*Boltz*Tk(Tair)^4 # Isothermal net radiation, Eq. D1

    ## 2.4 Penmon-Monteith equation for latent heat loss, isothermal form, Eq. 10
    GAMMA <- CPAIR*AIRMA*Patm*1000/LHV
    ET <- (1/LHV) * (SLOPE * Rnetiso + 1000*VPD * Gbh * CPAIR * AIRMA) / (SLOPE + GAMMA * Gbhr/gw)
    lambdaET <- LHV * ET # Latent heat loss
    
    ## 3. Energy balance and leaf temperature using Leuning 1995, Eq. 11 (both times)
    ## 3.1 Heat flux from Gradiation, Eq. 11
    Y <- 1/(1 + Gradiation/Gbh)                          # Psychometric constant
    H2 <- Y*(Rnetiso - lambdaET)                         # Sensible heat flux
    Tleaf2 <- Tair + H2/(CPAIR * AIRDENS * (Gbh/CMOLAR)) # Leaf-air temperature calculated from energy balance
    
    ## 3.2 Sensible heat flux from air-to-leaf temperature difference
    H <- -CPAIR * AIRDENS * (Gbh/CMOLAR) * (Tair - Tleaf) # (positive flux is heat loss from leaf)
    
    
    ## 4. Output for minimization optimum or fluxes:
    if (returnwhat == "sqrd") {out <- (Tleaf - Tleaf2)^2}  # For optimization via optimr(), more accurate than using uniroot()
    if (returnwhat == "abs")  {out <- abs(Tleaf - Tleaf2)} # For optimization via optimr(), needs more iteration than using ()^2
    if (returnwhat == "diff") {out <- (Tleaf - Tleaf2)}    # For optimization via uniroot()
    if (returnwhat == "verbose") {out <- list(tc_leaf = Tleaf, tc_leaf_star  = Tleaf2, eps = EnergyBal)}  # To investigate optimiziation
    if (returnwhat == "fluxes")  {out <- data.frame(ELEAFeb=1000*ET, Gradiation=Gradiation, Rsol=Rsol, Rnetiso=Rnetiso, Rlongup=Rlongup, H=H, lambdaET=lambdaET, gw=gw, Gbh=Gbh, H2=H2, Tleaf2=Tleaf2)}
    
    return(out)
}

## .................................................................................................
calc_tc_leaf_from_tc_air <- function(tc_air   = 25,   # input in degC
                                     gs       = 0.30, # input in mol/m2/d/Pa
                                     vpd      = 1000, # input in Pa
                                     patm     = 101325, # input in Pa
                                     ppfd     = 130   # input in mol/m2/d
                                     ) {
    
    
    # Use uniroots():
    sol_optimr <- try(uniroot(LeafEnergyBalance,
            interval = c(max(7, tc_air-30), tc_air+30),
            Tair      = tc_air,
            gs        = gs,        # input in mol/m2/d/Pa
            VPD       = vpd,       # input in Pa
            Patm      = patm,      # input in Pa
            PPFD      = ppfd,      # input in mol/m2/d
            returnwhat = "diff"))

    out <- sol_optimr$root
    
    return(out)
}


maximize_this_tc_leaf <- function(tc_leaf   = 25, # This gets optimized in optimr()
                                  tc_air    = 25,
                                  ppfd      = 130, # mol/m2/d
                                  fapar     = 1,
                                  co2       = 400,
                                  patm      = 101325,
                                  vpd       = 1000,
                                  kphio     = 0.05,
                                  method_jmaxlim_inst = "smith37",
                                  method_eb = "plantecophys",
                                  gs_water  = 2e-6, # mol/m^2/Pa/s
                                  gs_input  = "rpmodel", # rpmodel or prescribed
                                  method_opt = "diff") {  # sqrd, abs, diff
    
    
    ## Output: difference in tc_leaf assumed for gs and tc_leaf from energy balance
    
    ## 1: Get optimal gs, vcmax and jmax at given tc_leaf
    if (gs_input == "prescribed") {
        gs_water <- gs_water
        
    } else {
        varlist_opt_tcleaf <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = tc_leaf,
                                                             patm = patm,
                                                             co2 = co2,
                                                             vpd = vpd,
                                                             ppfd = ppfd,
                                                             fapar = fapar,
                                                             kphio = kphio,
                                                             method_jmaxlim_inst = method_jmaxlim_inst)$varlist
        
        ## Get optimal conductance to water from conductance to co2
        gs_water <- varlist_opt_tcleaf$gs_mine * 1.6        
    }
    

    ## 2. Call energy balance
    ## 2.1: Use energy balance from plantecophys package
    if (method_eb == "plantecophys") {
        ## 2.1: Via plantecophys energy balance
        tc_leaf_leb <- calc_tc_leaf_from_tc_air(tc_air = tc_air,
                                                gs     = gs_water,
                                                vpd      = vpd,
                                                patm     = patm,
                                                ppfd     = ppfd)
        
    } else if (method_eb == "tealeaves") {
        ## 2.2: Via tealeaves energy balance
        
        # Get relative humidity from vpd
        RH <- VPDtoRH(vpd/1000, tc_air, patm/1000) / 100
        
        # Get incident short-wave radiation flux density from ppfd
        S_sw <- ppfd / (24*3600) * 10^6 / 2.04 # mol/m2/d to umol/m2/s to J/s/m2 = W/m2, 4.6 from Bonan 2016 but Stocker 2020 used 2.04
        
        # Get leaf parameters:
        leaf_par <- make_leafpar(
            replace = list(
                g_sw = set_units(gs_water*10^6, "umol/m^2/s/Pa")))
        
        # Get environmental parameters:
        enviro_par <- make_enviropar(
            replace = list(
                T_air = set_units(tc_air + 273.15, "K"),
                S_sw  = set_units(S_sw, "W/m^2"), 
                RH    = set_units(RH),
                P     = set_units(patm, "Pa")))
        
        # Get physical constants:
        constants  <- make_constants()
        
        # Get tc_leaf:
        tc_leaf_leb <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)$T_leaf %>% 
            set_units("degree_Celsius") %>% drop_units()
    }
    
    # Get difference between tc_leaf and tc_leaf_x
    if (method_opt == "sqrd") {out <- (tc_leaf - tc_leaf_leb)^2}
    if (method_opt == "abs")  {out <- abs((tc_leaf - tc_leaf_leb))}
    if (method_opt == "diff") {out <- (tc_leaf - tc_leaf_leb)}
    
    return(out)
}


calc_tc_leaf_final <- function(tc_air    = 25,
                               ppfd      = 130,
                               fapar     = 1,
                               patm      = 101325,
                               co2       = 400,
                               vpd       = 1000,
                               kphio     = 0.05,
                               method_jmaxlim_inst = "smith37",
                               method_eb = "plantecophys",
                               gs_water  = 2e-6, # mol/m^2/Pa/s
                               gs_input  = "rpmodel") { # rpmodel or prescribed
    
    # Get optimized tc_leaf
    
    # Call optimize()
    sol_optimize <- tryCatch(
        {
            sol_optimize <- uniroot(maximize_this_tc_leaf,
                                     interval  = c(max(1, tc_air-30), tc_air+30),
                                     tc_air    = tc_air,
                                     ppfd      = ppfd,
                                     fapar     = fapar,
                                     co2       = co2,
                                     patm      = patm,
                                     vpd       = vpd,
                                     kphio     = kphio,
                                     method_jmaxlim_inst = method_jmaxlim_inst,
                                     method_eb = method_eb,
                                     gs_water  = gs_water,
                                     gs_input  = gs_input)
        },
        warning = function(cond){
            message("calc_tc_leaf_final(): Warning! Did not converge.")
            return(NA)
        },
        error = function(cond){
            message("calc_tc_leaf_final(): Error!")
            return(NA)
        },
        finally = {
            #pass
        })
    
    if (length(sol_optimize) == 1) {
        return (tc_leaf <- tc_air)
    } else {
        tc_leaf <- sol_optimize$root
    }
    
    return(tc_leaf)
}

## .................................................................................................
calc_tleaf_from_tair_with_fixed_gs <- function(tc_air, # degC
                                               gs_c,   # mol/m2/s
                                               vpd,    # Pa
                                               patm,   # Pa
                                               ppfd,   # mol/m2/s
                                               method_eb # model: plantecophys, tealeaves
                                               ) {
    
    ## .............................................................................................
    ## This function calls either the plantecophy or tealeaves energy balance and returns
    ## the leaf temperature that closes the energy balance.
    ## Function takes a fixed input of gs.
    ## .............................................................................................
    ## Call Check
    if (F) message("Energy balance is called: ", method_eb)
    ## .............................................................................................
    ## Function:
    if (method_eb == "plantecophys") {
        sol_optimize <- tryCatch(
            {
                sol_optimr <- uniroot(LeafEnergyBalance,
                                          interval = c(max(1, tc_air-30), tc_air+30),
                                          Tair      = tc_air,    
                                          gs        = gs_c*1.6,        
                                          VPD       = vpd,       
                                          Patm      = patm,      
                                          PPFD      = ppfd,      
                                          returnwhat = "diff")
            },
            warning = function(cond){
                message("plantecophys: Warning! Did not converge.")
                return(NA)
            },
            error = function(cond){
                message("plantecophys: Error! Did not converge")
                return(NA)
            },
            finally = {
                #pass
            })
        
        if (length(sol_optimize) == 1) {
            return (tc_leaf_out <- tc_air)
        } else {
            tc_leaf_out <- sol_optimize$root
        }
    }
    
    if (method_eb == "tealeaves") {
        # Get relative humidity from vpd
        RH <- VPDtoRH(vpd/1000, tc_air, patm/1000) / 100
        
        # Get incident short-wave radiation flux density from ppfd
        S_sw <- ppfd * 10^6 / 2.04 # umol/m2/s * umol/J = W/m2, 2.04 from Meek et al. (1984)
        
        # Get leaf parameters:
        leaf_par <- make_leafpar(replace = list(
            g_sw = set_units(gs_c*1.6*10^6, "umol/m^2/s/Pa")))
        
        # Get environmental parameters:
        enviro_par <- make_enviropar(replace = list(
            T_air = set_units(tc_air + 273.15, "K"),
            S_sw  = set_units(S_sw, "W/m^2"), 
            RH    = set_units(RH),
            P     = set_units(patm, "Pa")))
        
        # Get physical constants:
        constants  <- make_constants()
        
        # Get tc_leaf:
        tc_leaf_out <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)$T_leaf %>% 
            set_units("degree_Celsius") %>% drop_units()
    }
    
    if (!(method_eb %in% c("plantecophys", "tealeaves"))) {
        stop("> No energy balance model selected")
    }
    
    return(tc_leaf_out)
    
    ## Testcase:
    calc_tleaf_from_tair_with_fixed_gs(tc_air = 4, gs = 5e-6, ppfd = 800e-6, vpd = 1000,
                                       patm = 101325, method_eb = "plantecophys")
}

## .................................................................................................
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
    
    fit     <- lm(y ~ x, data = df_temp)
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
    crs_air <- round(-q/(m-1), 2) # From: x = mx*q
    
    ### For leaf temperatures
    sry     <- lm(tc_opt ~ tc_growth_leaf, data = df_temp) %>% summary()
    q       <- sry$coefficients[1, 1]
    m       <- sry$coefficients[2, 1]
    crs_leaf<- round(-q/(m-1), 2) # From: x = mx*q
    
    
    ## Create plots
    if (energy_balance_on) {
        p_temp <- df_temp %>%
            pivot_longer(cols = c(tc_growth_leaf, tc_growth_air), names_to = "condition", values_to = "tc_growth") %>% 
            ggplot(aes(x = tc_growth, y = tc_opt, color = condition, fill = condition, shape = condition)) +
            scale_color_manual(values = c(tc_growth_air = "dodgerblue", tc_growth_leaf = "seagreen")) +
            scale_fill_manual(values = c(tc_growth_air = "dodgerblue", tc_growth_leaf = "seagreen")) +
            labs(caption = paste("Intersection with one-to-one line: tc_growth_air at ", crs_air, "°C | tc_growth_leaf at ", crs_leaf, "°C"))
        
    } else {
        p_temp <- df_temp %>%
            ggplot(aes(x = tc_growth_air, tc_opt, fill = condition))
    }
    
    p_out <- p_temp +
        geom_abline(linetype = "dotted") +
        geom_point(alpha = 0.9, size = 1.75) +
        geom_smooth(method = "lm", fullrange = T)+
        xlim(0, 40) +
        ylim(0, 40) +
        xlab("T_growth [°C]") +
        ylab(ifelse(tc_opt == "tc_opt_sim", paste0("T_opt simulated [°C]"), paste0("T_opt observed [°C]"))) +
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
    
    standard_error <- function(x) sd(x, na.rm = T) / sqrt(length(x)) # Create own se function
    
    df_evaluation_p21_long <- df_in %>%
        mutate(tc_growth_air = purrr::map(forcing, ~ dplyr::select(., tc_growth_air)),
               tc_home = purrr::map(forcing, ~ dplyr::select(., tc_home))) %>% 
        unnest(c(tc_growth_air, tc_home)) %>% 
        mutate(tc_leaf   = purrr::map_dbl(data_raw, ~ mean(.$tleaf, na.rm = T)),
               jmax      = purrr::map_dbl(data_raw, ~ mean(.$jmax, na.rm = T))  * 10 ^-6,
               vcmax     = purrr::map_dbl(data_raw, ~ mean(.$vcmax, na.rm = T)) * 10 ^-6,
               jmax25    = jmax  / calc_ftemp_inst_jmax(tc_leaf, tc_growth_air, tc_home, method_ftemp = ftemp_method),
               vcmax25   = vcmax / calc_ftemp_inst_vcmax(tc_leaf, tc_growth_air, method_ftemp = ftemp_method),
               jmax_se   = purrr::map_dbl(data_raw, ~ mean(.$vcmax, na.rm = T)) * 10 ^-6,
               vcmax_se  = purrr::map_dbl(data_raw, ~ standard_error(.$vcmax)) * 10 ^-6,
               climate_zone = purrr::map_chr(data_raw, ~magrittr::extract(.$cl) %>% unique),
               climate_zone = as.factor(climate_zone)) %>% 
        pivot_longer(cols = c(jmax, vcmax, jmax25, vcmax25), names_to = "variable", values_to = "peng21") %>% 
        dplyr::select(sitename, date, tc_growth_air, tc_home, tc_leaf, variable, climate_zone, peng21, jmax_se, vcmax_se) 
    
    return(df_evaluation_p21_long)
}


make_long_df <- function(df_in,
                         v_unnest = c("rpm_accl", "forcing"),
                         v_vars   = c("vcmax", "jmax"),
                         dataset  = "dataset") {
    
    df_out <- df_in %>%
        unnest(v_unnest[1]) %>%
        unnest(v_unnest[2]) %>% # Done step-by-step to avoid "Error: Incompatible lengths: 4, 2."
        pivot_longer(cols = v_vars, names_to = "variable", values_to = "values") %>% 
        dplyr::select(sitename, date, variable, values, tc_growth_leaf, tc_growth_air, ppfd_growth, vpd_growth, patm_growth)
    
    names(df_out)[names(df_out)=="values"] <- dataset
    
    return(df_out)
}

plot_two_long_df <- function(df_x,
                             df_x_dataset,
                             df_y,
                             df_y_dataset,
                             kphio_vcmax_corr = T,
                             return_separate = F,
                             model = "add model!") {
    
    
        names(df_x)[names(df_x) == df_x_dataset] <- "x"
        names(df_y)[names(df_y) == df_y_dataset] <- "y"
        df_temp    <- left_join(df_x, df_y)
        
        # Scaling to umol:
        df_temp$x <- df_temp$x * 10^6
        df_temp$y <- df_temp$y * 10^6
        df_temp$vcmax_se <- df_temp$vcmax_se * 10^6
        df_temp$jmax_se <- df_temp$jmax_se * 10^6
        max <- max(df_x$x, df_y$y, na.rm = T) * 10^6
        
        if (kphio_vcmax_corr) {
            ## Linear model for vcmax values
            fit     <- lm(y ~ x, data = dplyr::filter(df_temp, variable == "vcmax"))
            sry     <- summary(fit)
            b0      <- sry$coefficients[1, 1]
            b1      <- sry$coefficients[2, 1]
            df_temp$x <- df_temp$x * b1
            
            message(">> Phi corrected based on fitted slope. kphio scaling factor: ", round(b1, 3))
        }
        
        
    if (return_separate) {
        ## VCMAX ...................................................................................
        ## Subset and Model
        df_subset <- df_temp %>% dplyr::filter(variable == "vcmax")
        df_subset$se <- df_subset$vcmax_se
        
        fit     <- lm(y ~ x, data = df_subset)
        sry     <- summary(fit)
        r2      <- sry$adj.r.squared %>% round(2)
        rmse    <- sqrt( ( c(crossprod(sry$residuals)) / length(sry$residuals) ) ) %>% round(2)
        
        ## Plot
        pvcmax <- df_subset %>%
            ggplot() +
            aes(x = x, y = y) +
            geom_errorbar(aes(ymin = y - se, ymax = y + se, color = climate_zone),
                          alpha = 0.5) +
            geom_point(aes(color = climate_zone),
                       shape = 16,
                       alpha = 0.8) +
            geom_smooth(method = "lm", fullrange = T, color = "black") +
            geom_abline(linetype = "dotted") +
            ylab(bquote("Observed" ~V[vcmax] ~ "["~µmol ~CO[2] ~ m^-2~s^-1~"]"))  +
            xlab(bquote("Predicted" ~V[vcmax] ~ "["~µmol ~CO[2] ~ m^-2~s^-1~"]"))  +
            ylim(0, 200) +
            xlim(0, 200) +
            scale_color_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                                         "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                                         "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                                         "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                                         "#64FFFF", "#F5FFFF"),
                              name = "Koeppen-Geiger Climate Zone",
                              guide = guide_legend(
                                  direction = "horizontal",
                                  title.position = "top",
                                  ncol = 10)) +
            labs(title = bquote(.(model) ~ "- Limitation"),
                 subtitle = bquote(R^2  ~  " = "  ~ .(r2)  ~  " | RMSE = "  ~ .(rmse) ~  " | n = "  ~ .(nrow(df_subset))),) +
            theme(legend.position = "bottom", legend.justification = "center")
            
        ## Jmax ...................................................................................
        ## Subset and Model
        df_subset <- df_temp %>% dplyr::filter(variable == "jmax")
        df_subset$se <- df_subset$jmax_se
        
        fit     <- lm(y ~ x, data = df_subset)
        sry     <- summary(fit)
        r2      <- sry$adj.r.squared %>% round(2)
        rmse    <- sqrt( ( c(crossprod(sry$residuals)) / length(sry$residuals) ) ) %>% round(2)
        
        ## Plot
        pjmax <- df_subset %>%
            ggplot() +
            aes(x = x, y = y) +
            geom_errorbar(aes(ymin = y - se, ymax = y + se, color = climate_zone),
                          alpha = 0.5) +
            geom_point(aes(color = climate_zone),
                       shape = 16,
                       alpha = 0.8) +
            geom_smooth(method = "lm", fullrange = T, color = "black") +
            geom_abline(linetype = "dotted") +
            ylab(bquote("Observed" ~J[max] ~ "["~µmol ~CO[2] ~ m^-2~s^-1~"]"))  +
            xlab(bquote("Predicted" ~J[max] ~ "["~µmol ~CO[2] ~ m^-2~s^-1~"]"))  +
            ylim(0, 200) +
            xlim(0, 200) +
            scale_color_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                                          "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                                          "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                                          "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                                          "#64FFFF", "#F5FFFF"),
                               name = "Koeppen-Geiger Climate Zone",
                               guide = guide_legend(
                                   direction = "horizontal",
                                   title.position = "top",
                                   ncol = 10)) +
            labs(title = bquote(.(model) ~ "- Limitation"),
                 subtitle = bquote(R^2  ~  " = "  ~ .(r2)  ~  " | RMSE = "  ~ .(rmse) ~  " | n = "  ~ .(nrow(df_subset))),) +
            theme(legend.position = "bottom", legend.justification = "center")
        
        
        out <- list(pvcmax = pvcmax,
                    pjmax  = pjmax)
        
        return(out)
        
    } else {
        ## FACET ...................................................................................
        p <- df_temp  %>%
            ggplot() +
            aes(x = x, y = y) +
            geom_point(aes(color = climate_zone)) +
            geom_smooth(method = "lm", fullrange = T) +
            geom_abline() +
            ggpmisc::stat_poly_eq(data = df_temp,
                                  formula = y ~ x,
                                  method = "lm",
                                  aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                                  parse = TRUE) +
            xlab(paste(df_x_dataset)) +
            ylab(paste(df_y_dataset)) +
            ylim(0, max) +
            xlim(0, max) +
            facet_wrap(~variable, scales = "free") +
            scale_color_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                                          "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                                          "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                                          "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                                          "#64FFFF", "#F5FFFF"),
                               name = "Koeppen-Geiger Climate Zone",
                               guide = guide_legend(
                                   direction = "horizontal",
                                   title.position = "top",
                                   ncol = 10))
        
        if (kphio_vcmax_corr) {
            p <- p + labs(caption = (paste0("phi correction factor: ", round(b1, 3))))
        }
        
        return(p)
        
    }
    
}

get_instant_vcmax_jmax <- function(df_in, ftemp_method) {

    names_v <- df_in$rpm_accl[[1]] %>% names()
    
    df_out <- df_in %>%
        unnest(c(rpm_accl, forcing)) %>%
        mutate(tc_leaf_dat = purrr::map_dbl(data_raw, ~ mean(.$tleaf, na.rm = T)),
               jmax    = jmax25  * calc_ftemp_inst_jmax(tc_leaf_dat, tc_growth_air, tc_home, method_ftemp = ftemp_method),
               vcmax   = vcmax25 * calc_ftemp_inst_vcmax(tc_leaf_dat, tc_growth_air, method_ftemp = ftemp_method),
               ) %>% 
        nest(rpm_accl = any_of(names_v)) %>% 
        dplyr::select(sitename, date, rpm_accl) %>% 
        left_join(df_in %>% dplyr::select(-rpm_accl))
    
    return(df_out)
}

# SENSITIVITY ANALYSIS -----------------------------------------------------------------------------
## > RPMODEL ####
rpmodel_env_response <- function(df_ref   = "no_input",
                                 settings = "no_input",
                                 targets  = "no_input",
                                 drivers  = "no_input",
                                 remove_non_convergence = F,
                                 scales_to_one          = T,
                                 p_output               = "per_target",
                                 df_full  = "no_input"){
    
    ## Input Check ####
    ## rpmodel setup
    if (settings[1] == "no_input") {settings <- get_settings()}
    if (targets[1] == "no_input") {targets   <- c("chi", "vcmax", "jmax", "gs", "a_gross")} # , "vcmax25", "jmax25"
    if (drivers[1] == "no_input") {drivers   <- c("tc", "vpd", "ppfd", "patm", "kphio")} # all: c("tc", "vpd", "co2", "ppfd", "patm", "kphio")
    
    ## Environmental setup
    if (df_ref[1] == "no_input") {
        df_ref <- get_df_ref(settings = settings)
    }
    
    
    
    ## Tibble building ####
    ## Environmental setup
    if (!(is_tibble(df_full))) {
        df_ana <- bind_rows(list =
                                tibble(
                                    tc = seq(from = 0, to = 40, length.out = df_ref$nsteps),
                                    tc_home   = rep(df_ref$tc_home, df_ref$nsteps),
                                    vpd       = rep(df_ref$vpd, df_ref$nsteps),
                                    co2       = rep(df_ref$co2, df_ref$nsteps),
                                    ppfd      = rep(df_ref$ppfd, df_ref$nsteps),
                                    patm      = rep(df_ref$patm, df_ref$nsteps),
                                    kphio     = rep(df_ref$kphio, df_ref$nsteps),
                                    env_var   = "tc"),
                            tibble(
                                tc        = rep(df_ref$tc, df_ref$nsteps),
                                tc_home   = seq(from = 0, to = 40, length.out = df_ref$nsteps),
                                vpd       = rep(df_ref$vpd, df_ref$nsteps),
                                co2       = rep(df_ref$co2, df_ref$nsteps),
                                ppfd      = rep(df_ref$ppfd, df_ref$nsteps),
                                patm      = rep(df_ref$patm, df_ref$nsteps),
                                kphio     = rep(df_ref$kphio, df_ref$nsteps),
                                env_var   = "tc_home"),
                            tibble(
                                tc        = rep(df_ref$tc, df_ref$nsteps),
                                tc_home   = rep(df_ref$tc_home, df_ref$nsteps),
                                vpd       = seq(from = 100, to = 3000, length.out = df_ref$nsteps),
                                co2       = rep(df_ref$co2, df_ref$nsteps),
                                ppfd      = rep(df_ref$ppfd, df_ref$nsteps),
                                patm      = rep(df_ref$patm, df_ref$nsteps),
                                kphio     = rep(df_ref$kphio, df_ref$nsteps),
                                env_var   = "vpd"),
                            tibble(
                                tc        = rep(df_ref$tc, df_ref$nsteps),
                                tc_home   = rep(df_ref$tc_home, df_ref$nsteps),
                                vpd       = rep(df_ref$vpd, df_ref$nsteps),
                                co2       = seq(from = 200, to = 1000, length.out = df_ref$nsteps),
                                ppfd      = rep(df_ref$ppfd, df_ref$nsteps),
                                patm      = rep(df_ref$patm, df_ref$nsteps),
                                kphio     = rep(df_ref$kphio, df_ref$nsteps),
                                env_var   = "co2"),
                            tibble(
                                tc        = rep(df_ref$tc, df_ref$nsteps),
                                tc_home   = rep(df_ref$tc_home, df_ref$nsteps),
                                vpd       = rep(df_ref$vpd, df_ref$nsteps),
                                co2       = rep(df_ref$co2, df_ref$nsteps),
                                ppfd      = seq(from = 100e-6, to = 1000e-6, length.out = df_ref$nsteps),
                                patm      = rep(df_ref$patm, df_ref$nsteps),
                                kphio     = rep(df_ref$kphio, df_ref$nsteps),
                                env_var   = "ppfd"),
                            tibble(
                                tc        = rep(df_ref$tc, df_ref$nsteps),
                                tc_home   = rep(df_ref$tc_home, df_ref$nsteps),
                                vpd       = rep(df_ref$vpd, df_ref$nsteps),
                                co2       = rep(df_ref$co2, df_ref$nsteps),
                                ppfd      = rep(df_ref$ppfd, df_ref$nsteps),
                                patm      = seq(from = calc_patm(4000), to = calc_patm(0), length.out = df_ref$nsteps),
                                kphio     = rep(df_ref$kphio, df_ref$nsteps),
                                env_var   = "patm"),
                            tibble(
                                tc        = rep(df_ref$tc, df_ref$nsteps),
                                tc_home   = rep(df_ref$tc_home, df_ref$nsteps),
                                vpd       = rep(df_ref$vpd, df_ref$nsteps),
                                co2       = rep(df_ref$co2, df_ref$nsteps),
                                ppfd      = rep(df_ref$ppfd, df_ref$nsteps),
                                patm      = rep(df_ref$patm, df_ref$nsteps),
                                kphio     = seq(from = 0.01, to = 0.1, length.out = df_ref$nsteps),
                                env_var   = "kphio"))
        
        
        ## Doubling tibble to compare analytical vs. numerical setup
        df_num       <- df_ana
        df_ana$model <- "analytical"
        df_num$model <- "numerical"
        df_full      <- bind_rows(list(df_ana, df_num))
        df_full$formulation <- settings$rpmodel_accl$method_jmaxlim
        
        
        ## Running rpmodel ####
        for (row in 1:nrow(df_full)) {
            message('\014')
            message("Testing environmental response... ", round(row / nrow(df_full) * 100), " %")
            
            if (df_full$model[row] == "analytical") {settings$rpmodel_accl$method_optim <- "analytical"}
            if (df_full$model[row] == "numerical") {settings$rpmodel_accl$method_optim  <- "numerical"}
            
            df_loop <- rpmodel(
                
                ## Inputs
                tc_growth_air = df_full$tc[row],
                tc_home       = df_full$tc_home[row],
                vpd           = df_full$vpd[row],
                co2           = df_full$co2[row],
                ppfd          = df_full$ppfd[row],
                patm          = df_full$patm[row],
                kphio         = df_full$kphio[row],
                
                ## Local parameters
                apar_soilm    = settings$rpmodel_accl$apar_soilm_calib,
                bpar_soilm    = settings$rpmodel_accl$bpar_soilm_calib,
                
                ## Calculation methods
                method_optim   = settings$rpmodel_accl$method_optim, 
                method_jmaxlim = settings$rpmodel_accl$method_jmaxlim,
                method_ftemp   = settings$rpmodel_accl$method_ftemp,
                method_eb      = settings$rpmodel_accl$method_eb
            )
            
            ## Definition of rpmodel output
            df_full$chi[row]     <- df_loop$chi
            df_full$ci[row]      <- df_loop$ci
            df_full$xi[row]      <- df_loop$xi
            df_full$gs[row]      <- df_loop$gs*10^6
            df_full$vcmax[row]   <- df_loop$vcmax*10^6
            df_full$vcmax25[row] <- df_loop$vcmax25*10^6
            df_full$jmax[row]    <- df_loop$jmax*10^6
            df_full$jmax25[row]  <- df_loop$jmax25*10^6
            # df_full$kphio[row]   <- df_loop$kphio
            df_full$a_gross[row] <- df_loop$a_gross*10^6
            df_full$tc_growth_leaf[row]   <- df_loop$tc_growth_leaf
            df_full$opt_convergence[row]   <- df_loop$opt_convergence
        }
    }
    
    ## Wrangling tibble for plotting
    df_temp <- df_full %>%
        mutate(env_var = as.factor(env_var),
               model = as.factor(model),
               formulation = as.factor(formulation),
               model = ifelse(model == "analytical", paste("Analytical"), paste("Numerical")),
               formulation = ifelse(formulation == "farquhar89", paste("Farquhar"), paste("Smith")),
               ppfd = ppfd * 10^6,
               patm = patm / 10^3,
               vpd  = vpd  / 10^3)
    
    ## Rescaling df_ref for display
    df_ref$ppfd <- df_ref$ppfd*10^6
    df_ref$patm <- df_ref$patm/10^3
    df_ref$vpd  <- df_ref$vpd/10^3
    
    if (remove_non_convergence) {
        df_temp <- df_temp %>% dplyr::filter(opt_convergence == 0 | is.na(opt_convergence))
    }
    
    ## For labeling facet wrap:
    vnames <-list("tc" = "T [°C]",
                  "vpd" = "D [kPa]",
                  "co2"  = bquote("CO"[2] ~ " [ppm]"), 
                  "ppfd" = bquote(I[abs] ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]"), 
                  "patm" = bquote(P[atm] ~ "[kPa]"),
                  "kphio"= bquote(Phi ~ "[-]"))
    
    vlabeller <- function(variable,value){
        return(vnames[value])
    }
    
    
    ## Farquhar corrections for kphio:
    if (settings$rpmodel_accl$method_jmaxlim == "farquhar89") {
        df_ref$kphio <- df_ref$kphio * 4
        df_temp$kphio <- df_temp$kphio * 4
        vnames$kphio <- bquote(Phi[j] ~ "[-]")
    }
    
    ## Generating plots ####
    ## Create empty list
    p_all <- list()
    
    
    ## Loop and append plots
    if (p_output == "per_driver") {
        ## Per driver ...........................................................................
        for (v in drivers) {
            
            
            ## Get reference value of driver
            vline <- df_ref[[v]]
            
            
            ## Rescale all target values to maximum across all env. conditions
            if (scales_to_one) {
                df_max <- tibble(max_val = NA, targets = NA)
                
                for (t in targets) {
                    ## Take max value of currently looped env. variable
                    max_val <- df_temp %>% dplyr::filter(env_var == v) %>% dplyr::select(!!t) %>% max()
                    
                    ## Take max value across all env. variables
                    # max_val <- max(df_temp[t])
                    
                    ## Divide by max value
                    df_temp[, t] <- df_temp[t] / max_val
                    
                    ## Add x coodrinate for plotting
                    x  <- (min(df_temp[v]) + max(df_temp[v]))  / 2
                    x <- min(df_temp[v])
                    
                    ## Add all to df
                    df_max <- bind_rows(list(df_max, tibble(max_val = max_val, targets = t, x = x)))
                }
                df_max <- df_max %>%
                    drop_na() %>%
                    mutate(targets = as.factor(targets),
                           max_val = ifelse(max_val < 0.001, round(max_val*10^6, 2), round(max_val, 2)),
                           label   = paste0("Max.: ", max_val),
                           model = "numerical")
            }
            
            
            ## Plotting
            p_append <- df_temp %>%
                dplyr::filter(env_var == v) %>%
                pivot_longer(cols = all_of(v), names_to = "driver", values_to = "driver_value") %>%
                pivot_longer(cols = c(targets), names_to = "targets", values_to = "targets_value") %>% 
                ggplot() +
                aes(x = driver_value, y = targets_value, color = model, shape = formulation) +
                geom_vline(xintercept = vline, alpha = 0.5, color = "grey50", linetype = "solid", size = 0.5) +
                geom_point(size = 2, alpha = 0.8) +
                ylab(paste(" ")) +
                xlab(paste(v)) +
                facet_wrap(~targets, ncol = 3) +
                labs(title = paste0("Sensitivity to ", v, ", using Jmax formulation by ", ifelse(settings$rpmodel_accl$method_jmaxlim == "Smith", "Smith", "Farquhar")),
                     caption = paste0("\n Reference conditions (vertical grey lines): \n",
                                      "tc = ", df_ref$tc, " [°C]",
                                      ", tc_home = ", df_ref$tc_home, " [°C]",
                                      ", vpd = ", df_ref$vpd, " [Pa]",
                                      ", co2 = ", df_ref$co2, " [ppm]",
                                      ", ppfd = ", df_ref$ppfd, " [mol/m2/s]",
                                      ", patm = ", df_ref$patm, " [Pa]",
                                      ", phi = ", df_ref$kphio, " [-]")) +
                theme(plot.caption = element_text(hjust = 0.5, size = 10, lineheight = 1.1),
                      legend.position = c(0.825, 0.25))
            
            ## Fix yscale from 0 to 1 for comparisons:
            if (scales_to_one) {
                p_append <- p_append +
                    ylim(0, 1) +
                    ylab("Relative to maximum value") +
                    geom_text(data = df_max, aes(x = x, y = 0.1, label = label), colour = "black",
                              hjust = 0,
                              size = 4,
                              inherit.aes = FALSE, parse = FALSE)
                
                ## Append plot to list
                p_all[[length(p_all)+1]] <- p_append
            }
        }
    } else if (p_output == "per_target") {
        ## Per per_target ...........................................................................
        for (t in targets) {
            
            ## Create base ggplot:
            p_append <- df_temp %>%
                pivot_longer(cols = all_of(drivers), names_to = "driver", values_to = "driver_value") %>%
                pivot_longer(cols = all_of(t), names_to = "targets", values_to = "targets_value") %>% 
                dplyr::filter(targets == t,
                              driver  == env_var) %>% 
                ggplot()
            
            ## Add vertical lines if only one formulation displayed:
            if (length(unique(df_full$formulation)) == 1) {
                p_append <- p_append +
                    geom_vline(data = df_ref %>% pivot_longer(cols = all_of(drivers), names_to = "driver", values_to = "x"),
                               aes(xintercept = x), alpha = 0.75, color = "grey50", linetype = "dotted", size = 0.75) +
                    labs(caption = paste0("\n Reference conditions (vertical grey lines): \n",
                                          "tc = ", df_ref$tc, " [°C]",
                                          ", tc_home = ", df_ref$tc_home, " [°C]",
                                          ", vpd = ", df_ref$vpd, " [Pa]",
                                          ", co2 = ", df_ref$co2, " [ppm]",
                                          ", ppfd = ", df_ref$ppfd, " [mol/m2/s]",
                                          ", patm = ", df_ref$patm, " [Pa]",
                                          ", phi = ", round(df_ref$kphio, 2), " [-]"))
            } else {
                p_append <- p_append + labs(caption = paste0("\n Reference conditions (vertical grey lines):\n",
                                                             "tc = ", df_ref$tc, " [°C]",
                                                             ", tc_home = ", df_ref$tc_home, " [°C]",
                                                             ", vpd = ", df_ref$vpd, " [Pa]",
                                                             ", co2 = ", df_ref$co2, " [ppm]",
                                                             ", ppfd = ", df_ref$ppfd, " [µmol /m2 / s]",
                                                             ", patm = ", df_ref$patm, " [Pa]"))
            }
            
            ## Add aesthetics
            p_append <- p_append +
                aes(x = driver_value, y = targets_value, color = formulation, shape = model, linetype = model) +
                # geom_point(size = 1, alpha = 0.6) +
                # scale_shape_manual(values = c("Farquhar" = "orange", "smitz37" = "brown")) +
                geom_line(size = 1, alpha = 0.75) +
                scale_color_manual(name = bquote(J[max] ~ "- Lim.:"),
                                   values = c("Farquhar" = "#00A1A0", "Smith" = "#FF8742")) +
                scale_linetype_manual(name = "Model: ",
                                      values = c("Analytical" = 6, "Numerical" = 1)) +
                facet_wrap(~driver, ncol = length(drivers), scales = "free", labeller = vlabeller) +
                
                ## Add text
                xlab(paste("Input Variable")) +
                ylab(paste(t)) +
                labs(title = paste0("Sensitivity of ", t, ", using Jmax formulation by ", ifelse(settings$rpmodel_accl$method_jmaxlim == "smith37", "Smith", "Farquhar")))+
                theme(axis.text.x = element_text(angle = 45, hjust=  1))
            # theme_bw() +
            # theme(plot.caption = element_text(hjust = 0.5, size = 10, lineheight = 1.1),
            #       strip.text.x = element_text(size = 10))
            
            ## Fixing scales for comparison
            # if ( t == "chi")     {p_append <- p_append + ylim(0.25, 1)}
            # if ( t == "vcmax")   {p_append <- p_append + ylim(0, 6e-4)}
            # if ( t == "jmax")    {p_append <- p_append + ylim(0, 6e-4)}
            # if ( t == "gs")      {p_append <- p_append + ylim(0, 8e-6)}
            # if ( t == "a_gross") {p_append <- p_append + ylim(0, 7e-5)}
            
            ## Append plot to list
            p_all[[length(p_all)+1]] <- p_append
        }
    }
    
    ## Define output ####
    out <- list(plot = p_all,
                data = df_full)   
    
    return(out)
}

## .................................................................................................
## > COST FUNCTION #####



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
    method_eb <- "tealeaves"
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
