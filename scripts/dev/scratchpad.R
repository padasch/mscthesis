¨# 03/08/2021 ####
# Attach vcmax25_pft from Kumarathunge2019
metainfo_k19 <- read_delim("~/data/mscthesis/final/metainfo_k19_2021-08-03.csv", ";", escape_double = FALSE, trim_ws = TRUE)
# write_csv(metainfo_k19, "~/data/mscthesis/final/metainfo_k19_2021-08-03.csv")

df_pft <- read.csv("~/data/photom/Tables/pft_specific_vcmax.csv") %>%
              mutate(vcmax25_pft = Vcmax25.mean * 10^-6, # Turning mumol/m2/s into mumol/m2/s
                     pft_alt = PFT,
                     pft = PFT,
                     pft = ifelse(pft == "BDT_TE", "DA-Te", pft),
                     pft = ifelse(pft == "BET_TE", "EA-Te", pft),
                     pft = ifelse(pft == "NET_TE", "EG-Te", pft),
                     pft = ifelse(pft == "NET_B",  "EG-Br", pft),
                     pft = ifelse(pft == "TET_TE", "EA-Tr", pft),
                     pft = ifelse(pft == "TUNDRA", "Arctic tundra", pft)) %>% 
        as_tibble() %>%
        dplyr::select(pft, pft_alt, vcmax25_pft)

metainfo_k19 <- metainfo_k19 %>% left_join(df_pft)
    
write_csv(metainfo_k19, "~/data/mscthesis/final/metainfo_k19.csv")

metainfo_k19 <- read_csv("~/data/mscthesis/final/metainfo_k19.csv")

# ......................................................................... ####
# 04/08/2021 ####
## Trying to get lm output into plots: ####
lm_eqn <- function(df){
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(unname(coef(m)[1]), digits = 2),
                          b = format(unname(coef(m)[2]), digits = 2),
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

p1 <- p + geom_text(x = 25, y = 300, label = lm_eqn(df), parse = TRUE)

# ......................................................................... ####
# 05/08/2021 ####
# Backup:

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

# ......................................................................... ####
# 06/08/2021 ####
## Fix optimization of jmax ####
## Get data:
df_drivers_p21    <- readRDS("~/data/mscthesis/final/df_drivers_p21.rds")
df_evaluation_p21 <- readRDS("~/data/mscthesis/raw/leaf_traits_peng2021/df_P21_clean.rds") %>%
    dplyr::select(-author) %>%
    rename_with(tolower) %>%
    dplyr::filter(!((is.na(jmax) | is.na(vcmax))),
                  !is.na(date),
                  !is.na(lat),
                  !is.na(lon)) %>%
    group_by(site, date) %>%
    nest() %>%
    ungroup() %>%
    rename(sitename = site,
           data_raw = data)

## Problem creation: jmax predictions are super off in numerical model:
settings <- get_settings()
settings$rpmodel_accl$method_optim <- "numerical"
df_acc   <- run_rpmodel_accl(settings, df_drivers_p21, df_evaluation_p21) 
df_x     <- df_acc %>%
    mutate(tc_growth_air = purrr::map(forcing, ~ dplyr::select(., tc_growth_air)),
           tc_home = purrr::map(forcing, ~ dplyr::select(., tc_home))) %>% 
    unnest(c(tc_growth_air, tc_home)) %>% 
    mutate(tc_leaf   = purrr::map_dbl(data_raw, ~ mean(.$tleaf)),
           jmax      = purrr::map_dbl(data_raw, ~ mean(.$jmax))  * 10 ^-6,
           vcmax     = purrr::map_dbl(data_raw, ~ mean(.$vcmax)) * 10 ^-6,
           jmax25    = jmax  / calc_ftemp_inst_jmax(tc_leaf, tc_growth_air, tc_home, method_ftemp = settings$rpmodel_accl$method_ftemp),
           vcmax25   = vcmax / calc_ftemp_inst_vcmax(tc_leaf, tc_growth_air, method_ftemp = settings$rpmodel_accl$method_ftemp))

jmax_obs <- df_x$jmax
jmax_rpm <- df_x %>% dplyr::select(-jmax, -jmax25, -vcmax, -vcmax25) %>% unnest(rpm_accl) %>% pull(jmax)

plot(jmax_rpm, jmax_obs, xlim = c(0, 6e-4), ylim = c(0, 6e-4))
abline(0, 1)

## Debugging optimization:
## Get input data:
df_1 <- df_x %>% dplyr::select(-tc_growth_air, -tc_home) %>% unnest(forcing) %>% slice(10)

tc_leaf <- df_1$tc_growth_air
patm <- df_1$patm_growth
co2 <- df_1$co2_growth
vpd <- df_1$vpd_growth
ppfd <- df_1$ppfd_growth * 3600 * 24
fapar <- df_1$fapar
kphio <- settings$rpmodel_accl$kphio_calib
beta <- 146
c_cost <- 0
vcmax_start <- df_x$vcmax %>% mean() * 3600 * 24
jmax_start <- df_x$jmax %>% mean() * 3600 * 24
gs_start <- 0.5 
method_jmaxlim_inst <- "smith37"

## Testing different optimizer
optimr::optimr(
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
    control    = list( maxit = 1000000000 )
)

optim(par        = c( vcmax_start,       gs_start,       jmax_start ), # starting values
      lower      = c( vcmax_start*0.001, gs_start*0.001, jmax_start*0.001 ),
      upper      = c( vcmax_start*1000,  gs_start*1000,  jmax_start*1000 ),
      fn         = optimise_this_tcleaf_vcmax_jmax,
      args       = c(tc_leaf, patm, co2, vpd),
      iabs       = (ppfd * fapar),
      kphio      = kphio,
      beta       = 146,
      c_cost     = 0,
      method_jmaxlim_inst = method_jmaxlim_inst,
      method     = "L-BFGS-B"
      # gr = NULL,
      # …,
      # method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"),
      # lower = -Inf,
      # upper = Inf,
      # control = list(),
      # hessian = FALSE
      )


(dat_1 <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf =tc_leaf,
                                        patm =patm,
                                        co2 =co2,
                                        vpd =vpd,
                                        ppfd =ppfd,
                                        fapar =fapar,
                                        kphio =kphio,
                                        beta =beta,
                                        c_cost =c_cost,
                                        vcmax_start =vcmax_start,
                                        gs_start =gs_start,
                                        jmax_start =jmax_start,
                                        method_jmaxlim_inst =method_jmaxlim_inst))


## Comparisons:
cat("Param | Org.| Opt. \n",
    "vcmax | ", vcmax_start, dat_1$vcmax_mine, "\n",
    "jmax | ",  jmax_start, dat_1$jmax_mine, "\n")

### Original function by beni: ####


(df_beni <- calc_optimal_gs_vcmax_jmax(
    kmm = calc_kmm(tc_leaf, patm),
    gammastar =calc_gammastar(tc_leaf, patm),
    ns_star=calc_density_h2o(tc_leaf, patm) / calc_density_h2o(25, 101325),
    ca=co2_to_ca(co2, patm),
    vpd=vpd,
    ppfd=ppfd,
    fapar=fapar,
    kphio=kphio,
    beta=beta,
    c_cost=c_cost,
    vcmax_start=vcmax_start,
    gs_start=gs_start,
    jmax_start=jmax_start,
    method_jmaxlim_inst = method_jmaxlim_inst
))

df_beni$par /3600/24 * 10^4





calc_optimal_gs_vcmax_jmax <- function( kmm, gammastar, ns_star, ca, vpd, ppfd, fapar, kphio, beta, c_cost, vcmax_start, gs_start, jmax_start, method_jmaxlim_inst ){
    
    optimise_this_gs_vcmax_jmax <- function( par, args, iabs, kphio, beta, c_cost, method_jmaxlim_inst, maximize=FALSE, return_all=FALSE ){
        
        kmm       <- args[1]
        gammastar <- args[2]
        ns_star   <- args[3]
        ca        <- args[4]
        vpd       <- args[5]
        
        vcmax <- par[1]
        gs    <- par[2]
        jmax  <- par[3]
        
        ## Electron transport is limiting
        ## Solve quadratic equation system using: A(Fick's Law) = A(Jmax Limitation)
        ## This leads to a quadratic equation:
        ## A * ci^2 + B * ci + C  = 0
        ## 0 = a + b*x + c*x^2
        
        ## Jmax Limitation following Smith (1937):
        if (method_jmaxlim_inst == "smith37") {
            ## A = gs (ca - ci)
            ## A = kphio * iabs (ci-gammastar)/ci+2*gammastar) * L
            ## L = 1 / sqrt(1 + ((4 * kphio * iabs)/jmax)^2)
            
            ## with
            L <- 1.0 / sqrt(1.0 + ((4.0 * kphio * iabs)/jmax)^2)
            A <- -gs
            B <- gs * ca - 2 * gammastar * gs - L * kphio * iabs
            C <- 2 * gammastar * gs * ca + L * kphio * iabs * gammastar
            
            ci_j <- QUADM(A, B, C)
            a_j  <- kphio * iabs * (ci_j - gammastar)/(ci_j + 2 * gammastar) * L  
        }
        
        ## Jmax Limitation following Farquhar (1989):
        if (method_jmaxlim_inst == "farquhar89") {
            ## A = gs (ca - ci)
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
        # 
        if (return_all){
            return(tibble(vcmax_mine = vcmax, jmax_mine = jmax, gs_mine = gs, ci_mine = ci, chi_mine = ci / ca, a_c_mine = a_c, a_j_mine = a_j, assim = assim, ci_c_mine = ci_c, ci_j_mine = ci_j, cost_transp = cost_transp, cost_vcmax = cost_vcmax, cost_jmax = cost_jmax, net_assim = net_assim, method_jmaxlim_inst = method_jmaxlim_inst))
        } else {
            return( net_assim )
        }
    }
    
    out_optim <- optimr::optimr(
        par        = c( vcmax_start, gs_start, jmax_start ), # starting values
        lower      = c( vcmax_start*0.001, gs_start*0.001, jmax_start*0.001 ),
        upper      = c( vcmax_start*1000, gs_start*1000, jmax_start*1000 ),
        fn         = optimise_this_gs_vcmax_jmax,
        args       = c(kmm, gammastar, ns_star, ca, vpd),
        iabs       = (ppfd * fapar),
        kphio      = kphio,
        beta       = beta,
        c_cost     = c_cost/4,
        method_jmaxlim_inst = method_jmaxlim_inst,
        method     = "L-BFGS-B",
        maximize   = TRUE,
        control    = list( maxit = 100000 )
    )
    
    varlist <- optimise_this_gs_vcmax_jmax( 
        par=out_optim$par, 
        args=c(kmm, gammastar, ns_star, ca, vpd), 
        iabs=(fapar*ppfd), 
        kphio, beta, c_cost/4,
        method_jmaxlim_inst,
        maximize=FALSE, 
        return_all=TRUE 
    )
    
    # return(varlist)
    return(out_optim)
}

## other

sol_optimize <- optimize(maximize_this_tc_leaf,
                         interval  = c(1, 40),
                         tc_air    = 25,
                         tc_growth = 25,
                         tc_home   = 25,
                         ppfd      = 130/3600/24,
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
                         c_cost    = 0.41)


test$con



# ......................................................................... ####
# 08/08/2021 ####
## Trying to get that 3d plot optimization of gs vcmax and jmax ####

## Testing plotly
library(plotly)
fig <- plot_ly(z = ~volcano)
fig <- fig %>% add_surface()

fig

fig <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
                                   yaxis = list(title = 'Gross horsepower'),
                                   zaxis = list(title = '1/4 mile time')))


## Testing my own plotly
df <- tibble(x = seq(1, 100, length.out = 100),
             y = sqrt(x) * x/2,
             z <- x^3/y/10,
             f <- (x+z+y)/x)

fig <- plot_ly(df, x = ~x, y = ~y, z = ~z, color = ~f)

x <- array(rep(1, 365*5*4), dim=c(365, 5, 4))



## Hyperbolix minimum formulaion ####
-QUADP(A = 1 - 1E-07, B = a_c + a_j, C = a_c*a_j)





# ......................................................................... ####
# 09/08/2021 ####
## Subset where convergence worked vs. all data ####
## Get data
## Analytical
settings <- get_settings()
settings$rpmodel_accl$kphio_calib <- settings$rpmodel_accl$kphio_calib * 1.5
df_ana_s37   <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax()
p_ana_s37    <- plot_two_long_df(df_x = make_long_df(df_in = df_ana_s37, dataset = "analytical_s37"), df_x_dataset = "analytical_s37", df_y = make_df_evaluation_p21_long(df_in = df_ana_s37), df_y_dataset = "peng21")


## Numerical - all data
settings$rpmodel_accl$method_optim   <- "numerical"
settings$rpmodel_accl$method_eb   <- "on"
df_num_s37 <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
p_num_s37  <- plot_two_long_df(df_x = make_long_df(df_in = df_num_s37, dataset = "numerical_s37"), df_x_dataset = "numerical_s37", df_y = make_df_evaluation_p21_long(df_in = df_num_s37), df_y_dataset = "peng21")


## Numerical - only where convergence worked
names_v <- df_num_s37$rpm_accl[[1]] %>% names()
df_num_s37_x <- df_num_s37 %>% unnest(rpm_accl) %>% dplyr::filter(opt_convergence == 0) %>% nest(rpm_accl = any_of(names_v))
p_num_s37_x  <- plot_two_long_df(
  df_x = make_long_df(df_in = df_num_s37_x, dataset = "numerical_s37"),
  df_x_dataset = "numerical_s37",
  df_y = make_df_evaluation_p21_long(df_in = df_num_s37_x),
  df_y_dataset = "peng21"
)


## Plots
p_ana_s37 + xlim(0, 4e-4) + ylim(0, 4e-4) + ggtitle("Analytic")
p_num_s37 + xlim(0, 4e-4) + ylim(0, 4e-4) + ggtitle("Including non-convergence")
p_num_s37_x + xlim(0, 4e-4) + ylim(0, 4e-4) + ggtitle("Without non-convergence")

## Results
#' Taking only the subset of converged data points does not show any significant
#' change in estimated parameters vcmax and jmax, thus this check can be left-out
#' but it should still be mentioned that the numerical solutions does not always
#' converge which hampers the predictions


## Trying heatmap plot from rbeni ####
df_t <- df_ana_s37 %>%
  mutate(tc_growth_air = purrr::map(forcing, ~ dplyr::select(., tc_growth_air)),
         tc_home = purrr::map(forcing, ~ dplyr::select(., tc_home))) %>% 
  unnest(c(tc_growth_air, tc_home)) %>% 
  mutate(tc_leaf   = purrr::map_dbl(data_raw, ~ mean(.$tleaf)),
         jmax      = purrr::map_dbl(data_raw, ~ mean(.$jmax))  * 10 ^-6,
         vcmax     = purrr::map_dbl(data_raw, ~ mean(.$vcmax)) * 10 ^-6,
         jmax25    = jmax  / calc_ftemp_inst_jmax(tc_leaf, tc_growth_air, tc_home, method_ftemp = "kumarathunge19"),
         vcmax25   = vcmax / calc_ftemp_inst_vcmax(tc_leaf, tc_growth_air, method_ftemp = "kumarathunge19")) %>% 
  rename(jmax_obs = jmax,
         vcmax_obs = vcmax,
         jmax25_obs = jmax25,
         vcmax25_obs = vcmax25) %>% 
  unnest(rpm_accl)
  # unnest(rpm_accl) %>% 
  # dplyr::filter(opt_convergence == 0)

df_t %>% ggplot(aes(vcmax*10^6, vcmax_obs*10^6)) +
  geom_point() + 
  geom_smooth(method = "lm") + 
  geom_abline() + 
  xlim(0, 4e-4*10^6) + 
  ylim(0, 4e-4*10^6)

library(rbeni)
heatscatter(df_t$vcmax*10^6, df_t$vcmax_obs*10^6, xlab = "Modelled", ylab = "Observed", ggplot = TRUE)  +
  geom_abline() +
  geom_smooth(data = df_t, aes(vcmax*10^6, vcmax_obs*10^6), method= "lm", fullrange = T) +
  xlim(0, 4e-4*10^6) + 
  ylim(0, 4e-4*10^6)


## Testing performance of EB ####
## EB via plantecophys
settings <- get_settings()
settings$rpmodel_accl$method_optim <- "numerical"
settings$rpmodel_accl$energy_balance <- "on"
df_eb_pl <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)

## EB via tealeaves
settings$rpmodel_accl$method_eb <- "tealeaves"
df_eb_te <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
p_eb_te  <- plot_two_long_df( df_x = make_long_df(df_in = df_eb_te, dataset = "eb_tealeaves"), df_x_dataset = "eb_tealeaves", df_y = make_df_evaluation_p21_long(df_in = df_eb_te), df_y_dataset = "peng21")

p_eb_te + xlim(0, 4e-4) + ylim(0, 4e-4) + ggtitle("tealeaves")


## Do analytic and numeric show same chi behavior? ####

df_1 <- df_ana_s37 %>% unnest(c(rpm_accl, forcing))
df_2 <- df_num_s37 %>% unnest(c(rpm_accl, forcing))
df_3 <- df_eb_te %>% unnest(c(rpm_accl, forcing))
df_4 <- df_eb_pl %>% unnest(c(rpm_accl, forcing))


plot(df_1$chi, df_2$chi, ylim = c(0, 1), xlim = c(0,1))
abline(0, 1)

plot(df_1$chi, df_4$chi, ylim = c(0, 1), xlim = c(0,1))
abline(0, 1)
  
plot(df_1$chi, df_3$chi, ylim = c(0, 1), xlim = c(0,1))
abline(0, 1)

plot(df_4$tc_growth_air, df_4$tc_growth_leaf, ylim = c(0, 40), xlim = c(0, 40))
abline(0, 1)

plot(df_3$tc_growth_air, df_3$tc_growth_leaf, ylim = c(0, 40), xlim = c(0, 40))
abline(0, 1)

ggplot(df_3, aes(tc_growth_air, tc_growth_leaf)) + geom_smooth(method = "lm", fullrange = T) + geom_point() + geom_abline() + ylim(0,40) + xlim(0,40)


# ......................................................................... ####
# 10/08/2021 ####
## Testing tealeaves ####
leaf_par   <- make_leafpar()   # leaf parameters
enviro_par <- make_enviropar() # environmental parameters
constants  <- make_constants() # physical constants

tair_v <- seq(0, 50, 0.5)
tleaf_v <- rep(NA, length(tair_v))

for (t in 1:length(tair_v)) {
  enviro_par <- make_enviropar(replace = list( T_air = set_units(tair_v[t] + 273.15, "K")))
  tleaf_v[t] <- tleaf(leaf_par, enviro_par, constants, quiet = TRUE)$T_leaf %>% set_units("degree_Celsius") %>% drop_units() %>% round(1)
}

plot(tair_v, tleaf_v, ylim = c(0, 40), xlim = c(0, 40))
abline(0, 1)

plot(tair_v, tleaf_v-tair_v, ylim = c(-20, 20), xlim = c(0, 40))
abline(0, 0)

tair_v[which(tair_v == tleaf_v)]


