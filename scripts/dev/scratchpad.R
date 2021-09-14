# setup packages ####
suppressMessages(source("~/projects/mscthesis/scripts/final/source_fct_pkg.R"))

# ......................................................................... ####
# 03/08/2021 ####
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
df_ana_s37   <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(df_in = ., ftemp_method = settings$rpmodel_accl$method_ftemp)
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


# ......................................................................... ####
# 14/08/2021 ####

tc_air <- 8

## Trying to fix calc_tc_leaf_from_eb() ####
### Optimr() - works fine for "sqrd" and "abs", "diff" & maximize = F, leads to: control$fnscale and control$maximize conflict

optimr::optimr(
  # Parameter boundaries to optimize within:
  par       = tc_air,
  lower     = tc_air - 30,
  upper     = tc_air + 30,
  
  # Function to optimize and its inputs:
  fn        = LeafEnergyBalance,
  Tair      = tc_air,    # input in degC
  returnwhat= "sqrd",
  
  # Optimr settings:
  method    = "L-BFGS-B",
  control   = list(maxit = 100)
)

### uniroot - works fine
uniroot(LeafEnergyBalance,
        interval = c(tc_air-30, tc_air+30), 
        Tair = tc_air,
        returnwhat = "diff")


### optimize - works fine for "abs" (slower) and "sqrd"
optimize(LeafEnergyBalance,
         interval = c(tc_air-30, tc_air+30), 
         Tair = tc_air,
         returnwhat = "sqrd")



## Trying to fix maximize_this_tc_leaf() ####
### optimize - works fine for "abs" (slower) and "sqrd" and both eb
optimize(maximize_this_tc_leaf,
         interval  = c(tc_air-30, tc_air+30),
         tc_air    = tc_air,
         method_opt = "sqrd",
         method_eb  = "plantecophys")

### uniroot - works fine for both eb
uniroot(maximize_this_tc_leaf,
        interval  = c(tc_air-30, tc_air+30),
        tc_air    = tc_air,
        method_opt = "diff",
        method_eb  = "plantecophys")

### optimr - "abs" and "sqrd" crashes, "diff" & maximize = F, leads to: control$fnscale and control$maximize conflict
optimr(fn = maximize_this_tc_leaf,
       par = tc_air,
       lower = tc_air - 30,
       upper = tc_air + 30,
       # control = list(maximize = T),
       method = "L-BFGS-B",
       method_opt = "abs",
       method_eb = "plantecophys")


# ......................................................................... ####
# 17/08/2020 ####

## Why does numerical farquhar fuck up at below 10 °C? ####
v1 <- seq(0, 40, length.out = 40)
v2 <- rep(NA, length(v1))

for (v in 1:length(v1)) {
  v2[v] <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = v1[v])$varlist$jmax_mine
}

plot(v1, v2)

# ......................................................................... ####
# 19/08/2020 ####

## Loop to find best starting conditions #####

## Fixed starting values
vcmax_fix <- 15
jmax_fix  <- 30
gs_fix    <- 0.5

## Vectors of starting values
vcmax_v <- c(1, 2, 6, 10, 15, 30, 50)
jmax_v  <- c(1, 2, 6, 10, 15, 30, 50)*2
gs_v    <- c(1, 2, 6, 10, 15, 30, 50)/10

df_start <- bind_rows(list =
                        tibble(vcmax_start = vcmax_v,
                               jmax_start  = rep(jmax_fix, length(vcmax_v)),
                               gs_start    = rep(gs_fix,   length(vcmax_v)),
                               opt_var     = "vcmax_start"),
                      tibble(vcmax_start = rep(vcmax_fix, length(vcmax_v)),
                             jmax_start  = jmax_v,
                             gs_start    = rep(gs_fix,   length(vcmax_v)),
                             opt_var     = "jmax_start"),
                      tibble(vcmax_start = rep(vcmax_fix, length(vcmax_v)),
                             jmax_start  = rep(jmax_fix, length(vcmax_v)),
                             gs_start    = gs_v,
                             opt_var     = "gs_start"))

df_start$data_n_s37 <- NA
df_start$data_n_f89 <- NA

p_list <- list()

for (i in 1:nrow(df_start)) {

  
  cat(i, "/", nrow(df_start), "\n")
  calc_optimal_gs_vcmax_jmax <- function(tc_air, # Input in: degC
                                         patm, # Input in: Pa
                                         co2, # Input in: Pa
                                         vpd, # Input in: Pa
                                         ppfd , # Input in: mol/m2/s
                                         fapar, # Input in: -
                                         kphio, # Input in: -
                                         beta = 146, # Input in: -
                                         vcmax_start = df_start$vcmax_start[i],# 50e-6  * 24 * 3600,  # Input in: mol/m2/h
                                         jmax_start  = df_start$jmax_start[i],# 100e-6 * 24 * 3600,  # Input in: mol/m2/h
                                         gs_start    = df_start$gs_start[i],# 5e-6  * 24 * 3600,  # Input in: mol/m2/h
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
    
    ## Test Case: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    calc_optimal_gs_vcmax_jmax(tc_air = 10,
                               patm = 101325,
                               co2 = 400,
                               vpd = 1000,
                               ppfd = 800e-6,
                               fapar = 1,
                               kphio = 0.05,
                               method_jmaxlim_inst = "smith37",
                               method_eb = "plantecophys",
                               which_optimizer = "optimr")
  }
  
  ## Comparing Numerical Outputs
  ## Smith:
  settings <- get_settings()
  settings$verbose <- F
  settings$rpmodel_accl$method_optim  <- "numerical"
  df_num_s37 <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
  p_num_s37  <-
    plot_two_long_df(
      df_x = make_long_df(df_in = df_num_s37, dataset = "numerical_s37"),
      df_x_dataset = "numerical_s37",
      df_y = make_df_evaluation_p21_long(df_in = df_num_s37),
      df_y_dataset = "peng21",
      kphio_vcmax_corr = F
    )
  
  ## Farquhar:
  settings <- get_settings()
  settings$verbose <- F
  settings$rpmodel_accl$method_jmaxlim <- "farquhar89"
  settings$rpmodel_accl$method_optim <- "numerical"
  df_num_f89 <-
    run_rpmodel_accl(settings = settings,
                     df_drivers = df_drivers_p21,
                     df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(ftemp_method = settings$rpmodel_accl$method_ftemp)
  p_num_f89  <-
    plot_two_long_df(
      df_x = make_long_df(df_in = df_num_f89, dataset = "numerical_f89"),
      df_x_dataset = "numerical_f89",
      df_y = make_df_evaluation_p21_long(df_in = df_num_f89),
      df_y_dataset = "peng21",
      kphio_vcmax_corr = F
    )
  
  ## Save data
  df_start$data_n_s37[i] <- list(p_num_s37$data %>% dplyr::filter(variable == "vcmax"))
  df_start$data_n_f89[i] <- list(p_num_f89$data %>% dplyr::filter(variable == "vcmax"))
  
  
  ## Plots
  p_num_s37 <- p_num_s37 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + theme(legend.position = "bottom") + ggtitle("Smith-Formulation")
  p_num_f89 <- p_num_f89 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + theme(legend.position = "bottom") + ggtitle("Farquhar-Formulation")
  p_new     <- p_num_s37 / p_num_f89
  
  p_list[[length(p_list)+1]] <- p_new
  
}

## Get metrics:
df_startingvalues <- df_start %>%
  mutate(sry_n_s37 = purrr::map(data_n_s37, ~summary(lm(y ~ x, data = .))),
         sry_n_f89 = purrr::map(data_n_f89, ~summary(lm(y ~ x, data = .))),
         r2_n_s37   = purrr::map_dbl(sry_n_s37, ~ .$r.squared),
         r2_n_f89   = purrr::map_dbl(sry_n_f89, ~ .$r.squared),
         rmse_n_s37 = purrr::map_dbl(sry_n_s37, ~ sqrt( ( c(crossprod(.$residuals)) / length(.$residuals) ) )),
         rmse_n_f89 = purrr::map_dbl(sry_n_f89, ~ sqrt( ( c(crossprod(.$residuals)) / length(.$residuals) ) )))

df_startingvalues_without_dfs <- df_startingvalues %>% dplyr::select(-c("data_n_s37"  ,"data_n_f89" , "sry_n_s37" ,"sry_n_f89"))

p_startingvalues <- p_list

# ......................................................................... ####
# 20/08/2020 ####

## Loop to find best acclimation timespan #####
vec <- c(7, 15, 30, 60, 90, 120)

df_start <- bind_rows(list =
                        tibble(tau_temp = vec,
                               tau_ppfd = rep(30, length(vec)),
                               opt_var  = "tau_temp"),
                      tibble(tau_temp = rep(30, length(vec)),
                             tau_ppfd = vec,
                             opt_var  = "tau_ppfd"))

df_start$data_a_s37 <- NA
df_start$data_a_f89 <- NA
df_start$data_n_s37 <- NA
df_start$data_n_f89 <- NA


p_list <- list()

for (i in 1:nrow(df_start)) {
  
  cat(i, "/", nrow(df_start), "\n")
  
  ## Analytical models:
  ## Smith:
  settings <- get_settings()
  settings$verbose <- F
  settings$rpmodel_accl$tau$tc_air <- df_start$tau_temp[i]
  settings$rpmodel_accl$tau$ppfd   <- df_start$tau_ppfd[i]
  df_ana_s37   <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(df_in = ., ftemp_method = settings$rpmodel_accl$method_ftemp)
  p_ana_s37    <- plot_two_long_df( df_x = make_long_df(df_in = df_ana_s37, dataset = "analytical_s37", v_vars = "vcmax"), df_x_dataset = "analytical_s37", df_y = make_df_evaluation_p21_long(df_in = df_ana_s37), df_y_dataset = "peng21", kphio_vcmax_corr = F)
  
  ## Farquhar:
  settings$rpmodel_accl$method_jmaxlim <- "farquhar89"
  df_ana_f89   <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(df_in = ., ftemp_method = settings$rpmodel_accl$method_ftemp)
  p_ana_f89    <- plot_two_long_df( df_x = make_long_df(df_in = df_ana_f89, dataset = "analytical_f89", v_vars = "vcmax"), df_x_dataset = "analytical_f89", df_y = make_df_evaluation_p21_long(df_in = df_ana_f89), df_y_dataset = "peng21", kphio_vcmax_corr = F)
  
  
  ## Numerical Models:
  ## Smith:
  settings$rpmodel_accl$method_jmaxlim <- "smith37"
  settings$rpmodel_accl$method_optim   <- "numerical"
  df_num_s37 <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(., ftemp_method = settings$rpmodel_accl$method_ftemp)
  p_num_s37  <- plot_two_long_df( df_x = make_long_df(df_in = df_num_s37, dataset = "numerical_s37", v_vars = "vcmax"), df_x_dataset = "numerical_s37", df_y = make_df_evaluation_p21_long(df_in = df_num_s37), df_y_dataset = "peng21", kphio_vcmax_corr = F)
  
  ## Farquhar:
  settings$rpmodel_accl$method_jmaxlim <- "farquhar89"
  settings$rpmodel_accl$method_optim   <- "numerical"
  df_num_f89 <- run_rpmodel_accl(settings = settings, df_drivers = df_drivers_p21, df_evaluation = df_evaluation_p21) %>% get_instant_vcmax_jmax(., ftemp_method = settings$rpmodel_accl$method_ftemp)
  p_num_f89  <- plot_two_long_df( df_x = make_long_df(df_in = df_num_f89, dataset = "numerical_f89", v_vars = "vcmax"), df_x_dataset = "numerical_f89", df_y = make_df_evaluation_p21_long(df_in = df_num_f89), df_y_dataset = "peng21", kphio_vcmax_corr = F)
  
  
  ## Save data:
  df_start$data_a_s37[i] <- list(p_ana_s37$data %>% dplyr::filter(variable == "vcmax"))
  df_start$data_a_f89[i] <- list(p_ana_f89$data %>% dplyr::filter(variable == "vcmax"))
  df_start$data_n_s37[i] <- list(p_num_s37$data %>% dplyr::filter(variable == "vcmax"))
  df_start$data_n_f89[i] <- list(p_num_f89$data %>% dplyr::filter(variable == "vcmax"))
  
  
  ## Plots
  p_ana_s37 <- p_ana_s37 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + ggtitle("Smith-Formulation - Analytic")
  p_ana_f89 <- p_ana_f89 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + ggtitle("Farquhar-Formulation - Analytic")
  p_num_s37 <- p_num_s37 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + ggtitle("Smith-Formulation - Numeric")
  p_num_f89 <- p_num_f89 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + ggtitle("Farquhar-Formulation - Numeric")
  
  p_new     <- p_num_s37 + p_num_f89 + p_ana_s37 + p_ana_f89
  
  p_list[[length(p_list)+1]] <- p_new
  
}

## Get model metrics:
df_timespan <- df_start %>%
  mutate(sry_a_s37 = purrr::map(data_a_s37, ~summary(lm(y ~ x, data = .))),
         sry_a_f89 = purrr::map(data_a_f89, ~summary(lm(y ~ x, data = .))),
         sry_n_s37 = purrr::map(data_n_s37, ~summary(lm(y ~ x, data = .))),
         sry_n_f89 = purrr::map(data_n_f89, ~summary(lm(y ~ x, data = .))),
         r2_a_s37   = purrr::map_dbl(sry_a_s37, ~ .$r.squared),
         r2_a_f89   = purrr::map_dbl(sry_a_f89, ~ .$r.squared),
         r2_n_s37   = purrr::map_dbl(sry_n_s37, ~ .$r.squared),
         r2_n_f89   = purrr::map_dbl(sry_n_f89, ~ .$r.squared),
         rmse_a_s37 = purrr::map_dbl(sry_a_s37, ~ sqrt( ( c(crossprod(.$residuals)) / length(.$residuals) ) )),
         rmse_a_f89 = purrr::map_dbl(sry_a_f89, ~ sqrt( ( c(crossprod(.$residuals)) / length(.$residuals) ) )),
         rmse_n_s37 = purrr::map_dbl(sry_n_s37, ~ sqrt( ( c(crossprod(.$residuals)) / length(.$residuals) ) )),
         rmse_n_f89 = purrr::map_dbl(sry_n_f89, ~ sqrt( ( c(crossprod(.$residuals)) / length(.$residuals) ) )))

p_timespan <- p_list
df_timespan_without_dfs <- df_timespan %>% dplyr::select(-c("data_a_s37" ,  "data_a_f89"  , "data_n_s37"  ,"data_n_f89" , "sry_a_s37" ,"sry_a_f89" ,"sry_n_s37" ,"sry_n_f89"))



## .....................................................................................................................
## Global Maps for Peng 2021 #### 
library(ggmap)
library(ggrepel)
# Global Maps Basis
kg <- readOGR(dsn = "~/data/climate_zones/1976-2000_GIS", layer = "1976-2000")
kg <- spTransform(kg, CRS("+proj=longlat +datum=WGS84"))
kg_f <- fortify(kg, region = "GRIDCODE")
key <- data.frame( id = c(11, 12, 13, 14, 21, 22, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 61, 62),
                   cl = c('Af', 'Am', 'As', 'Aw', 'BSh', 'BSk', 'BWh', 'BWk', 'Cfa', 'Cfb','Cfc', 'Csa', 'Csb', 'Csc', 'Cwa','Cwb', 'Cwc', 'Dfa', 'Dfb', 'Dfc','Dfd', 'Dsa', 'Dsb', 'Dsc', 'Dsd','Dwa', 'Dwb', 'Dwc', 'Dwd', 'EF','ET', 'Ocean'))
kg_final <- merge(kg_f, key)

p_base <- ggplot() +
  geom_polygon(data = kg_final, aes(x = long, y = lat, group = group, fill = cl), alpha = 0.8) +
  ylab("Latitude (Decimal Degree)") +
  xlab("Longitude (Decimal Degree)") +
  coord_cartesian(ylim=c(-80, 80)) +
  scale_fill_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC", "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64", "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00", "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF", "#64FFFF", "#F5FFFF")) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 12, direction = "horizontal", title.position = "top")) +
  labs(fill = "Köppen-Geiger Climate Zone")

## P21 Global Map
# Get df_evaluation_p21 from analysis of acclimate pmodel

siteinfo <- df_evaluation_p21 %>% mutate(lon = purrr::map_dbl(data_raw, . %>% pull(lon) %>% unique()),
                                         lat = purrr::map_dbl(data_raw, . %>% pull(lat) %>% unique()))

p <- p_base +
  geom_point(data = siteinfo, aes(x = lon, y = lat), color = "black", size = 3, pch = 21, fill = "white") +    
  ggtitle(bquote("Global distribution of " ~ V[cmax] ~ "observations"))

# ggsave("~/projects/mscthesis/docs/fig-global-map-p21.pdf", p, height = 7, width = 8)


## Attach KG Climate Zones
df_1 <- siteinfo
coordinates(df_1) <- ~ lon + lat
crs(df_1) <- crs(kg)
df_2 <- point.in.poly(df_1, kg, sp = TRUE, duplicate = TRUE)
df_2 <- as.data.frame(df_2)
df_2$id <- df_2$GRIDCODE
df_2 <- merge(df_2, key)
siteinfo_kg <- df_2 %>% dplyr::select(cl, sitename)


df_evaluation_p21_new <- df_evaluation_p21 %>%
  left_join(siteinfo_kg %>% rename(site = sitename)) %>% 
  distinct()

# saveRDS(df_evaluation_p21_new, "~/data/mscthesis/raw/leaf_traits_peng2021/df_P21_clean_with_cl.rds")

## .................................................................................................
## Global Map für K 19 ####
# K19 Data
siteinfo <- readRDS("~/Polybox/2_ETH/ONGOING/msc-thesis/rproject/data/mscthesis/siteinfo_for_cluster.rds") %>%
  left_join(read_csv("~/data/mscthesis/final/metainfo_k19.csv"))

p <- p_base +
  geom_point(data = siteinfo, aes(x = lon, y = lat), color = "black", size = 3, pch = 21, fill = "white") +    
  ggtitle(bquote("Global distribution of " ~ T[opt] ~ "observations")) 
  # geom_label_repel(data = siteinfo,
  #                  aes(x = lon, y = lat, label = source),
  #                  box.padding = 0.5,
  #                  point.padding = 0.5,
  #                  min.segment.length = 0, # draw all line segments
  #                  arrow = arrow(length = unit(0.015, "npc")))

# ggsave("~/projects/mscthesis/docs/fig-global-map-k19.pdf", p, height = 7, width = 8)

## Attaching Climate Zone to K19 Metadata ####
metainfo_k19 <- read_csv("~/data/mscthesis/final/metainfo_k19.csv")

df_1 <- siteinfo
coordinates(df_1) <- ~ lon + lat
crs(df_1) <- crs(kg)
df_2 <- point.in.poly(df_1, kg, sp = TRUE, duplicate = TRUE)
df_2 <- as.data.frame(df_2)
df_2$id <- df_2$GRIDCODE
df_2 <- merge(df_2, key)
siteinfo_kg <- df_2 %>% dplyr::select(cl, sitename)


metainfo_k19_new <- metainfo_k19 %>%
  left_join(siteinfo_kg %>% rename(dataset = sitename, climate_zone = cl)) %>% 
  distinct()

write_csv(metainfo_k19_new, "~/data/mscthesis/final/metainfo_k19.csv")



 # ......................................................................... ####
# 04/09/2020 ####
## Plots for thesis from scratch ####
## Prerequisites
ftemp_method <- settings$rpmodel_accl$method_ftemp
standard_error <- function(x) sd(x, na.rm = T) / sqrt(length(x)) # Create own se function
df_in <- df_ana_s37

df_temp <- df_in %>%
  mutate(tc_growth_air = purrr::map(forcing, ~ dplyr::select(., tc_growth_air)),
         tc_home = purrr::map(forcing, ~ dplyr::select(., tc_home))) %>% 
  unnest(c(tc_growth_air, tc_home)) %>% 
  mutate(
         # P21 Data
         p21_tc_leaf     = purrr::map_dbl(data_raw, ~ mean(.$tleaf, na.rm = T)),
         p21_jmax_mean   = purrr::map_dbl(data_raw, ~ mean(.$jmax, na.rm = T))  * 10 ^-6,
         p21_jmax_se     = purrr::map_dbl(data_raw, ~ standard_error(.$jmax)) * 10 ^-6,
         p21_vcmax_mean  = purrr::map_dbl(data_raw, ~ mean(.$vcmax, na.rm = T)) * 10 ^-6,
         p21_vcmax_se    = purrr::map_dbl(data_raw, ~ standard_error(.$vcmax)) * 10 ^-6,
         
         # Rpmodel predictions
         rpm_jmax        = purrr::map_dbl(rpm_accl, ~magrittr::extract(.$jmax25)),
         rpm_vcmax       = purrr::map_dbl(rpm_accl, ~magrittr::extract(.$vcmax25)),
         rpm_jmax        = rpm_jmax  * calc_ftemp_inst_jmax(p21_tc_leaf, tc_growth_air, tc_home, method_ftemp = ftemp_method),
         rpm_vcmax       = rpm_vcmax * calc_ftemp_inst_vcmax(p21_tc_leaf, tc_growth_air, method_ftemp = ftemp_method),
         
         # Climate Zones
         climate_zone    = purrr::map_chr(data_raw, ~ pull(., cl) %>% unique()),
         climate_zone    = as.factor(climate_zone)) 


## Vcmax plot
df_temp %>% 
  ggplot() +
  aes(x = rpm_vcmax * 10^6,
      y = p21_vcmax_mean * 10^6) +
  geom_point(aes(color = climate_zone)) + 
  geom_errorbar(aes(ymin=(p21_vcmax_mean-p21_vcmax_se)*1e6, ymax=(p21_vcmax_mean+p21_vcmax_se)*1e6, color = climate_zone), width=0.25) +
  geom_smooth(method = "lm", fullrange = T) +
  geom_abline() +
  ggpmisc::stat_poly_eq(data = df_temp,
                        formula = y ~ x,
                        method = "lm",
                        aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
                        parse = TRUE) +
  xlab("Observation") +
  ylab("Prediction") +
  ylim(0, 250) +
  xlim(0, 250) +
  scale_color_manual(name = c( "Af", "Am", "As", "Aw",
                              "BWk", "BWh", "BSk", "BSh",
                              "Cfa", "Cfb", "Cfc", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc",
                              "Dfa", "Dfb", "Dfc", "Dfd", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd",
                              "EF", "ET"),
                    values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                               "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                               "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                               "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                               "#64FFFF", "#F5FFFF"))

# ......................................................................... ####
# 06/09/2020 ####
## Testing GenSA for numerical model ####
df_ref <- get_df_ref()

library(GenSA)
out_optim <- GenSA(
  fn         = optimise_this_gs_vcmax_jmax,
  par        = c( df_ref$vcmax_start,       df_ref$gs_start,       df_ref$jmax_start ), # starting values
  lower      = c( df_ref$vcmax_start*0.001, df_ref$gs_start*0.001, df_ref$jmax_start*0.001 ),
  upper      = c( df_ref$vcmax_start*1000,  df_ref$gs_start*1000,  df_ref$jmax_start*1000 ),
  args       = c(df_ref$tc, df_ref$patm, df_ref$co2, df_ref$vpd),
  iabs       = (df_ref$ppfd * df_ref$fapar),
  kphio      = df_ref$kphio,
  beta       = 146,
  method_jmaxlim_inst = "smith37",
  method_eb  = "plantecophys",
  maximize   = T
)

out_optim$par



# ......................................................................... ####
# 07/09/2020 ####
## KOEPPEN GEIGER KEY ####
c("Af"="#960000", "Am"="#FF0000", "As"="#FF6E6E", "Aw"="#FFCCCC", "BSh"="#CC8D14", "BSk"="#CCAA54", "BWh"="#FFCC00", "BWk"="#FFFF64", "Cfa"="#007800", "Cfb"="#005000", "Cfc"="#003200", "Csa"="#96FF00", "Csb"="#00D700", "Csc"="#00AA00", "Cwa"="#BEBE00", "Cwb"="#8C8C00", "Cwc"="#5A5A00", "Dfa"="#550055", "Dfb"="#820082", "Dfc"="#C800C8", "Dfd"="#FF6EFF", "Dsa"="#646464", "Dsb"="#8C8C8C", "Dsc"="#BEBEBE", "Dsd"="#E6E6E6", "Dwa"="#6E28B4", "Dwb"="#B464FA", "Dwc"="#C89BFA", "Dwd"="#C8C8FF", "EF"="#6496FF", "ET"="#6496FF", "Ocean"="#6496FF")

# ......................................................................... ####
# 10/09/2020 ####
## A_gross curves ####
df_ref <- get_df_ref()
# df_ref$ppfd <- 1800e-6
df_ref$nsteps <- df_ref$nsteps
df_tib <- tibble(
  tc_air = seq(0, 40, length.out = df_ref$nsteps),
  ac     = rep(NA, df_ref$nsteps),
  aj_farq = rep(NA, df_ref$nsteps),
  aj_smith = rep(NA, df_ref$nsteps),
  a_lim  = rep(NA, df_ref$nsteps),
  lim    = rep(NA, df_ref$nsteps),
  rd     = rep(NA, df_ref$nsteps)
)

for (i in 1:nrow(df_tib)) {
  ## Get temperature dependent variables:
  gammastar <- calc_gammastar(df_tib$tc_air[i], df_ref$patm)
  kmm <- calc_kmm(df_tib$tc_air[i], df_ref$patm)
  vcmax <- df_ref$vcmax_start * calc_ftemp_inst_vcmax(df_tib$tc_air[i], tcgrowth = df_ref$tc_growth)
  jmax <- df_ref$jmax_start * calc_ftemp_inst_jmax(df_tib$tc_air[i], tcgrowth = df_ref$tc_growth, tchome = df_ref$tc_home)
  
  ## Rd
  df_tib$rd[i]       <- calc_rd(df_tib$tc_air[i], vcmax = vcmax)
  
  ## Agross:
  df_tib$ac[i]       <- calc_ac(ci = df_ref$ci/10, gammastar = gammastar, kmm = kmm, vcmax = vcmax)$ac
  df_tib$aj_farq[i]  <- calc_aj(j_method = "farquhar89", kphio = df_ref$kphio*4, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
  df_tib$aj_smith[i] <- calc_aj(j_method = "smith37", kphio = df_ref$kphio, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
  df_tib$a_lim[i]    <- min(df_tib$ac[i], df_tib$aj_farq[i], df_tib$aj_smith[i])

  # ## Anet:
  # df_tib$ac[i]       <- df_tib$ac[i]       - df_tib$rd[i]
  # df_tib$aj_farq[i]  <- df_tib$aj_farq[i]  - df_tib$rd[i]
  # df_tib$aj_smith[i] <- df_tib$aj_smith[i] - df_tib$rd[i]
  # df_tib$a_lim[i]    <- min(df_tib$ac[i], df_tib$aj_farq[i], df_tib$aj_smith[i])

  ## Define limiting rate
  if (df_tib$a_lim[i] == df_tib$ac[i]) df_tib$lim[i]        <- "ac"
  if (df_tib$a_lim[i] == df_tib$aj_farq[i]) df_tib$lim[i]   <- "aj_farq"
  if (df_tib$a_lim[i] == df_tib$aj_smith[i]) df_tib$lim[i]  <- "aj_smith"
}


## .................................................................................................
## Gross
maximum <- max(df_tib$aj_smith)

df_plot <- df_tib %>%
  pivot_longer(cols = c(ac, aj_farq, aj_smith), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         a     = a/maximum)

(p_gross <- ggplot() +
  geom_line(data = df_plot %>% dplyr::filter(rate == "ac"),       aes(tc_air, a, color = "Ac"),    linetype = 1,   alpha = 0.5, size = 1.25)  +
  geom_line(data = df_plot %>% dplyr::filter(rate == "aj_smith"), aes(tc_air, a, color = "Aj (Smith)"),    linetype = 3,   alpha = 0.5, size = 1.25)  +
  geom_line(data = df_plot %>% dplyr::filter(rate == "aj_farq"),  aes(tc_air, a, color = "Aj (Farquhar)"),    linetype = 1,   alpha = 0.5, size = 1.25)  +
  geom_line(data = df_plot ,                                      aes(tc_air, a_lim, color = "min(Ac, Aj)"), linetype = 2,   alpha = 1, size = 1.25) +
  geom_line(data = df_plot, aes(tc_air, rd, color = "Rd")) +
  scale_color_manual(values = c("Ac" = "#FFCD08",
                                "Aj (Farquhar)" = "#DE5207",
                                "Aj (Smith)" = "#DE5207",
                                "min(Ac, Aj)" = "#F58706",
                                "Rd" = "black"),
                     name = "Rate:",
                     guide = guide_legend(override.aes=list(alpha=1, linetype = c(1,2,3,2,1)))) +
  ylim(0, 1.01) +
  xlab("Temperature [°C]") +
  ylab(bquote("Normalized rate [-]")) +
  ggtitle("Gross assimilation and respiration rates"))

# ggsave("~/projects/mscthesis/docs/fig-comparison-gross-rates.pdf", p1, height = 4, width = 5.5)
## .................................................................................................
## Anet
df_plot2 <- df_tib %>%
  pivot_longer(cols = c(ac, aj_farq, aj_smith), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim - rd,
         a     = a - rd,
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         a     = a/maximum)

(p_net<- ggplot() +
    geom_line(data = df_plot2 %>% dplyr::filter(rate == "ac"),       aes(tc_air, a, color = "Ac"),    linetype = 1,   alpha = 0.5, size = 1.25)  +
    geom_line(data = df_plot2 %>% dplyr::filter(rate == "aj_smith"), aes(tc_air, a, color = "Aj (Smith)"),    linetype = 3,   alpha = 0.5, size = 1.25)  +
    geom_line(data = df_plot2 %>% dplyr::filter(rate == "aj_farq"),  aes(tc_air, a, color = "Aj (Farquhar)"),    linetype = 1,   alpha = 0.5, size = 1.25)  +
    geom_line(data = df_plot2 ,                                      aes(tc_air, a_lim, color = "min(Ac, Aj)"), linetype = 2,   alpha = 1, size = 1.25) +
    scale_color_manual(values = c("Ac" = "#FFCD08",
                                  "Aj (Farquhar)" = "#DE5207",
                                  "Aj (Smith)" = "#DE5207",
                                  "min(Ac, Aj)" = "#F58706"),
                       name = "Rate:",
                       guide = guide_legend(override.aes=list(alpha=1, linetype = c(1,1,3,1)))) +
    ylim(0, 1.01) +
    xlab("Temperature [°C]") +
    ylab(bquote("Normalized rate [-]")) +
    ggtitle("Net assimilation rates"))

# ggsave("~/projects/mscthesis/docs/fig-comparison-net-rates.pdf", p2, height = 4, width = 5.5)

## .................................................................................................
## Comparison plot:
pnet_2 <- p_net + ylab("") + scale_color_manual(values = c("Ac" = "#FFCD08",
                                                        "Aj (Farquhar)" = "#DE5207",
                                                        "Aj (Smith)" = "#DE5207",
                                                        "min(Ac, Aj)" = "#F58706",
                                                        "Rd" = "black"),
                                             name = "Rate:",
                                             guide = guide_legend(override.aes=list(alpha=1, linetype = c(1,1,3,2,1)))) 

p_gross_2 <- p_gross + guides(color = "none")

(p3 <- p_gross_2 + pnet_2 + plot_layout(guides = "collect") & theme(legend.position = "bottom"))

# ggsave("~/projects/mscthesis/docs/fig-comparison-both-rates.pdf", p3, height = 4, width = 8)

# ......................................................................... ####
# 10/09/2020 ####
## Aj Sensitivities ####

### Setup ####
df_ref <- get_df_ref()
v_ppfd <- c(250, 500, 1000, 1500, 2000)*10^-6
v_jmax25 <- c(25, 50, 75, 100, 150, 200)*10^-6
v_vcmax25 <- c(25, 50, 75, 100, 150, 200)*10^-6

df_tib <- tibble(
  tc_air = seq(0, 40, length.out = df_ref$nsteps),
  ppfd   = rep(NA, df_ref$nsteps),
  jmax25   = rep(NA, df_ref$nsteps),
  ac     = rep(NA, df_ref$nsteps),
  aj_farq = rep(NA, df_ref$nsteps),
  aj_smith = rep(NA, df_ref$nsteps),
  a_lim  = rep(NA, df_ref$nsteps),
  lim    = rep(NA, df_ref$nsteps),
  rd     = rep(NA, df_ref$nsteps)
)

## Labeller:
vnames <-list("ac" = "Ac",
              "aj_smith"  = "Aj (Smith)",
              "aj_farq"  = "Aj (Farquhar)")



## .................................................................................................
## INSTANT PPFD ####
df_out_ppfd <- tibble()

for (p in v_ppfd) {
    df_loop <- df_tib
    df_loop$ppfd <- p
    
  for (i in 1:nrow(df_tib)) {
    ## Get temperature dependent variables:
    gammastar <- calc_gammastar(df_loop$tc_air[i], df_ref$patm)
    kmm <- calc_kmm(df_loop$tc_air[i], df_ref$patm)
    vcmax <- df_ref$vcmax_start * calc_ftemp_inst_vcmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth)
    jmax <- df_ref$jmax_start * calc_ftemp_inst_jmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth, tchome = df_ref$tc_home)
    
    ## Rd
    df_loop$rd[i]       <- calc_rd(df_loop$tc_air[i], vcmax = df_ref$vcmax_start)
    
    ## Agross:
    df_loop$ac[i]       <- calc_ac(ci = df_ref$ci/10, gammastar = gammastar, kmm = kmm, vcmax = vcmax)$ac
    df_loop$aj_farq[i]  <- calc_aj(j_method = "farquhar89", kphio = df_ref$kphio*4, jmax = jmax, ppfd = df_loop$ppfd[i], ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
    df_loop$aj_smith[i] <- calc_aj(j_method = "smith37", kphio = df_ref$kphio, jmax = jmax, ppfd = df_loop$ppfd[i], ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
    df_loop$a_lim[i]    <- min(df_loop$ac[i], df_loop$aj_farq[i], df_loop$aj_smith[i])
    
    # ## Anet:
    # df_loop$ac[i]       <- df_loop$ac[i]       - df_loop$rd[i]
    # df_loop$aj_farq[i]  <- df_loop$aj_farq[i]  - df_loop$rd[i]
    # df_loop$aj_smith[i] <- df_loop$aj_smith[i] - df_loop$rd[i]
    # df_loop$a_lim[i]    <- min(df_loop$ac[i], df_loop$aj_farq[i], df_loop$aj_smith[i])
    
    ## Define limiting rate
    if (df_loop$a_lim[i] == df_loop$ac[i]) df_loop$lim[i]        <- "ac"
    if (df_loop$a_lim[i] == df_loop$aj_farq[i]) df_loop$lim[i]   <- "aj_farq"
    if (df_loop$a_lim[i] == df_loop$aj_smith[i]) df_loop$lim[i]  <- "aj_smith"
  }
  df_out_ppfd <- bind_rows(list(df_out_ppfd, df_loop))
}

## .................................................................................................
## PPFD:
maximum <- max(df_out_ppfd$aj_smith, df_out_ppfd$aj_farq, df_out_ppfd$ac)

df_plot_ppfd <- df_out_ppfd %>%
  pivot_longer(cols = c(ac, aj_farq, aj_smith), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         a     = a/maximum,
         ppfd  = ifelse(rate == "ac", 200e-6, ppfd),
         ppfd  = as.factor(ppfd*10^6))

(p_ppfd <- df_plot_ppfd %>% ggplot() +
    aes(tc_air, a, color = ppfd, group = ppfd) +
    geom_line(alpha = 1, size = 1.5) +
    facet_wrap(~rate, labeller = vlabeller) +
    scale_color_viridis_d(name = bquote("PPFD [µmol" ~ m^-2 ~ s ^-1 ~ "]:  ")) +
    guides(color = guide_legend(ncol = 8, direction = "horizontal")) +
    ylab(bquote("Normalized rate [-]")) +
    xlab("Temperature [°C]") +
    theme(legend.position = "bottom") +
    ggtitle(bquote("Sensitivity of assimilation rates to light")))

# ggsave("~/projects/mscthesis/docs/fig-", p, height = 7, width = 8)

## .................................................................................................
## JMAX25 ####
df_out_jmax <- tibble()

for (j in v_jmax25) {
  df_loop <- df_tib
  df_loop$jmax25 <- j
  
  for (i in 1:nrow(df_tib)) {
    ## Get temperature dependent variables:
    gammastar <- calc_gammastar(df_loop$tc_air[i], df_ref$patm)
    kmm <- calc_kmm(df_loop$tc_air[i], df_ref$patm)
    vcmax <- df_ref$vcmax_start * calc_ftemp_inst_vcmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth)
    jmax  <- df_loop$jmax25[i] * calc_ftemp_inst_jmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth, tchome = df_ref$tc_home)
    
    ## Rd
    df_loop$rd[i]       <- calc_rd(df_loop$tc_air[i], vcmax = df_ref$vcmax_start)
    
    ## Agross:
    df_loop$ac[i]       <- calc_ac(ci = df_ref$ci/10, gammastar = gammastar, kmm = kmm, vcmax = vcmax)$ac
    df_loop$aj_farq[i]  <- calc_aj(j_method = "farquhar89", kphio = df_ref$kphio*4, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
    df_loop$aj_smith[i] <- calc_aj(j_method = "smith37", kphio = df_ref$kphio, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
    df_loop$a_lim[i]    <- min(df_loop$ac[i], df_loop$aj_farq[i], df_loop$aj_smith[i])
    
    # ## Anet:
    df_loop$ac[i]       <- df_loop$ac[i]       - df_loop$rd[i]
    df_loop$aj_farq[i]  <- df_loop$aj_farq[i]  - df_loop$rd[i]
    df_loop$aj_smith[i] <- df_loop$aj_smith[i] - df_loop$rd[i]
    df_loop$a_lim[i]    <- min(df_loop$ac[i], df_loop$aj_farq[i], df_loop$aj_smith[i])
    
    ## Define limiting rate
    if (df_loop$a_lim[i] == df_loop$ac[i]) df_loop$lim[i]        <- "ac"
    if (df_loop$a_lim[i] == df_loop$aj_farq[i]) df_loop$lim[i]   <- "aj_farq"
    if (df_loop$a_lim[i] == df_loop$aj_smith[i]) df_loop$lim[i]  <- "aj_smith"
  }
  df_out_jmax <- bind_rows(list(df_out_jmax, df_loop))
}

## .................................................................................................
## Plot
maximum <- max(df_out_jmax$aj_smith, df_out_jmax$aj_farq, df_out_jmax$ac)

df_plot_jmax <- df_out_jmax %>%
  pivot_longer(cols = c(ac, aj_farq, aj_smith), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         a     = a/maximum,
         jmax25  = ifelse(rate == "ac", v_jmax25[1], jmax25),
         jmax25  = as.factor(jmax25*10^6))

(p_jmax <- df_plot_jmax %>% ggplot() +
    aes(tc_air, a, color = jmax25, group = jmax25) +
    geom_line(alpha = 1, size = 1.5) +
    facet_wrap(~rate, labeller = vlabeller) +
    scale_color_viridis_d(name = bquote(J[max]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]:  "),
                          option = "magma") +
    guides(color = guide_legend(ncol = 8, direction = "horizontal")) +
    ylab(bquote("Normalized rate [-]")) +
    xlab("Temperature [°C]") +
    theme(legend.position = "bottom") +
    ggtitle(bquote("Sensitivity of assimilation rates to " * J[max]^25)))

# ggsave("~/projects/mscthesis/docs/fig-", p, height = 7, width = 8)

## .................................................................................................
## Combine plots JMAX25 and INSTANT PPFD ####
p <- p_ppfd / p_jmax

# ggsave("~/projects/mscthesis/docs/fig-sens-aj-ppfd-jmax.pdf", p, height = 7, width = 8)

## .................................................................................................
## Ac sensitivities ####
## .................................................................................................
## VMAX:
df_out_vcmax <- tibble()

for (v in v_vcmax25) {
  df_loop <- df_tib
  df_loop$vcmax25 <- v
  
  for (i in 1:nrow(df_tib)) {
    ## Get temperature dependent variables:
    gammastar <- calc_gammastar(df_loop$tc_air[i], df_ref$patm)
    kmm <- calc_kmm(df_loop$tc_air[i], df_ref$patm)
    vcmax <- df_loop$vcmax25[i] * calc_ftemp_inst_vcmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth)
    jmax  <- df_ref$jmax_start * calc_ftemp_inst_jmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth, tchome = df_ref$tc_home)
    
    ## Rd
    df_loop$rd[i]       <- calc_rd(df_loop$tc_air[i], vcmax = df_loop$vcmax25[i])
    
    ## Agross:
    df_loop$ac[i]       <- calc_ac(ci = df_ref$ci/10, gammastar = gammastar, kmm = kmm, vcmax = vcmax)$ac
    df_loop$aj_farq[i]  <- calc_aj(j_method = "farquhar89", kphio = df_ref$kphio*4, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
    df_loop$aj_smith[i] <- calc_aj(j_method = "smith37", kphio = df_ref$kphio, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar)$aj
    df_loop$a_lim[i]    <- min(df_loop$ac[i], df_loop$aj_farq[i], df_loop$aj_smith[i])
    
    # ## Anet:
    # df_loop$ac[i]       <- df_loop$ac[i]       - df_loop$rd[i]
    # df_loop$aj_farq[i]  <- df_loop$aj_farq[i]  - df_loop$rd[i]
    # df_loop$aj_smith[i] <- df_loop$aj_smith[i] - df_loop$rd[i]
    # df_loop$a_lim[i]    <- min(df_loop$ac[i], df_loop$aj_farq[i], df_loop$aj_smith[i])
    
    ## Define limiting rate
    if (df_loop$a_lim[i] == df_loop$ac[i]) df_loop$lim[i]        <- "ac"
    if (df_loop$a_lim[i] == df_loop$aj_farq[i]) df_loop$lim[i]   <- "aj_farq"
    if (df_loop$a_lim[i] == df_loop$aj_smith[i]) df_loop$lim[i]  <- "aj_smith"
  }
  df_out_vcmax <- bind_rows(list(df_out_vcmax, df_loop))
}

maximum <- max(df_out_vcmax$aj_smith, df_out_vcmax$aj_farq, df_out_vcmax$ac)

df_plot_vcmax <- df_out_vcmax %>%
  pivot_longer(cols = c(ac, aj_farq, aj_smith), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         a     = a/maximum,
         vcmax25  = ifelse(str_detect(rate, "aj"), v_vcmax25[1], vcmax25),
         vcmax25  = as.factor(vcmax25*10^6))

(p_vcmax <- df_plot_vcmax %>% ggplot() +
    aes(tc_air, a, color = vcmax25, group = vcmax25) +
    geom_line(alpha = 1, size = 1.5) +
    facet_wrap(~rate, labeller = vlabeller) +
    scale_color_viridis_d(name = bquote(V[cmax]^25 ~ "[µmol" ~ m^-2 ~ s ^-1 ~ "]:  "),
                          option = "magma") +
    guides(color = guide_legend(ncol = 8, direction = "horizontal")) +
    ylab(bquote("Normalized rate [-]")) +
    xlab("Temperature [°C]") +
    theme(legend.position = "bottom") +
    ggtitle(bquote("Sensitivity of assimilation rates to " * V[cmax]^25)))

# ......................................................................... ####
# 13/09/2020 ####
## Rogers et al. 2019 - theta data ####
temp <- c(5.0, 4.8 , 4.9 , 5.0 , 5.0 , 4.9 , 14.9, 14.9, 14.9, 14.9, 14.9, 14.9, 24.9, NA, 24.8, 24.8, 24.9, NA)
theta <- c(0.42, 0.51 , 0.27 , 0.34 , 0.37 , 0.59 , 0.43 , 0.51 , 0.38 , 0.52 , 0.60 , 0.56 , 0.67 , NA, 0.56 , 0.72, 0.71 , NA)
lm(theta ~ temp) %>% summary()

## .................................................................................................
## rpmodel sensitivity ####
### Setup ####
df_ref <- get_df_ref()
# v_ppfd <- c(250, 500, 1000, 1500, 2000)*10^-6
v_ppfd <- c(250, 500, 750,  1000, 1250, 1500)*10^-6

df_tib <- tibble(
  tc_air = seq(0, 40, length.out = df_ref$nsteps),
  ppfd   = rep(NA, df_ref$nsteps),
  limitation   = rep(NA, df_ref$nsteps),
  ac     = rep(NA, df_ref$nsteps),
  aj = rep(NA, df_ref$nsteps),
  a_lim  = rep(NA, df_ref$nsteps),
  lim    = rep(NA, df_ref$nsteps),
  rd     = rep(NA, df_ref$nsteps),
  vcmax25 = rep(NA, df_ref$nsteps),
  jmax25 = rep(NA, df_ref$nsteps)
)

### Loop ####
## GROWTH PPFD:
# df_ref$ppfd <- 1500 * 10^-6
# df_ref$tc_home <- 15
# df_ref$tc_growth <- 15

df_out_ppfd <- tibble()

for (limitation in c("smith37", "farquhar89")) {
  df_loop <- df_tib
  df_loop$limitation <- limitation
  
  for (p in v_ppfd) {
    
    ## Get acclimated variables:
    df_loop$ppfd <- p
    
    out <- rpmodel(ppfd = p, # Looped variable
                   method_jmaxlim = limitation, 
                   tc_growth_air = df_ref$tc_growth, tc_home = df_ref$tc_home, vpd = df_ref$vpd, co2 = df_ref$co2, fapar = 1, patm = df_ref$patm, # Forcings
                   kphio = df_ref$kphio, beta = 146.0, soilm = stopifnot(!do_soilmstress), meanalpha = 1.0, apar_soilm = 0.0, bpar_soilm = 0.73300, c4 = FALSE, soilmstress = 1, # Local parameters
                   method_optim   = "analytical", method_optci   = "prentice14", method_ftemp   = "kumarathunge19", method_eb      = "off", # Calculation methods
                   do_ftemp_kphio = TRUE, do_ftemp_theta = F, do_soilmstress = FALSE, returnvar = NULL, verbose = FALSE # Other settings
    )
    
    df_loop$jmax25  <- out$jmax25
    df_loop$vcmax25 <- out$vcmax25
    
    for (i in 1:nrow(df_tib)) {
      ## Get instantaneous reaction:
      
      ## Get temperature dependent variables:
      gammastar <- calc_gammastar(df_loop$tc_air[i], df_ref$patm)
      kmm <- calc_kmm(df_loop$tc_air[i], df_ref$patm)
      vcmax <- df_loop$vcmax25[i] * calc_ftemp_inst_vcmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth)
      jmax <-  df_loop$jmax25[i] * calc_ftemp_inst_jmax(df_loop$tc_air[i], tcgrowth = df_ref$tc_growth, tchome = df_ref$tc_home)
      
      ## Rd
      df_loop$rd[i] <- calc_rd(df_loop$tc_air[i], vcmax = df_loop$vcmax25[i])
      
      ## Agross:
      df_loop$ac[i] <- calc_ac(ci = df_ref$ci/10, gammastar = gammastar, kmm = kmm, vcmax = vcmax)$ac
      
      if (df_loop$limitation[i] == "farquhar89") {
        kphio <- df_ref$kphio * 4
      } else {
        kphio <- df_ref$kphio
      }
      df_loop$aj[i]    <- calc_aj(j_method = df_loop$limitation[i], kphio = kphio, jmax = jmax, ppfd = df_ref$ppfd, ci = df_ref$ci/10, fapar = 1, gammastar = gammastar, theta = 0.85)$aj
      df_loop$a_lim[i] <- min(df_loop$ac[i], df_loop$aj[i])
      
      # ## Anet:
      df_loop$ac[i]  <- df_loop$ac[i]  - df_loop$rd[i]
      df_loop$aj[i]  <- df_loop$aj[i]  - df_loop$rd[i]
      df_loop$a_lim[i] <- min(df_loop$ac[i], df_loop$aj[i])
      
      ## Define limiting rate
      if (df_loop$a_lim[i] == df_loop$ac[i]) df_loop$lim[i]   <- "ac"
      if (df_loop$a_lim[i] == df_loop$aj[i]) df_loop$lim[i]   <- "aj"
    }
    
    df_out_ppfd <- bind_rows(list(df_out_ppfd, df_loop))
  }
}

## .................................................................................................
### PLOTS ####
vnames <-list("ac" = "Ac", "aj"  = "Aj")

#### FARQUHAR
df_farq <- df_out_ppfd %>% dplyr::filter(limitation == "farquhar89")

maximum <- max(df_farq$aj, df_farq$ac)

df_plot_farq <- df_farq %>%
  pivot_longer(cols = c(ac, aj), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         # a     = a/maximum,
         a     = a*10^6,
         ppfd  = as.factor(ppfd*10^6))


p_farq <- df_plot_farq %>% ggplot() +
    aes(tc_air, a, color = ppfd, group = ppfd) +
    geom_line(alpha = 1, size = 1) +
    facet_wrap(~rate, labeller = vlabeller) +
    scale_color_viridis_d(name = bquote("Growth PPFD [µmol" ~ m^-2 ~ s ^-1 ~ "]:  ")) +
    guides(color = guide_legend(ncol = 4, direction = "horizontal")) +
    ylab(bquote("Rate [µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
    xlab("Temperature [°C]") +
    theme(legend.position = "bottom") +
    ggtitle(bquote("Farquhar - Formulation")) +
    ylim(0,35)


#### SMITH
df_smith <- df_out_ppfd %>% dplyr::filter(limitation == "smith37")

maximum <- max(df_smith$aj, df_smith$ac)

df_plot_smith <- df_smith %>%
  pivot_longer(cols = c(ac, aj), values_to = "a", names_to = "rate") %>% 
  mutate(rate = as.factor(rate),
         a_lim = a_lim/maximum,
         rd    = rd/maximum,
         # a     = a/maximum,
         a     = a*10^6,
         ppfd  = as.factor(ppfd*10^6))


p_smith <- df_plot_smith %>% ggplot() +
    aes(tc_air, a, color = ppfd, group = ppfd) +
    geom_line(alpha = 1, size = 1) +
    facet_wrap(~rate, labeller = vlabeller) +
    scale_color_viridis_d(name = bquote("Growth PPFD [µmol" ~ m^-2 ~ s ^-1 ~ "]:  ")) +
    guides(color = guide_legend(ncol = 4, direction = "horizontal")) +
    ylab(bquote("Rate [µmol" ~ m^-2 ~ s ^-1 ~ "]")) +
    xlab("Temperature [°C]") +
    xlab("") +
    theme(legend.position = "bottom") +
    ggtitle(bquote("Smith-Formulation")) +
    ylim(0,35)


(p_both <- p_smith / p_farq +
    plot_layout(guides = "collect") &
    plot_annotation(title = "Sensitivity of assimilation rates to growth light conditions") &
    theme(legend.position = "bottom"))

# ggsave("~/projects/mscthesis/docs/fig-sens-acaj-growth-ppfd.pdf", p_both, height = 7, width = 5.5)









#___________________________________________________________________________________________________ ####
# . ####
# . ####
# . ####
# . ####
# . ####
# . ####
# . ####
# . ####
# . ####.
