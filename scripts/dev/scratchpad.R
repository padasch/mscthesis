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
vcmax_fix <- 2
jmax_fix  <- 8
gs_fix    <- 0.6

## Vectors of starting values
vcmax_v <- c(1, 2, 4, 6, 8, 10, 12, 15)
jmax_v  <- c(1, 2, 4, 6, 8, 10, 12, 15)
gs_v    <- c(1, 2, 4, 6, 8, 10, 12, 15)/10

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
  calc_optimal_tcleaf_vcmax_jmax <- function(tc_leaf = 25,
                                             patm = 101325,
                                             co2 = 400,
                                             vpd = 1000,
                                             ppfd = 130,
                                             fapar = 1,
                                             kphio = 0.05,
                                             beta = 146,
                                             c_cost = 0.41,
                                             vcmax_start = df_start$vcmax_start[i],
                                             gs_start = df_start$gs_start[i],
                                             jmax_start = df_start$jmax_start[i],
                                             method_jmaxlim_inst = "smith37") {
    
    out_optim <- optimr::optimr(
      par        = c( vcmax_start,       gs_start,       jmax_start ), # starting values
      lower      = c( vcmax_start*0.001, gs_start*0.001, jmax_start*0.001 ),
      upper      = c( vcmax_start*1000,  gs_start*1000,  jmax_start*2 ),
      fn         = optimise_this_tcleaf_vcmax_jmax,
      args       = c(tc_leaf, patm, co2, vpd),
      iabs       = (ppfd * fapar),
      kphio      = kphio,
      beta       = beta,
      c_cost     = c_cost/4,
      method_jmaxlim_inst = method_jmaxlim_inst,
      method     = "L-BFGS-B",
      maximize   = TRUE,
      control    = list(maxit=1000)
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
    
    out <- list(out_optim = out_optim,
                varlist   = varlist)
    
    return(out)
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
  p_num_s37 <- p_num_s37 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + ggtitle("Smith-Formulation")
  p_num_f89 <- p_num_f89 + ylim(0, 250) + xlim(0, 250) + ylab("Observed Values [µmol/m2/s]") + xlab("Predicted Values [µmol/m2/s]") + ggtitle("Farquhar-Formulation")
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


## P21 Global Map
# Get df_evaluation_p21 from analysis of acclimate pmodel

siteinfo <- df_evaluation_p21 %>% mutate(lon = purrr::map_dbl(data_raw, . %>% pull(lon) %>% unique()),
                                         lat = purrr::map_dbl(data_raw, . %>% pull(lat) %>% unique()))

kg <- readOGR(dsn = "~/data/climate_zones/1976-2000_GIS", layer = "1976-2000")
kg <- spTransform(kg, CRS("+proj=longlat +datum=WGS84"))
kg_f <- fortify(kg, region = "GRIDCODE")
key <- data.frame( id = c(11, 12, 13, 14, 21, 22, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 61, 62),
                   cl = c( "Af", "Am", "As", "Aw", "BWk", "BWh", "BSk", "BSh", "Cfa", "Cfb", "Cfc", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Dfa", "Dfb", "Dfc", "Dfd", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "EF", "ET"))
kg_final <- merge(kg_f, key)

(p <- siteinfo %>%
  ggplot() +
  # borders("world", colour="black", fill="gray50") +
  geom_polygon(data = kg_final, aes(x = long, y = lat, group = group, fill = cl), alpha = 0.8) +
  geom_point(aes(x = lon, y = lat), color = "black", size = 3, pch = 4) +
  ylab("Latitude (Decimal Degree)") +
  xlab("Longitude (Decimal Degree)") +
  coord_cartesian(ylim=c(-80, 80)) +
  scale_fill_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                               "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                               "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                               "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                               "#64FFFF", "#F5FFFF")) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 12, direction = "horizontal", title.position = "top")) +
  labs(fill = "Köppen-Geiger Climate Zone") +
  ggtitle(bquote(("Global distribution of " ~ V[cmax] ~ "observations"))))

ggsave("~/projects/mscthesis/docs/fig-global-map-p21.pdf", p, height = 7, width = 8)


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

kg <- readOGR(dsn = "~/data/climate_zones/1976-2000_GIS", layer = "1976-2000")
kg <- spTransform(kg, CRS("+proj=longlat +datum=WGS84"))
kg_f <- fortify(kg, region = "GRIDCODE")
key <- data.frame( id = c(11, 12, 13, 14, 21, 22, 26, 27, 31, 32, 33, 34, 35, 36, 37, 38, 39, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 61, 62),
                   cl = c( "Af", "Am", "As", "Aw", "BWk", "BWh", "BSk", "BSh", "Cfa", "Cfb", "Cfc", "Csa", "Csb", "Csc", "Cwa", "Cwb", "Cwc", "Dfa", "Dfb", "Dfc", "Dfd", "Dsa", "Dsb", "Dsc", "Dsd", "Dwa", "Dwb", "Dwc", "Dwd", "EF", "ET"))
kg_final <- merge(kg_f, key)


(p <- siteinfo %>%
    ggplot() +
    # borders("world", colour="black", fill="gray50") +
    geom_polygon(data = kg_final, aes(x = long, y = lat, group = group, fill = cl), alpha = 0.8) +
    geom_point(aes(x = lon, y = lat), color = "black", size = 3, pch = 4) +
    ylab("Latitude (Decimal Degree)") +
    xlab("Longitude (Decimal Degree)") +
    coord_cartesian(ylim=c(-80, 80)) +
    scale_fill_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                                 "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                                 "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                                 "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                                 "#64FFFF", "#F5FFFF")) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 12, direction = "horizontal", title.position = "top")) +
    labs(fill = "Köppen-Geiger Climate Zone") +
    ggtitle(bquote(("Global distribution of " ~ V[cmax] ~ "observations"))))

ggsave("~/projects/mscthesis/docs/fig-global-map-k19.pdf", p, height = 7, width = 8)


siteinfo %>%
  ggplot() +
  # borders("world", colour="black", fill="gray50") +
  # geom_line(data = kg_2, aes(x = Lon, y = Lat, color = Cls)) +
  # geom_polygon(kg_f, aes(x = long, y = lat, ))
  geom_polygon(data = kg_final, aes(x = long, y = lat, group = group, fill = cl, alpha = 0.9)) +
  geom_point(aes(x = lon, y = lat), color = "black", size = 3, shape = 18) +
  geom_label_repel(aes(x = lon, y = lat, label = source),
                   box.padding = 0.5,
                   point.padding = 0.5,
                   min.segment.length = 0, # draw all line segments
                   arrow = arrow(length = unit(0.015, "npc"))) +  # draw arrows
  ylab("Latitude (Decimal Degree)") +
  xlab("Longitude (Decimal Degree)") +
  # coord_cartesian(ylim=c(-85, 80)) +
  scale_fill_manual(values = c("#960000", "#FF0000", "#FF6E6E", "#FFCCCC",
                               "#CC8D14", "#CCAA54", "#FFCC00", "#FFFF64",
                               "#007800", "#005000", "#003200", "#96FF00", "#00D700", "#00AA00", "#BEBE00", "#8C8C00", "#5A5A00",
                               "#550055", "#820082", "#C800C8", "#FF6EFF", "#646464", "#8C8C8C", "#BEBEBE", "#E6E6E6", "#6E28B4", "#B464FA", "#C89BFA", "#C8C8FF", "#6496FF",
                               "#64FFFF", "#F5FFFF"),
                    name = "Koeppen-Geiger Climate Zone",
                    guide = guide_legend(
                      direction = "horizontal",
                      title.position = "top",
                      ncol = 12))  +
  theme(legend.position = "bottom")

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

## Removing duplicate entries in P21 data ####

















































