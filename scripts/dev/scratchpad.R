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

df_x_unnest <- df_x %>% dplyr::select(-tc_growth_air, -tc_home) %>% unnest(forcing)

df_1  <- df_drivers$forcing[[1]]
dat_1 <- calc_optimal_tcleaf_vcmax_jmax(tc_leaf = df_1$temp,
                                        patm = df_1$patm,
                                        co2 = df_1$co2,
                                        vpd = df_1$vpd,
                                        ppfd = df_1$ppfd,
                                        fapar = df_1$fapar,
                                        kphio = df_1$kphio,
                                        beta = df_1$beta,
                                        c_cost = df_1$c_cost,
                                        vcmax_start = df_1$vcmax_start,
                                        gs_start = df_1$gs_start,
                                        jmax_start = df_1$jmax_start,
                                        method_jmaxlim_inst = "smith37")



























































