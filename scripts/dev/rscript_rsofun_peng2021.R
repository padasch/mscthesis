# PENG21 DATA LOAD
# PENG 2021 PREPARATIONS ####
## Date formatting ####
df_P21_clean <- read_csv("~/data/mscthesis/raw/leaf_traits_peng2021/combined_leaf_traits_updated.csv",
												 col_types = cols(Jmax = col_double(),
												 								 Jmax25 = col_double(), Date = col_character(),
												 								 Year = col_character(), carea = col_double(),
												 								 C_percent = col_double(), Author = col_character(),
												 								 author = col_character(), Aarea = col_double(),
												 								 rep = col_character())) %>%

	# To clean up date formatting: (instead of paste0(...), just use [.] to fin "." pattern)
	mutate(Date = ifelse(grepl("/", Date), as_date(Date) %>% paste0(), Date),
				 Date = ifelse(grepl(paste0("."), Date, fixed = T), dmy(Date) %>% paste0(), Date),
				 Date = ifelse(!(grepl(paste0("."), Date, fixed = T) & grepl("/", Date)), as_date(Date) %>% paste0(), Date),
				 Date = as_date(Date))

saveRDS(df_P21_clean, "~/data/mscthesis/raw/leaf_traits_peng2021/df_P21_clean.rds")
df_P21_clean <- readRDS("~/data/mscthesis/raw/leaf_traits_peng2021/df_P21_clean.rds")

## Siteinfo for cluster (distinct sites but not dates!) ####
df_P21_cluster <- df_P21_clean %>%
	drop_na(Date, Vcmax25) %>%
	dplyr::select(lon, lat, z, start_yr, end_yr, site) %>%
	distinct() %>%
	rowid_to_column("ID") %>% # Shortcut for non-distinct sitenames (small diff in lat, lon)
	rename(sitename   = ID,
				 year_start = start_yr,
				 year_end   = end_yr) %>%
	mutate(whc = 170)

saveRDS(df_P21_cluster, "~/data/mscthesis/raw/leaf_traits_peng2021/siteinfo_for_cluster_PENG21.rds")

# ...................................................................................................................
# RSOFUN LOOP FOR P21 DATA ####
all_sites <- readRDS("~/data/mscthesis/raw/leaf_traits_peng2021/siteinfo_for_cluster_PENG21.rds") %>%
	mutate(sitename = as.character(sitename))

## LOOP FOR OUTPUT ####

for (i in 1:nrow(all_sites)) {

	if( i == 44 | i == 70) {
		i <- i + 1
	}

	cat("\014")
	cat("Site-ID:", i, "/", nrow(all_sites), "\n")

	# Get siteinfo
	siteinfo <- all_sites[i, ] %>%
		# Add elevation data
		left_join(readRDS(paste0("~/data/mscthesis/raw/leaf_traits_peng2021/df_etopo_PENG21_", i, ".rds"))) %>%
		# Add tc_home data
		left_join(readRDS(paste0("~/data/mscthesis/raw/leaf_traits_peng2021/df_tc_home_PENG21_", i, ".rds")))

	# Get meteorological data
	ddf_meteo <- readRDS(paste0("~/data/mscthesis/raw/leaf_traits_peng2021/ddf_meteo_PENG21_", i, ".rds"))

	# Get uniformal 1.0 fapar dataframe
	ddf_fapar_unity <- ingest(
		siteinfo  = siteinfo,
		source    = "fapar_unity"
	)

	# CO2 data has to be loaded here instead via cluster
	df_co2 <- ingestr::ingest(
		siteinfo,
		source  = "co2_mlo",
		verbose = FALSE
	)

	## Simulation settings
	params_siml <- list(
		spinup             = TRUE,      # to bring soil moisture to steady state
		spinupyears        = 10,        # number of spinup years. 10 is enough for soil moisture.
		recycle            = 1,         # number of years recycled during spinup
		soilmstress        = FALSE,     # boolean for whether soil moisture stress function is included
		tempstress         = FALSE,     # boolean for whether temperature stress function is included
		calc_aet_fapar_vpd = FALSE,     # set to FALSE - should be dropped again
		in_ppfd            = TRUE,      # if available from forcing files, set to TRUE
		in_netrad          = FALSE,     # if available from forcing files, set to TRUE
		outdt              = 1,
		ltre               = FALSE,
		ltne               = FALSE,
		ltrd               = FALSE,
		ltnd               = FALSE,
		lgr3               = TRUE,
		lgn3               = FALSE,
		lgr4               = FALSE
	)

	# Has to be updated when changing rsofun/pmodel
	params_modl <- list(
		kphio           = 0.09423773,
		soilm_par_a     = 0.33349283,
		soilm_par_b     = 1.45602286
	)

	df_soiltexture <- bind_rows(
		top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
		bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1)
	)

	## Collect all drivers
	df_drivers <- collect_drivers_sofun(
		siteinfo       = siteinfo,
		params_siml    = params_siml,
		meteo          = ddf_meteo,
		fapar          = ddf_fapar_unity,
		co2            = df_co2,
		df_soiltexture = df_soiltexture
	)

	# Drop leap years cause of issues with rsofun
	df_drivers <- df_drivers %>%
		dplyr::select(sitename, forcing) %>%
		unnest(forcing) %>%
		dplyr::filter(!(month(date)==2 & mday(date)==29)) %>%
		group_by(sitename) %>%
		nest() %>%
		rename(forcing = data) %>%
		right_join(
			df_drivers %>%
				dplyr::select(-forcing),
			by = "sitename"
		) %>%
		ungroup()

	# saveRDS(df_drivers, "~/data/mscthesis/df_drivers.rds")
	# df_drivers <- readRDS("~/data/mscthesis/df_drivers.rds")

	## Run rsofun
	# detach("package:rsofun", unload = TRUE)
	# library(rsofun)

	df_output <- runread_pmodel_f(
		df_drivers,
		params_modl = params_modl,
		makecheck = TRUE,
		parallel = FALSE
	)

	if(TRUE){
		df_output <- df_output %>%

			# Attaching tc_leaf from df_drivers
			left_join(df_drivers %>%
									unnest(c("forcing", "siteinfo")))
	}

	saveRDS(df_output, paste0("~/data/mscthesis/raw/leaf_traits_peng2021/df_output_PENG21_", i, ".rds"))

}

## LOOP FOR BINDING ####
df_output <- readRDS("~/data/mscthesis/raw/leaf_traits_peng2021/df_output_PENG21_1.rds")
df_output_P21_all <- df_output %>% slice(0)

for (i in 1:nrow(all_sites)) {

	cat("\014")
	cat("Site-ID:", i, "/", nrow(all_sites), "\n")

	if( i == 44 | i == 70) {i <- i + 1}
	df_output_P21_all <- bind_rows(list(df_output_P21_all, readRDS(paste0("~/data/mscthesis/raw/leaf_traits_peng2021/df_output_PENG21_", i, ".rds"))))
}

# saveRDS(df_output_P21_all, "~/data/mscthesis/raw/leaf_traits_peng2021/df_rsofun_wang2017_tau30_P21.rds") # Adapt based on rsofun settings
# df_output_P21_all <- readRDS("~/data/mscthesis/raw/leaf_traits_peng2021/df_rsofun_wang2017_tau30_P21.rds")

df_drivers_p21 <- df_output_P21_all %>% 
    nest(forcing = c(date, ppfd, temp, patm, vpd, fapar, co2, rain, snow, prec, qair, vapr, ccov_int, ccov, doy, elv, z, myvar)) %>%
    nest(siteinfo = c(lon, lat, year_end, whc, tc_home, date_start, date_end)) %>% 
    mutate(sitename = site) %>% 
    rename(rsofun_out = data) %>% 
    dplyr::select(-site) 

saveRDS(df_drivers_p21, "~/data/mscthesis/final/df_drivers_p21.rds")


# RSOFUN DATES WITH MEASUREMENTS ####
df_output_P21_red <- left_join(
    df_P21_clean %>%
        drop_na(lon, lat, Date, Vcmax25) %>%
        dplyr::select(lon, lat, Date, site, Vcmax25) %>%
        rename(
            date        = Date,
            vcmax25_p21 = Vcmax25,
            sitename    = site
        ) %>%
        mutate(vcmax25_p21  = vcmax25_p21 * 10 ^ -6),
    # from umol to mol
    
    df_output_P21_all %>%
        dplyr::select(
            site,
            date,
            vcmax25,
            tc_growth,
            ppfd,
            tc_leaf,
            tc_home,
            elv,
            vpd,
            patm
        ) %>%
        rename(sitename = site,
               vcmax25_rso = vcmax25)
) %>% 
	mutate(sitename = as.factor(sitename)) %>%
	dplyr::filter(!(is.na(tc_growth))) # Dropping values with dates in the future

saveRDS(df_output_P21_red, "~/data/mscthesis/df_output_P21_red.rds") # Adapt based on rsofun settings
df_output_P21_red <- readRDS("~/data/mscthesis/df_output_P21_red.rds")


# CODE CHECKS ####
# Exchange "df_P21_clean" with "df_output_P21_red" to check if
# distinct sites and dates are kept correctly. For this, change "lon, late" into "sitename",
# "Date" into "agg_date" and "Vcmax25" into "vcmax25_p21".

# Number of total entries
df_P21_clean %>% nrow()

# Number of entries with distinct date, location and vcmax25
df_P21_clean %>% drop_na(lon, lat, Date, Vcmax25) %>% nrow()

# Number of entries with distinct sites and dates
df_P21_clean %>% drop_na(lon, lat, Date, Vcmax25) %>%
	dplyr::select(Date, site) %>% distinct() %>% nrow()

# Number of entries with distinct sites
df_P21_clean %>% drop_na(lon, lat, Date, Vcmax25) %>%
	dplyr::select(site) %>% distinct() %>% nrow()






































