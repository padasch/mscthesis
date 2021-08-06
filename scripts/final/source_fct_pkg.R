# Scripts used by Pascal

## Setup
library(rsofun)
library(ingestr)
library(rpmodel)
library(lme4)
library(tidyverse)
library(patchwork)
library(ggridges)
library(plantecophys)
library(beepr)
library(optimr)
library(units)
library(tealeaves)
library(magrittr)
library(ggmap)
library(ggrepel)
library(spatialEco)
library(ggpmisc)

source("~/projects/mscthesis/scripts/final/rpmodel_mainroutines.R")
source("~/projects/mscthesis/scripts/final/rpmodel_subroutines.R")

if (F) {
    source('~/Polybox/2_ETH/ONGOING/msc-thesis/rproject/scripts/fct_inst_rpmodel.R')
    source('~/Polybox/2_ETH/ONGOING/msc-thesis/rproject/scripts/fct_non_linear_fit.R')
    source('~/Polybox/2_ETH/ONGOING/msc-thesis/rproject/scripts/fct_leaf_energy_balance.R')  
    
    path = "~/data/photom"
    #source(paste0(path,"/R/Functions/functions_for_analysis.R"))
    #source(paste0(path,"/R/Functions/metdataprocess.R"))
    #source(paste0(path,"/R/Functions/rsq_mm.R"))
    pkgs <- c("doBy", "magicaxis", "RColorBrewer", "propagate", "gplots", "readxl", "maps", "mapdata", "rgeos", "sp", "raster", "nlstools", "rgl", "mgcv", "scales", "data.table", "dplyr", "dismo", "AICcmodavg", "png", "plantecophys", "nlme", "lme4", "broom", "car", "lubridate", "tidyverse")
    lapply(pkgs, library, character.only = TRUE)
    
    # Call are functions for instantaneous pmodel:
    x <- getwd()
    setwd("~/Polybox/2_ETH/ONGOING/msc-thesis/rproject/scripts/000_numerical_pmodel")
    files.sources <- list.files()[-1]
    sapply(files.sources, source)
    setwd(x)
    remove(x)
    
    # K19 dataset key:
    # Siteinfo key
    site_data_key <- function(){
        tib <- tibble(sitename = c("alida",
                                   "anna",
                                   "arctic",
                                   "ellsworth",
                                   "han_cypress",
                                   "han_pine",
                                   "hikosaka",
                                   "kelly",
                                   "kelsey",
                                   "medlyn_eucalypt",
                                   "medlyn_pine",
                                   "onoda",
                                   "savanna",
                                   "slot",
                                   "tarvainen_pine",
                                   "tarvainen_spruce",
                                   "togashi",
                                   "tribuzi",
                                   "wang",
                                   "wtc"),
                      dataset = sitename,
                      dataset_yr = c("Mau (unpub.)",
                                     "Jensen (2015)",
                                     "Rogers (2017)",
                                     "Ellsworth (unpub.)",
                                     "Han (2006)",
                                     "Han (2004)",
                                     "Hikosaka (2007)",
                                     "Kelly (2014)",
                                     "Carter (unpub.)",
                                     "Medlyn (2007)",
                                     "Medlyn (2002)",
                                     "Onoda (2005)",
                                     "Cernusak (2011)",
                                     "Slot (2017)",
                                     "Tarvainen (2018)",
                                     "Tarvainen (2013)",
                                     "Togashi (2017)",
                                     "Tribuzy (2005)",
                                     "Wang (1996)",
                                     "UWS (various)"),
                      ppfd_data = c(mean(c(600, 800))/10^6, # Data from Kumarathunge et al. (2019) Table S1
                                    1700/10^6,
                                    2000/10^6,
                                    1800/10^6,
                                    2004/10^6,
                                    1100/10^6,
                                    1000/10^6,
                                    1000/10^6,
                                    800/10^6,
                                    1500/10^6,
                                    1400/10^6,
                                    mean(c(1500, 2000))/10^6,
                                    2000/10^6,
                                    1500/10^6,
                                    1500/10^6,
                                    1050/10^6,
                                    1800/10^6,
                                    1000/10^6,
                                    1200/10^6,
                                    1800/10^6))
        
        out <-tib %>% mutate(sitename = as.factor(sitename),
                             dataset_yr = as.factor(dataset_yr))
        
        return(tib)
    }
}


