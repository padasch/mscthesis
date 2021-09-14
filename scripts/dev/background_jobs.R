# START OF SCRIPT ----

# SOURCE ALL ----
suppressMessages(source("~/projects/mscthesis/scripts/final/source_fct_pkg.R"))

# ADD SCRIPT ----
all_tleaf_obs_smith_lowlight   <- plot_topt_tgrowth(y = "tc_opt_obs", forcing = "growth", limitation = "smith37")
all_tleaf_sim_smith_lowlight   <- plot_topt_tgrowth(y = "tc_opt_sim", forcing = "growth", limitation = "smith37")

## HIGH LIGHT TC_GROWTH_LEAF PREDICTIONS
all_tleaf_obs_smith_highlight  <- plot_topt_tgrowth(y = "tc_opt_obs", forcing = "metainfo", limitation = "smith37")
all_tleaf_sim_smith_highlight  <- plot_topt_tgrowth(y = "tc_opt_sim", forcing = "metainfo", limitation = "smith37")

# all_inst_smith_meas_with_eb   <- xxx_fct(smith_or_farq = "smith37",    ppfd_growth_or_meas = "metainfo", with_eb_models = T)
# all_inst_smith_meas_with_eb   <- xxx_fct(smith_or_farq = "farquhar89",    ppfd_growth_or_meas = "metainfo", with_eb_models = T)
# 
# all_inst_smith_growth_with_eb   <- xxx_fct(smith_or_farq = "smith37",    ppfd_growth_or_meas = "growth", with_eb_models = T)
# all_inst_smith_growth_with_eb   <- xxx_fct(smith_or_farq = "farquhar89",    ppfd_growth_or_meas = "growth", with_eb_models = T)
# END OF SCRIPT ----
beep(3)
