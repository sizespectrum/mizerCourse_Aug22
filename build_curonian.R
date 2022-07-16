library(tidyverse)
library(mizerExperimental)

# Test that ordering of species is consistent
all(colnames(curonian_interaction) == curonian_params$species)

sp <- curonian_params |>
    select(species, a, b, w_inf, w_mat, k_vb, beta, sigma, biomass_cutoff,
           Scaled_biomass, mat_age, max_age, resilience, latinName, funcgr) |>
    rename(biomass_observed = Scaled_biomass)

p <- newMultispeciesParams(sp, curonian_interaction, no_w = 200)

p <- tuneParams(p)

plotSpectra(p)
p <- steady(p)
