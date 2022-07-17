library(tidyverse)
remotes::install_github("sizespectrum/mizer")
remotes::install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)
source("helpers.R")

curonian_params <- read_csv("curonian_params.csv")
curonian_interaction <- read_csv("curonian_interaction.csv")

# Test that ordering of species is consistent
all(colnames(curonian_interaction) == curonian_params$species)

sp <- curonian_params |>
    select(species, a, b, w_inf, w_mat, k_vb, beta, sigma, biomass_cutoff,
           Scaled_biomass, mat_age, max_age, resilience, latinName, funcgr) |>
    rename(biomass_observed = Scaled_biomass)

p <- newMultispeciesParams(sp, curonian_interaction, no_w = 200,
                           initial_effort = 0.3)
plotSpectra(p)

p <- matchBiomasses(p)
plotSpectra(p)

p <- alignResource(p)
plotSpectra(p)

p <- steady(p)
plotSpectra(p)

p <- p |> calibrateBiomass() |> matchBiomasses() |> steady()
plotSpectra(p)

p <- tuneGrowth(p)

p <- tuneParams(p)



plotFeedingLevel(p, include_critical = TRUE)
plotGrowthCurves(p, species_panel = TRUE)

