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
                           initial_effort = 0.3, lambda = 2.1)
plotSpectra(p)

p <- matchBiomasses(p)
plotSpectra(p)

# Rescale resource to be in line with fish community
p <- alignResource(p)
plotSpectra(p)

p <- steady(p)
plotSpectra(p)

# Run the following several times until steady state is reached quickly
p <- p |> calibrateBiomass() |> matchBiomasses() |> steady()
plotSpectra(p)

# Now tune the growth rates
p <- tuneGrowth(p)

# Finally look at other aspects of the model
p <- tuneParams(p)

saveRDS(p, "curonian_model_2.rds")
