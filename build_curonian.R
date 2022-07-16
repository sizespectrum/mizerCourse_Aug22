library(tidyverse)
remotes::install_github("sizespectrum/mizerExperimental")
library(mizerExperimental)

curonian_params <- read_csv("curonian_params.csv")
curonian_interaction <- read_csv("curonian_interaction.csv")

# Test that ordering of species is consistent
all(colnames(curonian_interaction) == curonian_params$species)

sp <- curonian_params |>
    select(species, a, b, w_inf, w_mat, k_vb, beta, sigma, biomass_cutoff,
           Scaled_biomass, mat_age, max_age, resilience, latinName, funcgr) |>
    rename(biomass_observed = Scaled_biomass)

p <- newMultispeciesParams(sp, curonian_interaction, no_w = 200)

plotSpectra(p)

# The following will be integrated into mizer 2.4 and released before the school

# Set initial abundances to the solutions in a fixed power-law background
p@initial_n[] <- 0
pt <- p
pt@initial_n_pp <- p@resource_params$kappa * 
    p@w_full ^ (-p@resource_params$lambda)
pt@interaction[] <- 0
income <- getEReproAndGrowth(pt) + p@metab
for (i in seq_len(nrow(p@species_params))) {
    # At small sizes the income should be A w^n. Determine A
    # Use w_min_idx + 1 in case user has implemented reduced growth
    # for the smallest size class (see e.g. #241)
    iw <- p@w_min_idx[i] + 1
    A <- income[i, iw] / (p@w[iw] ^ p@species_params$n[i])
    
    mort <- 0.4 * A * p@w ^ (p@species_params$n[i] - 1) + getFMort(p)
    growth <- getEGrowth(pt)[i, ]
    
    idxs <- p@w_min_idx[i]:(min(which(c(growth, 0) <= 0)) - 1)
    idx <- idxs[1:(length(idxs) - 1)]
    # Steady state solution of the upwind-difference scheme used in project
    p@initial_n[i, idxs] <- 
        c(1, cumprod(growth[idx] / ((growth + mort * p@dw)[idx + 1])))
}
plotSpectra(p)

p <- matchBiomasses(p)
plotSpectra(p)

# Rescale resource to be in line with fish community
# First get average fish abundance
sc <- colSums(p@initial_n) * 
    p@w ^ p@resource_params$lambda
min_w_mat = min(p@species_params$w_mat)
max_w_mat = max(p@species_params$w_mat)
sel <- (p@w > min_w_mat) & (p@w < max_w_mat)
factor <- sum(sc[sel]) / sum(sel) / resource_params(p)$kappa 

p@cc_pp <- p@cc_pp * factor
p@resource_params$kappa <- p@resource_params$kappa * factor
initialNResource(p) <- p@initial_n_pp * factor
p@search_vol = p@search_vol / factor
p@species_params$gamma <- p@species_params$gamma / factor

plotSpectra(p)
plotGrowthCurves(p, species_panel = TRUE)
plotFeedingLevel(p, include_critical = TRUE)

p <- setBevertonHolt(p, reproduction_level = 0.5)

p <- tuneParams(p)

s <- project(p)
animateSpectra(s)

ps <- steady(p)

plotSpectra(ps)

p <- tuneParams(p)

plotSpectra(p)
p <- steady(p)
