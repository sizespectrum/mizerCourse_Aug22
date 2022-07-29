library(tidyverse)
#remotes::install_github("sizespectrum/mizer")
remotes::install_github("sizespectrum/mizerExperimental")
remotes::install_github("sizespectrum/mizerMR")
library(mizerExperimental)
library(mizerMR)
source("experiment/helpers.R")

curonian_params <- read_csv("experiment/curonian_params.csv")
curonian_interaction <- read_csv("experiment/curonian_interaction.csv")

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

#saveRDS(p, "experiment/curonian_model_2.rds")


p <- readRDS("experiment/curonian_model_2.rds")
### Asta continues

# #number of size groups 
# no_size_groups = 200
# #timestep used in the integration 
# dt = 0.2

# p <- newMultispeciesParams(sp, curonian_interaction, no_w = 200,
#                            initial_effort = 0.3, lambda = 2.1)

## resource params
kappa =  p@resource_params$kappa*3 #2#1 #2 # 20 # 20 # intercept assuming g/m2
lambda = p@resource_params$lambda #2.1 # 
w_pp_cutoff = p@resource_params$w_pp_cutoff/5 # this is now 10, should be lower
r_pp = p@resource_params$r_pp #default is 10 and still gets depleted 
min_w_pp =  1e-10 #g

kappa_ben = (kappa *9) #2 #8#4 # 8 #80 #80 # intercept assuming g/m2  
lambda_ben = 1.9 # 1.85 #this slope does not include urchins and lobsters
w_bb_cutoff = 10 #
r_bb = 5 # 1.5 # something to be calibrated. Default mizer option is 10
min_w_bb = 0.001 # 0.01

#p <- newMultispeciesParams(sp, curonian_interaction, no_w = 200,
#                           initial_effort = 0.3, lambda = 2.1)

resource_params(p) <- data.frame(
    resource = c("pl", "bb"),
    lambda = c(lambda, lambda_ben),
    kappa = c(kappa, kappa_ben),
    r_pp = c(r_pp, r_bb),
    w_min = c(min_w_pp, min_w_bb),
    w_max = c(w_pp_cutoff, w_bb_cutoff)
)



#now our default availabilities are 1 for everything 

# smallfish      1  1
# ruffe          1  1
# breams         1  1
# roach          1  1
# vimba          1  1
# carassius      1  1
# perch          1  1
# pikeperch      1  1
# burbot         1  1
# predator_fish  1  1

plankton_avail <- c(0.7, 0.7, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3)
benthos_avail <- c(0.2, 0.5, 0.7, 0.7, 0.7, 0.7, 0.3, 0.3, 0.3, 0.3)

resource_interaction(p)[, 1] <- plankton_avail
resource_interaction(p)[, 2] <- benthos_avail

plotSpectra(p)

#p <- matchBiomasses(p)
#plotSpectra(p)

# Rescale resource to be in line with fish community
#p <- alignResource(p)
#plotSpectra(p)

p <- steady(p)
plotSpectra(p)

# Run the following several times until steady state is reached quickly
p <- p |> calibrateBiomass() |> matchBiomasses() |> steady()
plotSpectra(p)

# Now tune the growth rates
p <- tuneGrowth(p)

# Finally look at other aspects of the model
p <- tuneParams(p)

saveRDS(p, "experiment/curonian_model_2res_a.rds")

#predators still have too much benthos, reduce availability further 



