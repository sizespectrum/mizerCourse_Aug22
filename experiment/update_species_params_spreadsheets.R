sp <- read.csv(paste0("build/curonian_params.csv")) |>
    select(species, w_inf, w_mat, k_vb, beta, sigma, biomass_cutoff,
           Scaled_biomass, a, b, mat_age, max_age, resilience, latinName, funcgr) |>
    rename(biomass_observed = Scaled_biomass) |>
    mutate(sigma = 2)

write.csv(sp, row.names = FALSE, quote = FALSE, 
          paste0("build/curonian_species_params.csv"))

sp[sp$species == "roach", c("w_inf", "w_mat", "k_vb", "a", "b", "mat_age", "max_age")] <- NA
write.csv(sp, row.names = FALSE, quote = FALSE, 
          paste0("build/curonian_species_params_exercise.csv"))

sp <- sp[-4, ]
write.csv(sp, row.names = FALSE, quote = FALSE, 
          paste0("build/curonian_species_params_noRoach.csv"))
