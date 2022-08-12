load("experiment/growth_data/breamsCombined_sizeAge.RData")
load("experiment/growth_data/perch_sizeAge.RData")
load("experiment/growth_data/pikeperch_sizeAge.RData")
load("experiment/growth_data/roach_sizeAge.RData")
load("experiment/growth_data/ruffe_sizeAge.RData")

# combine into a single data frame
size_at_age <- rbind(bream, perch, pikeperch, roach, ruffe) |>
    select(species, age, TL) |>
    rename(length = TL)
# aggregate all breams
size_at_age$species[size_at_age$species == "breams_2"] <- "breams"

# Load species parameters
sp <- read.csv("build/curonian_species_params.csv")
# If the von Bertalanffy t0 parameter is missing, set to 0
if (!hasName(sp, "t0")) {
    sp$t0 <- 0
}
# convert weights to lengths
sp <- sp |>
    mutate(l_mat = (w_mat / a) ^ (1/b),
           l_inf = (w_inf / a) ^ (1/b)) |>
    select(species, l_mat, l_inf, mat_age, max_age, k_vb, t0) 

# create data frame with von Bertalanffy curves
vb <- expand.grid(species = sp$species, age = 0:20) |>
    left_join(sp, by = "species") |>
    mutate(length = l_inf * (1 - exp(-k_vb * (age - t0))))

ggplot() +
    geom_point(aes(x = age, y = length), data = size_at_age,
               position = "jitter", alpha = 0.2) +
    geom_line(aes(x = age, y = length), data = vb) +
    geom_hline(aes(yintercept = l_inf), data = sp) +
    geom_hline(aes(yintercept = l_mat), data = sp) +
    geom_vline(aes(xintercept = max_age), data = sp) +
    geom_vline(aes(xintercept = mat_age), data = sp) +
    facet_wrap(vars(species), scales = "free_y")
