
library(tidyverse)

load("experiment/catch_size.RData")

catch_size <- catch_size |>
    rename(species = species_name) |>
    drop_na()
head(catch_size)

summary(catch_size$total_length)
breaks <- seq(14, 68, 2)
length <- breaks[-length(breaks)]

catch_with_bins <- catch_size |>
    mutate(bin = cut(total_length, breaks = breaks,
                     include.lowest = TRUE,
                     labels = FALSE))
head(catch_with_bins)

catch <- catch_with_bins |>
    group_by(bin, species) |>
    summarise(catch = n()) |>
    mutate(dl = 2,
           length = dl * bin + 12) |>
    ungroup() |>
    select(species, catch, length, dl)
head(catch)

ggplot(catch) +
    geom_line(aes(x = length, y = catch, colour = species))

saveRDS(catch, "build/catch.rds")
