source("experiment/getFCurves.R")

cur_model <- readRDS("build/cur_model_refined.rds") |>
    steady()

cur_model <- setBevertonHolt(cur_model, reproduction_level = 0.2)
plotFCurves(cur_model, species = "breams", F_max = 0.5, no_steps = 20)

cur_model <- setBevertonHolt(cur_model, reproduction_level = 0.5)
plotFCurves(cur_model, species = "predator_fish", F_max = 0.5, no_steps = 20)

F_range <- c(seq(0, 0.3, 0.02), seq(0.3, 1, 0.1))
plotFCurves(cur_model, species = "predator_fish", F_range = F_range)

gear_params(cur_model)
