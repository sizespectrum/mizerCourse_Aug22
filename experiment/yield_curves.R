plotYieldVsSize(cur_model, species = "pikeperch", catch = catch_lengths,
                x_var = "Length")

gear_params(cur_model)["pikeperch, Main", "l25"] <- 33
gear_params(cur_model)["pikeperch, Main", "l50"] <- 42

cur_model <- steady(cur_model)

gear_params(cur_model)

gear_params(cur_model)$catchability <- c(0,15,2,15,8,3,1,1,0.5,2)

cur_model <- tuneParams(cur_model, catch = catch_lengths)

getReproductionLevel(cur_model)
cur_model <- setBevertonHolt(cur_model, reproduction_level = 0.7)
plotYieldVsF(cur_model, species = "predator_fish", F_max = 1, no_steps = 20, tol = 0.001)
plotYieldVsF(cur_model, species = "roach", F_max = 10, no_steps = 20)


cur_model <- setBevertonHolt(cur_model, reproduction_level = 0.9)
plotYieldVsF(cur_model, species = "breams", F_max = 2)
plotYieldVsF(cur_model, species = "roach", F_max = 10, no_steps = 20)

cur_model <- setBevertonHolt(cur_model, reproduction_level = 0.1)
plotYieldVsF(cur_model, species = "breams", F_max = 2)
plotYieldVsF(cur_model, species = "roach", F_max = 10, no_steps = 20)
