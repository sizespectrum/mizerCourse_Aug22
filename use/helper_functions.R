plotBevertonHolt <- function(params, species) {
    select <- species_params(params)$species == species
    erepro <- species_params(params)$erepro[select]
    w0 <- params@w[params@w_min_idx[select]]
    E_R_ss <- getRDI(params)[select] / erepro * 2 * w0
    R_dd_ss <- getRDD(params)[select]
    R_max  <- species_params(params)$R_max[select]
    E_R <- seq(0, 2 * E_R_ss, length.out = 50)
    R_di = erepro * E_R / 2 / w0
    R_dd <- R_di / (1 + R_di / R_max)
    df <- melt(data.frame(E_R, R_dd, R_di, R_max), id.vars = "E_R")
    ggplot(df) +
        geom_line(aes(x = E_R, y = value, linetype = variable)) +
        geom_point(aes(x = E_R_ss, y = R_dd_ss), size = 3, color = "red") +
        ylim(NA, 1.1 * R_max) +
        ylab("Reproduction rate [eggs/year]") +
        xlab("Energy invested [g/year]")
}

plotBevertonHolt2 <- function(params, params2, species) {
    select <- species_params(params)$species == species
    erepro <- species_params(params)$erepro[select]
    w0 <- params@w[params@w_min_idx[select]]
    E_R_ss <- getRDI(params)[select] / erepro * 2 * w0
    R_dd_ss <- getRDD(params)[select]
    E_R <- seq(0, 2 * E_R_ss, length.out = 50)
    
    R_max  <- species_params(params)$R_max[select]
    R_di = erepro * E_R / 2 / w0
    R_dd <- R_di / (1 + R_di / R_max)
    df <- melt(data.frame(E_R, R_dd, R_di, R_max), id.vars = "E_R")
    df$Model <- "Model 1"
    
    erepro <- species_params(params2)$erepro[select]
    R_max  <- species_params(params2)$R_max[select]
    R_di = erepro * E_R / 2 / w0
    R_dd <- R_di / (1 + R_di / R_max)
    df2 <- melt(data.frame(E_R, R_dd, R_di, R_max), id.vars = "E_R")
    df2$Model <- "Model 2"
    
    ggplot(rbind(df, df2)) +
        geom_line(aes(x = E_R, y = value, linetype = variable,
                      colour = Model, size = Model)) +
        geom_point(aes(x = E_R_ss, y = R_dd_ss), size = 3, color = "red") +
        ylim(NA, 1.1 * R_max) +
        ylab("Reproduction rate [eggs/year]") +
        xlab("Energy invested [g/year]") +
        labs(linetype = "", size = "R_max", colour = "R_max") +
        scale_size_manual(values = c(0.5, 1)) +
        scale_colour_manual(values = c("blue", "black")) +
        scale_linetype_manual(values = c("solid", "dashed", "dotted"))
}