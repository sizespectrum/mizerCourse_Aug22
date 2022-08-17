library(tidyverse)
library(assertthat)
library(mizerExperimental)
library(mizerMR)

#' Calculate F curves
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' Calculates SSB, RDI, RDD, and the yield of a species for a range of
#' fishing mortalities for that species while the fishing mortalities for the
#' other species are held fixed.
#'
#' @param params An object of class `MizerParams`.
#' @param species Name of the target species
#' @param F_range A sequence of fishing mortalities at which to evaluate the
#'   SSB. If missing, it is set to
#'   `seq(from = 0, to = F_max, length.out = no_steps)`.
#' @param F_max The maximum fishing mortality. Used only if `F_range` is
#'   missing.
#' @param no_steps The number of steps to use. Only used if `F_range` is
#'   missing.
#' @param tol The `projectToSteady` function stops when the relative change in
#'   the egg production RDI over t_per years is less than tol for every species.
#' @param t_max The longest time to run project to find steady state.
#' 
#' @return A data frame with columns `F`, `SSB`, `RDI`, `RDD`, `Yield`.
#' @export
#' @family summary functions
#' @concept summary_function
#' @seealso plotSSBVsF, getYieldVsF
#' @md
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' y <- getSSBVsF(params, "Cod", F_max = 1, no_steps = 5)
#' }
getFCurves <- function(params,
                       species,
                       F_range,
                       no_steps = 10,
                       F_max,
                       tol = 0.001,
                       t_max = 100) {
    # Check parameters
    params <- validParams(params)
    species <- valid_species_arg(params, species)
    if (length(species) > 1) {
        stop("You can only make this plot for one species at a time.")
    }
    idx_species <- which(params@species_params$species == species)
    if (length(idx_species) != 1) {
        stop("Invalid species specification")
    }
    
    # Project to steady with the current fishing
    params <- projectToSteady(params, t_max = t_max, progress_bar = FALSE,
                              tol = tol)
    
    # First make a new gear for that specific species
    sp_name <- params@species_params$species[[idx_species]]
    gp <- gear_params(params)
    gp$gear <- as.character(gp$gear)
    gps <- gp$species == sp_name
    gp_extra <- gp[gps, ]
    if (nrow(gp_extra) > 1) {
        stop("This function only works in the case where the target species ",
             "is selected by a single gear only")
    }
    current_FMort <- params@initial_effort * gp_extra$catchability
    gp_extra$gear <- "tmp"
    gp_extra$catchability <- 1
    gp$catchability[gps] <- 0
    gear_params(params) <- rbind(gp, gp_extra)
    initial_effort(params)["tmp"] <- initial_effort(params)[gp$gear[gps]] # setting initial effort same as original gear
    effort <- getInitialEffort(params)
    
    if (!missing(F_max)) {
        F_range = seq(0, F_max, length.out = no_steps)
    }
    
    assert_that(is.numeric(F_range))
    sel <- F_range < current_FMort
    F_range1 <- rev(F_range[sel])
    F_range2 <- F_range[!sel]
    
    df1 <- calc_vals(params = params, F_range = F_range1, 
                     idx_species = idx_species,
                     tol = tol, t_max = t_max)
    df2 <- calc_vals(params = params, F_range = F_range2, 
                     idx_species = idx_species,
                     tol = tol, t_max = t_max)
    
    return(rbind(df1, df2))
}

#' Plot normalised SBB, RDD and Yield versus F
#'
#' @inherit getFCurves
#'
#' @inheritParams getYieldVsF
#'
#' @return A ggplot object
#' @export
#' @family plotting functions
#' @seealso getYieldVsF
#' @md
#' @examples
#' \dontrun{
#' params <- newMultispeciesParams(NS_species_params_gears, inter)
#' plotFCurves(params, "Cod")
#' }
plotFCurves <- function(params,
                        species,
                        no_steps = 10,
                        F_max,
                        F_range,
                        tol = .001,
                        t_max = 100) {
    
    curve <- getFCurves(params,
                        species = species,
                        F_range = F_range,
                        no_steps = no_steps,
                        F_max = F_max,
                        tol = tol,
                        t_max = t_max) |>
        mutate(SSB = SSB / max(SSB),
               RDI = RDI / max(RDI),
               RDD = RDD / max(RDD),
               Yield = Yield / max(Yield)) |>
        pivot_longer(!F , names_to = "Quantity", values_to = "Quantiles")
    
    ggplot(curve) +
        geom_line(aes(x = F, y = Quantiles, linetype = Quantity)) +
        xlab("Fishing mortality (1/yr)") +
        ggtitle(species)
}

#' Calculates SSB, RDD and yield at a range of fishing mortalities
#'
#' @description
#' This function replaces a loop used multiple times within
#' `getFCurves`
#'
#' @inheritParams getYieldVsF
#'
#' @return A data frame with columns `F`, `SSB`, `RDI`, `RDD`, `Yield`
#'
calc_vals <- function(params, F_range, idx_species,
                      tol = 0.001, t_max = 100) {
    SSB_vec <- F_range # To get the right length
    RDI_vec <- F_range
    RDD_vec <- F_range
    Yield_vec <- F_range
    for (i in seq_along(F_range)) {
        params@initial_effort["tmp"] <- F_range[i]
        params <- projectToSteady(params, t_max = t_max, progress_bar = FALSE,
                                  tol = tol)
        
        SSB_vec[i] <- getSSB(params)[idx_species]
        RDI_vec[i] <- getRDI(params)[idx_species]
        RDD_vec[i] <- getRDD(params)[idx_species]
        Yield_vec[i] <- getYield(params)[idx_species]
    }
    
    data.frame(F = F_range, SSB = SSB_vec, RDD = RDD_vec, RDI = RDI_vec,
               Yield = Yield_vec)
}
