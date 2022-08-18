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


##########------------------------------------------------------
########## Calculate error function for time series calibration 
##########------------------------------------------------------

getErrorTime <- function(vary,params,dat,env=state,tol = 0.1) {
    
    params@species_params$R_max[1:12]<-10^vary[1:12]
    params@species_params$erepro[1:12]<-vary[13:24]
    params@resource_params$kappa<-10^vary[25]
    params@resource_params$r_pp<-vary[26]
    
    params <- setParams(params)
    # run to steady state and update params
    # env$params<- projectToSteady(env$params, distance_func = distanceSSLogN,
    #                 tol = tol, t_max = 200,return_sim = F)
    
    params_steady<- projectToSteady(params, distance_func = distanceSSLogN,
                                    tol = tol, t_max = 200,return_sim = F)
    
    #run time-varying effort model tthough time with new erepro
    
    simt <- project(params_steady, effort = effort,initial_n =  params_steady@initial_n, initial_n_pp = params_steady@initial_n_pp)
    
    # get biomass through time
    biomass <- sweep(simt@n, 3, simt@params@w * simt@params@dw, "*")
    
    #get yield through time from model:
    
    f_gear<-getFMortGear(params,effort)
    yield_species_gear <- apply(sweep(f_gear, c(1, 3, 4), biomass, "*"),
                                c(1, 2, 3), sum)
    yield_species_gear
    
    yield_species <-apply(yield_species_gear, c(1, 3), sum)
    
    yield_frame <- melt(yield_species)
    
    # leave out spin up and change units to tonnes    
    y<-yield_frame[yield_frame$time >= 1947,]
    
    # disregard zeroes - these were NAs only filled in to run the model   
    
    obs<-dat$value[which(dat$value>0)]/1e3   
    pred<-y$value[which(dat$value>0)]/1e9
    
    # sum of squared errors, could use  log-scale of predictions and data (could change this or use other error or likelihood options)
    
    error <- sum((log(pred) - log(obs))^2,na.rm=T)
    
    # can use a strong penalty on the error to ensure we reach a minimum of 10% of the data (biomass or catch) for each species
    # if(any(pred < 0.1*dat)) discrep <- discrep + 1e10
    
    return(error)
    
}

#vary<-c(log10(simt@params@species_params$R_max),simt@params@species_params$erepro,log10(5e11),4)

# ## test it
# err<-getErrorTime(vary = vary, params = params, dat = obsy)
# 
# 
# err

##########------------------------------------------------------
########## Plot the outputs of the time series calibration 
##########------------------------------------------------------


plotFittedTime<-function(sim=simt,obsy=obsy,allSpecies=T,plotSpecies=NULL,startyr=1947,endyr=2019){
    
    biomass <- sweep(sim@n, 3, sim@params@w * sim@params@dw, "*")
    params<-sim@params
    effort<-sim@effort
    
    f_gear<-getFMortGear(params,effort)
    yield_species_gear <- apply(sweep(f_gear, c(1, 3, 4), biomass, "*"),
                                c(1, 2, 3), sum)
    yield_species_gear
    
    yield_species <-apply(yield_species_gear, c(1, 3), sum)
    
    yield_frame <- melt(yield_species)
    
    
    # output modelled yields and reshape for plotting - dont know why built-in getYield function doesn't woprk
    
    # y <- getYield(simt)
    # y <- reshape2::melt(y)
    
    y<-yield_frame[yield_frame$time >= startyr,]
    
    
    # plot these
    
    if (allSpecies ==T) { 
        p<-ggplot(y) + geom_line(data=y, aes(x = time, y = (value)/1e6, 
                                             colour = sp)) +
            geom_point(data=obsy,aes(x = time, y = (value), 
                                     colour = sp),size=0.1) +
            facet_wrap(~sp,scales="free_y") +
            scale_y_continuous(name = "yield [g/year]")  +
            scale_colour_manual(values = sim@params@linecolour) +
            xlim(startyr, endyr)
    }
    
    # look only at  one species at a time and examine on linear 
    if (allSpecies ==F){
        p<-ggplot(y) + geom_line(data=filter(y,sp==plotSpecies), aes(x = time, y = value/1e6,colour = sp)) +
            geom_point(data=filter(obsy,sp=="Cod"),aes(x = time, y = value, 
                                                       colour = sp),size=0.6) +
            #facet_wrap(~sp) +
            scale_y_continuous(name = "Yield [g/year]")  +
            scale_colour_manual(values = sim@params@linecolour) +
            xlim(startyr, endyr)
    }
    
    return(p)
}

### RF: Set of functions using size selectivity of the data


# this function adds a lower boundary to selected size
plotBiomassObservedVsModelCustom <- function (object, species = NULL, ratio = FALSE, log_scale = TRUE, 
                                              return_data = FALSE, labels = TRUE, show_unobserved = FALSE) 
{
    if (is(object, "MizerSim")) {
        params = object@params
        n <- finalN(object)
    }
    else if (is(object, "MizerParams")) {
        params = object
        n <- initialN(params)
    }
    else {
        stop("You have not provided a valid mizerSim or mizerParams object.")
    }
    sp_params <- params@species_params
    species = valid_species_arg(object, species)
    if (length(species) == 0) 
        stop("No species selected, please fix.")
    row_select = match(species, sp_params$species)
    if (!"biomass_observed" %in% names(sp_params)) {
        stop("You have not provided values for the column 'biomass_observed' ", 
             "in the mizerParams/mizerSim object.")
    }
    else if (!is.numeric(sp_params$biomass_observed)) {
        stop("The column 'biomass_observed' in the mizerParams/mizerSim object", 
             " is not numeric, please fix.")
    }
    else {
        biomass_observed = sp_params$biomass_observed
    }
    
    cutoffLow <- sp_params$biomass_cutoffLow[row_select]
    if (is.null(cutoffLow)) {
        cutoffLow = rep(0, length(species))
    }
    else if (!is.numeric(cutoffLow)) {
        stop("params@species_params$biomass_cutoffLow is not numeric, \",\n                 \"please fix.")
    }
    cutoffLow[is.na(cutoffLow)] <- 0
    
    cutoffHigh <- sp_params$biomass_cutoffHigh[row_select]
    if (is.null(cutoffHigh)) {
        cutoffHigh = rep(0, length(species))
    }
    else if (!is.numeric(cutoffHigh)) {
        stop("params@species_params$biomass_cutoffHigh is not numeric, \",\n                 \"please fix.")
    }
    cutoffHigh[is.na(cutoffHigh)] <- 0
    
    sim_biomass = rep(0, length(species))
    for (j in 1:length(species)) {
        sim_biomass[j] = sum((n[row_select[j], ] * params@w * 
                                  params@dw)[params@w >= cutoffLow[j] & cutoffHigh[j] >= params@w])
    }
    dummy = data.frame(species = species, model = sim_biomass, 
                       observed = biomass_observed[row_select]) %>% mutate(species = factor(species, 
                                                                                            levels = species), is_observed = !is.na(observed) & observed > 
                                                                               0, observed = case_when(is_observed ~ observed, !is_observed ~ 
                                                                                                           model), ratio = model/observed)
    if (sum(dummy$is_observed) == 0) {
        cat(paste("There are no observed biomasses to compare to model,", 
                  "only plotting model biomasses.", sep = "\n"))
    }
    if (!show_unobserved) {
        dummy <- filter(dummy, is_observed)
    }
    if (return_data == TRUE) 
        return(dummy)
    tre <- round(sum(abs(1 - dummy$ratio)), digits = 3)
    caption <- paste0("Total relative error = ", tre)
    if (any(!dummy$is_observed)) {
        caption <- paste(caption, "\n Open circles represent species without biomass observation.")
    }
    if (ratio == FALSE) {
        gg <- ggplot(data = dummy, aes(x = observed, y = model, 
                                       colour = species, shape = is_observed)) + geom_abline(aes(intercept = 0, 
                                                                                                 slope = 1), colour = "purple", linetype = "dashed", 
                                                                                             size = 1.3) + geom_point(size = 3) + labs(y = "model biomass [g]") + 
            coord_cartesian(ylim = range(dummy$model, dummy$observed))
    }
    else {
        gg <- ggplot(data = dummy, aes(x = observed, y = ratio, 
                                       colour = species, shape = is_observed)) + geom_hline(aes(yintercept = 1), 
                                                                                            linetype = "dashed", colour = "purple", 
                                                                                            size = 1.3) + geom_point(size = 3) + labs(y = "model biomass / observed biomass") + 
            coord_cartesian(ylim = range(dummy$ratio))
    }
    gg <- gg + labs(x = "observed biomass [g]", caption = caption) + 
        scale_colour_manual(values = getColours(params)[dummy$species]) + 
        scale_shape_manual(values = c(`TRUE` = 19, `FALSE` = 1)) + 
        guides(shape = "none")
    if (log_scale == TRUE & ratio == FALSE) {
        gg = gg + scale_x_log10() + scale_y_log10()
    }
    if (log_scale == TRUE & ratio == TRUE) {
        gg = gg + scale_x_log10()
    }
    if (labels == TRUE) {
        gg = gg + ggrepel::geom_label_repel(aes(label = species), 
                                            box.padding = 0.35, point.padding = 0.5, segment.color = "grey50", 
                                            show.legend = FALSE, max.overlaps = Inf, seed = 42)
    }
    gg
}




# adapting cutoff here too
calibrateBiomassCustom <- function (params) 
{
    if ((!("biomass_observed" %in% names(params@species_params))) || 
        all(is.na(params@species_params$biomass_observed))) {
        return(params)
    }
    no_sp <- nrow(params@species_params)
    
    cutoffLow <- params@species_params$biomass_cutoffLow
    if (is.null(cutoffLow)) 
        cutoffLow <- rep(0, no_sp)
    cutoffLow[is.na(cutoffLow)] <- 0
    
    cutoffHigh <- params@species_params$biomass_cutoffHigh
    if (is.null(cutoffHigh)) 
        cutoffHigh <- rep(0, no_sp)
    cutoffHigh[is.na(cutoffHigh)] <- 0
    
    observed <- params@species_params$biomass_observed
    observed_total <- sum(observed, na.rm = TRUE)
    sp_observed <- which(!is.na(observed))
    model_total <- 0
    for (sp_idx in sp_observed) {
        model_total <- model_total + sum((params@initial_n[sp_idx, 
        ] * params@w * params@dw)[params@w >= cutoffLow[sp_idx] & cutoffHigh[sp_idx] >= params@w])
    }
    scaleModel(params, factor = observed_total/model_total)
}

# same as above
matchBiomassCustom <- function (params, species = NULL) 
{
    if (!("biomass_observed" %in% names(params@species_params))) {
        return(params)
    }
    species <- valid_species_arg(params, species = species, return.logical = TRUE) & 
        !is.na(params@species_params$biomass_observed) & params@species_params$biomass_observed > 
        0
    for (sp in (1:nrow(params@species_params))[species]) {
        cutoffLow <- params@species_params$biomass_cutoffLow[[sp]]
        if (is.null(cutoffLow) || is.na(cutoffLow)) {
            cutoffLow <- 0
        }
        cutoffHigh <- params@species_params$biomass_cutoffHigh[[sp]]
        if (is.null(cutoffHigh) || is.na(cutoffHigh)) {
            cutoffHigh <- 0
        }
        
        total <- sum((params@initial_n[sp, ] * params@w * params@dw)[params@w >= cutoffLow & cutoffHigh >= params@w])
        factor <- params@species_params$biomass_observed[[sp]]/total
        params@initial_n[sp, ] <- params@initial_n[sp, ] * factor
    }
    setBevertonHolt(params)
}

####### Compare ecosystem states
plot_relative_biomass = function(sim0, sim1, ratio = FALSE) {
    
    # assume sim0 is steady state, sim1 is some kind of variation such as fishing
    # mean of last five years
    fish_sim <- apply(N(sim1)[(dim(N(sim1))[1]-5):(dim(N(sim1))[1]),,],2:3,mean)
    unfish_sim <- apply(N(sim0)[(dim(N(sim0))[1]-5):(dim(N(sim0))[1]),,],2:3,mean)
    
    if (ratio) {
        relative_n <- melt(fish_sim/unfish_sim) # Julia's original calculation
    } else {
        relative_n <- melt((fish_sim - unfish_sim) / (fish_sim + unfish_sim)) # Gustav's suggested calculation
    }
    colnames(relative_n)[1] <- "Species"
    legend_levels <- intersect(names(sim0@params@linecolour), relative_n$Species)
    p <- ggplot(relative_n) +
        geom_line(aes(x = w, y = value, colour = Species), size = 1) +
        scale_x_continuous(trans = "log10", name = "Weight [g]") +
        scale_color_manual(values = sim0@params@linecolour[legend_levels]) +
        theme(legend.key = element_rect(fill = "white"))
    
    if (ratio == T) {
        p = p + scale_y_continuous(trans = "log10") +
            geom_hline(yintercept = 1, linetype = 1, colour="dark grey", size=0.75) +
            labs(y="Relative abundance")
    } else {
        p = p + geom_hline(yintercept = 0, linetype = 1, colour="dark grey", size=0.75) +
            labs(y="Relative difference")
    }
    
    print(p)
    
}



plot_relative_biomass<-function(sim0, sim1, ratio = FALSE) {
    
    # assume sim0 is steady state, sim1 is also steady state
    fish_sim <- N(sim1)[2,,]
    unfish_sim <- N(sim0)[2,,]
    if (ratio) {
        relative_n <- melt(fish_sim/unfish_sim) # Julia's original calculation
    } else {
        relative_n <- melt((fish_sim - unfish_sim) / (fish_sim + unfish_sim)) # Gustav's suggested calculation
    }
    colnames(relative_n)[1] <- "Species"
    legend_levels <- intersect(names(sim0@params@linecolour), relative_n$Species)
    p <- ggplot(relative_n) +
        geom_line(aes(x = w, y = value, colour = Species), size = 1) +
        scale_x_continuous(trans = "log10", name = "Weight [g]") +
        scale_color_manual(values = sim0@params@linecolour[legend_levels]) +
        theme(legend.key = element_rect(fill = "white"))
    
    if (ratio == T) {
        p = p + scale_y_continuous(trans = "log10") +
            geom_hline(yintercept = 1, linetype = 1, colour="dark grey", size=0.75) +
            labs(y="Relative abundance")
    } else {
        p = p + geom_hline(yintercept = 0, linetype = 1, colour="dark grey", size=0.75) +
            labs(y="Relative difference")
    }
    
    print(p)
    
}
