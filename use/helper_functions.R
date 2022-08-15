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