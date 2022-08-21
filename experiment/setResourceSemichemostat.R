setResourceSemichemostat <- function(params, resource_params) {
    if (!is.null(getComponent(params, "MR"))) {
        mizerMR::resource_params(params) <- resource_params
    
        rr <- mizerMR::resource_rate(params)
        cc <- (rr + getResourceMort(params))/rr * mizerMR::initialNResource(params)
        cc[rr == 0] <- 0
        mizerMR::resource_capacity(params) <- cc
    } else {
        stop("This function is not yet implemented for single-resource models.")
    }
    params
}