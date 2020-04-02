
get_giMoments <- function(params) {
    ## FIXME: assumes ICU1 model. Consider adding a test in case this changes?
    ##  (will have to rethink this once we have a structured model)
	 Rv <- get_R0(params, components=TRUE)
	 R <- sum(Rv)
	 irates <- with(as.list(params), c(lambda_a, lambda_p, lambda_m, lambda_s))

	 Gbar <- sum(Rv/irates)/R + 1/gamma
	 iww <- sum(Rv/irates^2)/R + 1/gamma^2
	 kappa <- iww/Gbar^2 - 1
}

