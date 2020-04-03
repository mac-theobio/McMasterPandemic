##' given a log-linear intercept and slope and target values for R0 or GI,
##' return initial conditions and parameters that will (approximately?)
##' recapitulate the observed log-linear dynamics
##' @param int log-linear intercept (expected value on date1)
##' @param slope log-linear slope (growth rate)
##' @param R0 desired/calibrated R0
##' @param GI desired/calibrated mean generation interval (stub)
##' @param date0 date for initial conditions
##' @param date1 initial date for regression
##' @param var variable for regression (H, D, ICU)
##' @export
calibrate <- function(int,slope,pop,params,R0=3,GI=NULL,
                            date0=ldmy("1-Mar-2020"),
                            date1=ldmy("25-Mar-2020"),
                            var="H") {
    ldmy <- lubridate::dmy  ## sugar) {
    ## first adjust base params to set r and R0 (or GI) to specified value
    params[["N"]] <- pop
    state <- make_state(N=pop,E0=1)
    p2 <- fix_pars(params,state,target=c(r=slope,R0=R0))
    ## now get initial conds
    list(state=get_init(date0,date1,p2,int,slope,var),
         params=p2)
}    

