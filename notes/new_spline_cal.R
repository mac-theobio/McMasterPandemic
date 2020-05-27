## make a CSV file (or a data frame somehow) that has all of the knob settings you want to use
## leave defaults to calibrate_comb() out
read.csv(text="
region      ,use_zeta,use_mobility,use_spline,spline_df,vars
NewYorkCity ,TRUE    ,TRUE        ,TRUE      ,3        ,hosp/report
NewYorkCity ,TRUE    ,TRUE        ,FALSE     ,NA       ,report
NewYorkCity ,TRUE    ,TRUE        ,TRUE      ,3        ,hosp
HudsonValley,FALSE   ,FALSE       ,TRUE      ,7        ,report
")

## x is a row of the data frame
calib_fun <- function(x) {
    dat <- full_data %>% filter(region==x$region,var %in% strsplit(x$vars,"/")[[1]])
    do.call(calibrate_comb,
            c(nlist(
                ## whatever data you need to have modified
                data=dat
            ),
            x  ## parameters getting passed straight through from the data frame
               ## to calibrate_comb()
            ))
}

    
res_list <- mclapply(seq(nrow(inputs)), ...)

# rdsave(inputs, res_list)
