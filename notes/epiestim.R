library(EpiEstim)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
load("notes/ontario_clean.RData")
incid <- (ont_recent
    %>% filter(var=="newConfirmations")
    %>% select(date,value)
    %>% rename(dates=date,I=value)
)

## see ontario_calibration.Rmd
schoolClose <- "2020-03-17"
countryClose <- "2020-03-23"
socialClose <- "2020-03-28"
bd <- as.Date(c(schoolClose,countryClose,socialClose))
e1 <- estimate_R(incid=incid,
           method="uncertain_si",
           config=make_config(list(mean_si=6,
                                   std_mean_si = 1,
                                   min_mean_si = 3,
                                   max_mean_si = 8,
                                   std_si=1.5,
                                   std_std_si=0.5,
                                   min_std_si=0.5,
                                   max_std_si=2.5,
                                   n1=100,n2=100)))
plot(e1)

epiestim_fit <- (e1$R
    %>% rename_at(vars(starts_with("Q")),
                  ~gsub("Quantile\\.0\\.([0-9]+)\\(R\\)","q\\1",.))
    %>% rename(med="Median(R)")
    ## FIXME: not 100% sure how dates line up?
    %>% mutate(date=head(sort(unique(ont_recent$date)),nrow(e1$R)))
)

# rdsave("epiestim_fit")

