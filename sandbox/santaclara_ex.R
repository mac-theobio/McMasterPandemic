library(tidyverse)
kain_gh_url <- "https://raw.githubusercontent.com/morgankain/COVID_interventions/master/santa_clara/covid19_sccph.csv"
is_pct <- function(x) {
    all(grepl("%",na.omit(x)))
}
make_pct <- function(x) {
    as.numeric(gsub("%","",x))
}
sc_data <- (read_csv(kain_gh_url,na=c("","-"))
    %>% mutate_if(is_pct,make_pct)
    %>% mutate_at("Date",lubridate::mdy)
)
