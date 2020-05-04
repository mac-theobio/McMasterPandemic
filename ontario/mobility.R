load("../data/open_mobility.RData")
load("ontario_clean.RData")
## comb_sub
## restrict to dates with data
comb_sub2 <- (comb_sub
    %>% right_join(ont_all %>% select(date),by="date")
    %>% na.omit()
    %>% mutate_at("rel_activity",pmin,1)
)

plot(rel_activity ~ date, comb_sub2)
