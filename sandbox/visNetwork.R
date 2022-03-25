library(visNetwork)

state_nms = names(model$state)
grps = sub(".*_", "", state_nms)
epi = sub("_.*", "", state_nms)
grp_nums = case_when(grps == 'unvax' ~ 1, grps == 'vaxdose1' ~ 2, grps == 'vaxprotect1' ~ 3, grps == 'vaxdose2' ~ 4, grps == 'vaxprotect2' ~ 5)
shape = c("ellipse", "circle", "database", "box", "text")[grp_nums]
nn = visNetwork(
  data.frame(id = state_nms,
             label = epi,
             title = state_nms,
             group = grps,
             y = grp_nums,
             physics = FALSE),
  select(rate_summary(model), from, to) %>%
    mutate(arrows = "to")) %>%
  visLegend() %>%
  visPhysics(stabilization = FALSE)
  #visHierarchicalLayout(direction = "UD")
nn
