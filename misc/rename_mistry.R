
base_path <- file.path(here::here(), 'inst', 'params', 'mistry-cmats')
from_paths <- list.files(base_path)
from_paths <- from_paths[stringr::str_detect(from_paths, 'subnational')]

province_dict <- data.frame(
  full = c('Alberta','British_Columbia','Manitoba','New_Brunswick','Newfoundland_and_Labrador',
        'Northwest_Territories','Nova_Scotia','Nunavut','Ontario','Prince_Edward_Island',
        'Quebec','Saskatchewan','Yukon'),
  abbrev = c('AB', 'BC', 'MB', 'NB', 'NL', 'NT', 'NS', 'NU', 'ON', 'PE', 'QC', 'SK', 'YT')
)

lookup_fn <- function(long){
  return(province_dict[province_dict$full==long,]$abbrev)
}

to_paths <- c()
for(path in from_paths){
  for(province in province_dict$full){
    if (stringr::str_detect(path, province)){
      to_paths <- append(to_paths, stringr::str_replace(path, province, lookup_fn(province)))
    }
  }
}


file.rename(file.path(base_path, from_paths), file.path(base_path,to_paths))
