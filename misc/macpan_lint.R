library(lintr)
library(stringr)
library(dplyr)
library(rex)

old_regexes <- lintr:::style_regexes
regexes <- old_regexes
regexes$macpan_symbols <- rex(start, one_or_more(alpha), zero_or_more(number),
                              zero_or_more(alpha), end)

## Maybe .
## 1 or more alphanumeric, we don't care about case

macpan_snake <- rex(start,
                    one_or_more(alnum),
                    '_',
                    one_or_more(one_of('_', alnum)))

regexes$macpan_snake <- macpan_snake


assignInNamespace("style_regexes", regexes, ns = "lintr", pos = "package:lintr")

chosen_linters <- with_defaults(
  line_length_linter = line_length_linter(120),
  commented_code_linter = NULL,
  object_name_linter = object_name_linter(styles = c("macpan_symbols", "macpan_snake", "snake_case"))
)

lints <- lintr::lint_package('.', linters=chosen_linters)
assignInNamespace("style_regexes", old_regexes, ns = "lintr", pos = "package:lintr")

df <- as.data.frame(lints)
write.csv(df, file="misc/lints.csv")

print("Lints have been written to misc/lints.csv")

