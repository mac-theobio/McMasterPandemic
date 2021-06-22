library(styler)
library(dplyr)
library(stringr)

replace_comment_with_dbl <- function(text) {
    ## this looks safe (only replace the first # in a comment?)
    str_replace(text, "#+", "##")
}
set_double_comments <- function(pd_flat) {
    op <- pd_flat$token %in% "COMMENT"
    pd_flat$text[op] <- sapply(pd_flat$text[op], replace_comment_with_dbl)
    return(pd_flat)
}

new_pipes <- function(pd_flat) {
    op <- Reduce("&", list(pd_flat$token %in% "SPECIAL-PIPE", pd_flat$lag_newlines > 0L))
    if (!all(op == FALSE)) {
        pos <- c(1:length(pd_flat$token))[op]
        pd_flat$indent[pos] <- pd_flat$indent[pos + 1]
    }
    return(pd_flat)
}


macpan_style <- tidyverse_style(indent_by = 4)
macpan_style$token$set_double_comments <- set_double_comments
## WORKAROUND: this allows the indentation of newline pipes to be automatically corrected.
macpan_style$token$new_pipes <- new_pipes
style_pkg(".", transformers = macpan_style)
style_dir("inst", transformers = macpan_style)
style_dir("misc", transformers = macpan_style)
