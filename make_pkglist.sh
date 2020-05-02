#!/bin/bash
## find packages referred to by library()/require() in R/Rmd/Rnw code
find . -type f -name "*.R*" -exec egrep -H "library\(|require\(" {} \; | \
    sed -e 's/^.*\(library\|require\)//' | \
    sed -e 's/[,)"].*$//' | \
    sed -e 's/^[(]//' | \
    sort | uniq
## omit 'pkg' (false positive)
## add suggests/imports/depends from DESCRIPTION ?
