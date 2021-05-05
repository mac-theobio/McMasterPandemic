library("covr")
# For some reason, there seems to be no programmatic way that works for everyone to set
# working directories to the source file location (the location of this file)

base_dir <- paste(getwd(), "/../../", sep = "")
print(base_dir)

cov <- package_coverage(base_dir)

df <- as.data.frame(cov)
write.csv(df, paste(base_dir, "/refactor/coverage/coverage_details.csv", sep = ""))

no_cov <- zero_coverage(cov)
no_cov_df <- as.data.frame(no_cov)
write.csv(no_cov_df, paste(base_dir, "/refactor/coverage/no_coverage_details.csv", sep = ""))

# For some reason, sink doesn't work with ... colors? I just copy-pasted the print into coverage_report.txt
print(cov)
