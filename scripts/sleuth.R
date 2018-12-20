suppressMessages({
  library("sleuth")
})


sample_id <- dir(file.path(".", "results"))

kal_dirs <- file.path(".", "results", sample_id, "kallisto")

s2c = read.table(file.path(".","data", "samples.tsv"), 
                 header = TRUE, 
                 stringsAsFactors = FALSE)

s2c = dplyr::mutate(s2c, path = kal_dirs)

so <- sleuth_prep(s2c, extra_bootstrap_summary = TRUE)

