suppressMessages({
  library("sleuth")
})

sample_id <- dir(file.path(".", "results"))

kal_dirs <- file.path(".", "results", sample_id, "kallisto")

s2c = read.table(file.path(".","data", "samples.tsv"), 
                 header = TRUE, 
                 stringsAsFactors = FALSE)

s2c = dplyr::mutate(s2c, path = kal_dirs)

so = sleuth_prep(s2c, ~ condition, extra_bootstrap_summary = TRUE)

so = sleuth_fit(so)
so = sleuth_fit(so, ~1, 'reduced')

so = sleuth_lrt(so, 'reduced', 'full')

normalized = kallisto_table(so)

write.table(normalized, file='test.tsv', quote=FALSE, sep='\t', col.names = NA)

for (id in sample_id) {
  write.table(dplyr::filter(normalized, sample == id),
              file = paste(id, '.tsv', sep=''), 
              quote = FALSE, 
              sep = '\t', 
              col.names = NA)
}
