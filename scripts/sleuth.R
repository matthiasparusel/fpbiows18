suppressMessages({
  library("sleuth")
})

# Get the direcories for the samples and the corresponding kallisto results
samples_dir <- file.path(snakemake@input['sam'])
kal_dirs <- file.path('.', unlist(snakemake@input['kal']))

s2c <- read.table(samples_dir, header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c,
                     path = kal_dirs,
                     output_file = file.path('.', unlist(snakemake@output$'sleuth')))

so <- sleuth_prep(s2c, ~ condition, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

# Normalize and write in tsb file
normalized <- kallisto_table(so)
write.table(normalized,
            file = snakemake@output$'counts_norm',
            quote = FALSE,
            sep = '\t',
            col.names = NA)

# Write in tsv file for use with Python
sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table,
            file = snakemake@output$'sleuth_res',
            quote = FALSE,
            sep = '\t',
            col.names = NA)
