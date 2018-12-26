suppressMessages({
  library("sleuth")
})

# Get the direcories for the samples and the corresponding kallisto results
samples_dir <- file.path('.', snakemake@input['sam'])
kal_dirs <- file.path('.', unlist(snakemake@input['kal']))

s2c <- read.table(samples_dir, header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c,
                     path = kal_dirs,
                     output_file = file.path('.', unlist(snakemake@output$'sleuth')))

so <- sleuth_prep(s2c, ~ condition, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')

so <- sleuth_lrt(so, 'reduced', 'full')

normalized <- kallisto_table(so)

# First write a table containing all samples
write.table(normalized,
            file = snakemake@output$'complete',
            quote = FALSE,
            sep = '\t',
            col.names = NA)

# Then write one table for each sample only containing its data
for (row in 1:nrow(s2c)) {
    write.table(dplyr::filter(normalized[, 1:ncol(normalized)-1],
                              sample == s2c[row, ]$'sample'),
                file = s2c[row, ]$'output_file',
                quote = FALSE,
                sep = '\t',
                col.names = NA)
}
