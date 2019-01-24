suppressMessages({
  library("sleuth")
})
library(ggplot2)

# Get the direcories for the samples and the corresponding kallisto results
so <- sleuth_load(file(snakemake@input$'sleuth_object'))
pca_plot <- plot_pca(so, pc_x = 1L, pc_y = 2L, use_filtered = TRUE,
                     units = 'est_counts', text_labels = TRUE,
                     color_by = NULL)

ggsave(snakemake@output$'graph_pca', dpi = 300)

# Normalize and write in tsb file
# normalized <- kallisto_table(so)
