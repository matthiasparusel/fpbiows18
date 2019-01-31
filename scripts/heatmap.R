# ComplexHeatmap installieren
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap", version = "3.8")

library(ComplexHeatmap)
library(circlize)

# pvalues einlesen
data_pval = read.table(snakemake@input$'sle', header=TRUE)
data_pval = data_pval[c("target_id", "pval")]

# nach Relevanz filtern
data_pval = dplyr::filter(data_pval, pval<0.0000001)

# Counts einlesen
data_counts = read.table(snakemake@input$'kal', header=TRUE)
data_counts = data_counts[c("target_id", "est_counts")]

# dataframes zu counts & pvalues mergen
merged = merge(x = data_pval, y = data_counts, by = "target_id")

# Matrix, die als Heatmap dargestellt werden soll
mat = merged["est_counts"]
mat = data.matrix(mat)
mat = matrix(mat, ncol=6)

svg("results/heatmap.svg")

# Heatmap plotten
Heatmap(mat, cluster_rows = TRUE, cluster_columns=TRUE, clustering_distance_rows = "canberra")#, top_annotations=ha1, bottom_annotations=ha2)

dev.off()
