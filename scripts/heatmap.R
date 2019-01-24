# ComplexHeatmap installieren
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ComplexHeatmap", version = "3.8")

library(ComplexHeatmap)
library(circlize)

data_pval = read.table("temp/sleuth_table.tsv", header=TRUE)
data_pval = data_pval[c("target_id", "pval")]

data_pval = dplyr::filter(data_pval, pval<0.0000001)
#data_pval

data_counts = read.table("temp/counts_normalized.tsv", header=TRUE)
data_counts = data_counts[c("target_id", "est_counts")]

merged = merge(x = data_pval, y = data_counts, by = "target_id")

mat = merged["est_counts"]

mat = data.matrix(mat)
#mat = matrix(mat, nrow=sqrt(), ncol=sqrt())

svg("results/heatmap.svg")

Heatmap(mat, cluster_rows = TRUE, cluster_columns=TRUE, clustering_distance_rows = "canberra")

dev.off()
