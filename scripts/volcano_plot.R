library("sleuth")

so <-sleuth_load(file.path(snakemake@input[1]))
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <-sleuth_wt(so, colnames(so$design_matrix)[2])
tests(so)
plot_volcano(so,colnames(so$design_matrix)[2],"wt")

svg(file.path(snakemake@output[1]))

plot_volcano(so,colnames(so$design_matrix)[2],"wt")
dev.off()
