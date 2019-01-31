#r = getOption("repos")
#r["CRAN"] = "http://cran.us.r-project.org"
#options(repos = r)
#if(!requireNamespace("BiocManager",quietly=TRUE))
#  install.packages("BiocManager")
#if(!require("AnnotationDbi"))
#  BiocManager::install("AnnotationDbi", version = "3.8")
#if(!require("EnsDb.Hsapiens.v86"))
#  BiocManager::install("EnsDb.Hsapiens.v86", version = "3.8")
#if(!require("#"))
#  BiocManager::install("pathview", version = "3.8")
#if(!require("gage"))
#  BiocManager::install("gage", version = "3.8")
#if(!require("gageData"))
#  BiocManager::install("gageData", version = "3.8")


library("gage")
library("AnnotationDbi")
library("gageData")
library("pathview")
library("EnsDb.Hsapiens.v86")
columns(EnsDb.Hsapiens.v86)
# load the sleuth object and get the wlad test result
sleuth_table <- read.csv(as.character(snakemake@input[1]), header=TRUE, sep="\t")

#get the transcript with significant difference
sleuth_table<- na.omit(sleuth_table)
result_table <- subset(sleuth_table,select=c(target_id,pval,b))
result_table<- subset(result_table,pval<0.05)
target_id <- result_table$target_id
valid_id <- gsub("\\..*","",target_id)

#map the transcript id to entrezid
result_table$EntrezId <- mapIds(EnsDb.Hsapiens.v86, valid_id, column ="ENTREZID",keytype="TXID",multiVals="list")

#make the table with entrezid and beta value
foldChange <- result_table$b
names(foldChange)<-result_table$EntrezId


#get gene set test
data(kegg.sets.hs)
data(sigmet.idx.hs)
kegg.sets.hs <- kegg.sets.hs[sigmet.idx.hs]
gageResult <- gage(foldChange,gsets=kegg.sets.hs, same.dir=TRUE) 

# Get the pathways
upregulation10 <- data.frame(id=rownames(gageResult$greater), gageResult$greater) 
upregulation10 <- head(upregulation10,n = 10)
downregulation10 <-data.frame(id=rownames(gageResult$less),gageResult$less)
downregulation10 <-head(downregulation10,n=10)
up10Id <- substr(as.character(upregulation10$id), start=1, stop=8)
down10Id <- substr(as.character(downregulation10$id),start = 1,stop = 8)

#Define plotting function for applying later
plot_pathway = function(pid) pathview(gene.data=foldChange, pathway.id=pid, species="hsa", new.signature=FALSE,gene.idtype ="entrez" )

#plot multiple pathways (plots saved to disk and returns a throwaway list object)
snakemake@output[1]
typeof(snakemake@output[1])
as.character(snakemake@output[1])
snakemake@output[2]
upPath = sapply(up10Id, function(pid) pathview(gene.data=foldChange, pathway.id=pid, species="hsa",kegg.dir = as.character(snakemake@output[1])))
downPath=sapply(down10Id, function(pid) pathview(gene.data=foldChange, pathway.id=pid, species="hsa",kegg.dir = as.character(snakemake@output[2])))