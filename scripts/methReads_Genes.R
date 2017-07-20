library(jsonlite)
library(RColorBrewer)
library(rtracklayer)
library("ComplexHeatmap")
library("reshape")
resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}
araport.gff <- import("~/ARAPORT11/Araport11_GFF3_genic_regions.genes.TEs.201606.gff")
araport.genes <- read.table("~/ARAPORT11/Araport11_GFF3_genes_201606.bed", as.is = T)


output_fol <- "/lustre/scratch/projects/cegs/rahul/016.bshap/005.mutants.GSE39901/001.methReads_aragenes/"
setwd(output_fol)

init_cls = c("chr", "start", "end", "(0,0,0)","(1,0,0)","(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")

### Checking each single file......
plot_heatmap <- function(methreads.dir, methreads.file, i, cluster_wt=NULL, req_cols = NULL, top = NULL){
  test.mut <- read.csv(file = file.path(methreads.dir, as.character(methreads.file)[i]), header = F)
  colnames(test.mut) <- get("init_cls", envir = .GlobalEnv)
  #req.mat <- na.omit(test.mut[,c(5,6,7,8,9,10,11)]/rowSums(test.mut[,c(5,6,7,8,9,10,11)]))  ## Gene fractions filtering non methylated
  #### Just taking the number of reads as logarthimic number
  req.mat <- log10(na.omit(test.mut[,c(4,5,6,7,8,9,10,11)]))  ## Gene fractions filtering non methylated
  req.mat[req.mat == -Inf] = 0
  ## Taking fracion now with the total reads and logarithmic
  req.mat <- log10(na.omit(test.mut[,c(4,5,6,7,8,9,10,11)]/rowSums(test.mut[,c(4,5,6,7,8,9,10,11)])))  ## Gene fractions filtering non methylated
  req.mat[req.mat == -Inf] = 0
  
  #req.mat <- na.omit((sum(rowSums(test.mut[,c(4,5,6,7,8,9,10,11)]))*test.mut[,c(5,6,7,8,9,10,11)])/rowSums(test.mut[,c(4,5,6,7,8,9,10,11)]))  ### Conditional probability on the number of reads matching
  #
  if(!is.null(top)){
    req.mat <- req.mat[1:top,]
  } else{
    req.mat <- req.mat[,]
  }
  if(!is.null(req_cols)){
    req.mat <- req.mat[,req_cols]
  }
  
  if(!is.null(cluster_wt)){
    draw(Heatmap(req.mat, cluster_rows=cluster_wt,cluster_columns =F, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_rows = "binary", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = "Log10 \nnumber of reads", color_bar = "continuous"), row_title = paste("Reads mapped to genes in", names(methreads.file)[i])))
  } else{
    draw(Heatmap(req.mat, cluster_rows=T,cluster_columns =F, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_rows = "binary", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = "Log10 \nnumber of reads", color_bar = "continuous"), row_title = paste("Reads mapped to genes in", names(methreads.file)[i])))
  }
  
}

mutant.data$Run_s[which(mutant.data$genotype_s == "cmt2")]

methreads.dir <- "/projects/cegs/rahul/016.bshap/005.mutants.GSE39901/001.methReads_aragenes/"

methreads.file <- c("meths.SRR534177.summary.txt", "meths.SRR534239.summary.txt", "meths.SRR534240.summary.txt", "meths.SRR534182.summary.txt", "meths.SRR534215.summary.txt", "meths.SRR534223.summary.txt", "meths.SRR869314.summary.txt")
methreads.file <- as.list(methreads.file)
names(methreads.file) <- c("WT (Col-0)", "met1 (Col-0)", "met1 cmt3 (Col-0)", "nrpe1 (Col-0)", "ddm1 (Col-0)", "drm1/2 (Col-0)", "cmt2 (Col-0)")

top <- 1500
req_cols <- c(1,2,8)
i = 1
names(methreads.file)[i]
test.mut <- read.csv(file = file.path(methreads.dir, as.character(methreads.file)[i]), header = F)
colnames(test.mut) <- get("init_cls", envir = .GlobalEnv)
req.mat <- log10(na.omit(test.mut[,c(4,5,6,7,8,9,10,11)]))  ## Gene fractions filtering non methylated
req.mat[req.mat == -Inf] = 0
req.mat <- req.mat[1:top,]

cluster_wt <- hclust(dist(req.mat, method = "binary"))

plot_heatmap(methreads.dir = methreads.dir, methreads.file = methreads.file, i = i, top = top, cluster_wt = cluster_wt)

i = 6
names(methreads.file)[i]
plot_heatmap(methreads.dir = methreads.dir, methreads.file = methreads.file, i = i, top = top, cluster_wt = cluster_wt)


i = 2
names(methreads.file)[i]
plot_heatmap(methreads.dir = methreads.dir, top = top, methreads.file = methreads.file, i = i, cluster_wt = cluster_wt)


i = 7
names(methreads.file)[i]
plot_heatmap(methreads.dir = methreads.dir, top = top, methreads.file = methreads.file, i = i, cluster_wt = cluster_wt)

#### __________________________________________
#### 
araport.genes.stend <- paste(araport.genes$V1,araport.genes$V2,araport.genes$V3, sep = ",")
init_cls = c("chr", "start", "end", "(0,0,0)","(1,0,0)","(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")

mutant.data <- read.table("/projects/cegs/rahul/013.alignMutants_GSE39901/SraRunTable_StroudMutants.txt", sep = "\t", header = T, as.is = T)

methreads.fol <- "/lustre/scratch/projects/cegs/rahul/016.bshap/005.mutants.GSE39901/100.ARCHIVE.methReads_aragenes/"
methreads.files <- list.files(methreads.fol, pattern = "summary.txt")

required.frac.ind <- 11
init_cls[required.frac.ind]

required.frac.methreads <- c()
for(mut in methreads.files){
  t.mut <- read.csv(file = file.path(methreads.fol, mut), header = F)
  t.id <- mutant.data$genotype_s[which(mutant.data$Run_s == unlist(strsplit(mut, "[.]"))[2])]
  t.req.frac <- rep(NA, length(araport.genes.stend))
  req.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")
  ### taking the fraction of required column
  #req.sums <- rowSums(t.mut[,c(5,6,7,8,9,10,11)])
  #t.req.inds <- which(req.sums > 0)
  #t.req.frac[match(req.ids[t.req.inds], araport.genes.stend)] = t.mut[t.req.inds,required.frac.ind] / req.sums[t.req.inds]
  ### Now plotting the absolute numbers
  req.nums <- log10(t.mut[,required.frac.ind])
  req.nums[req.nums == -Inf] = 0
  t.req.frac[match(req.ids, araport.genes.stend)] = req.nums
  ###
  required.frac.methreads <- cbind(required.frac.methreads, t.req.frac)
  colnames(required.frac.methreads)[length(colnames(required.frac.methreads))] <- t.id
}

Heatmap(na.omit(required.frac.methreads), cluster_rows=T,cluster_columns =T, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_rows = "binary", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = paste("%", init_cls[required.frac.ind]), color_bar = "continuous"))



