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
araport.genes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_genes_201606.bed", as.is = T)


output_fol <- "/projects/cegs/rahul/016.bshap/005.mutants.GSE39901/001.methReads_aragenes/"
#output_fol <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0/001.methReads_aragenes/"
setwd(output_fol)

init_cls = c("chr", "start", "end", "(0,0,0)","(1,0,0)","(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")

getReqMat <- function(test.mut, case = 1){
  if (case == 1){
    ### Taking the fraction of each cluster with total reads mapped to the gene
    req.mat <- test.mut[,c(4,5,6,7,8,9,10,11)]/rowSums(test.mut[,c(4,5,6,7,8,9,10,11)])  ## Gene fractions filtering non methylated
    req.mat <- na.omit(req.mat)
    req.mat.description <- "Fraction of reads\nin each cluster"
  } else if (case == 2) {
    #### Just taking the number of reads as logarthimic number
    req.mat <- log10(na.omit(test.mut[,c(4,5,6,7,8,9,10,11)]))  ## Gene fractions filtering non methylated
    req.mat[req.mat == -Inf] = 0
    req.mat.description <- "Log10 \nnumber of reads"
  } else if (case == 3){
    ## Taking fracion now with the total reads and logarithmic
    req.mat <- log10(na.omit(test.mut[,c(4,5,6,7,8,9,10,11)]/rowSums(test.mut[,c(4,5,6,7,8,9,10,11)])))  ## Gene fractions filtering non methylated
    req.mat[req.mat == -Inf] = 0
    req.mat.description <- "Log10 \n fraction of reads in each cluster"
  } else if (case == 4) {
    ## Conditional probability on the number of reads matching
    req.mat <- na.omit((sum(rowSums(test.mut[,c(4,5,6,7,8,9,10,11)]))*test.mut[,c(5,6,7,8,9,10,11)])/rowSums(test.mut[,c(4,5,6,7,8,9,10,11)]))
    req.mat.description <- "Conditional probability on the number of reads matching"
  }
  return(list(req.mat, req.mat.description))
}


### Checking each single file......
plot_heatmap <- function(methreads.dir, methreads.file, i, cluster_wt=NULL, req_cols = NULL, top = NULL, case = 2){
  test.mut <- read.csv(file = file.path(methreads.dir, as.character(methreads.file)[i]), header = F)
  colnames(test.mut) <- get("init_cls", envir = .GlobalEnv)
  req.mat.list <- getReqMat(test.mut, case = case)
  req.mat <- req.mat.list[[1]]
  req.mat.description <- req.mat.list[[2]]
  
  if(!is.null(top)){
    req.mat <- req.mat[1:top,]
  } else{
    req.mat <- req.mat[,]
  }
  if(!is.null(req_cols)){
    req.mat <- req.mat[,req_cols]
  }
  
  if(!is.null(cluster_wt)){
    draw(Heatmap(req.mat, cluster_rows=cluster_wt,cluster_columns =F, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_rows = "binary", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = req.mat.description, color_bar = "continuous"), row_title = paste("Reads mapped to genes in", names(methreads.file)[i])))
  } else{
    draw(Heatmap(req.mat, cluster_rows=T,cluster_columns =F, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_rows = "binary", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = req.mat.description, color_bar = "continuous"), row_title = paste("Reads mapped to genes in", names(methreads.file)[i])))
  }
}

mutant.data$Run_s[which(mutant.data$genotype_s == "met1")]

methreads.dir <- "/projects/cegs/rahul/016.bshap/005.mutants.GSE39901/001.methReads_aragenes/"

methreads.file <- c("meths.SRR534177.summary.txt", "meths.SRR534239.summary.txt", "meths.SRR534240.summary.txt", "meths.SRR534182.summary.txt", "meths.SRR534215.summary.txt", "meths.SRR534223.summary.txt", "meths.SRR869314.summary.txt")
methreads.file <- as.list(methreads.file)
names(methreads.file) <- c("WT (Col-0)", "met1 (Col-0)", "met1 cmt3 (Col-0)", "nrpe1 (Col-0)", "ddm1 (Col-0)", "drm1/2 (Col-0)", "cmt2 (Col-0)")

top <- 1500
req_cols <- c(1,2,8)
i = 2 
names(methreads.file)[i]
test.mut <- read.csv(file = file.path(methreads.dir, as.character(methreads.file)[i]), header = F)
colnames(test.mut) <- get("init_cls", envir = .GlobalEnv)


req.mat.list <- getReqMat(test.mut, case = 1)
req.mat <- req.mat.list[[1]][1:top,]

hist(req.mat.list[[1]][,8][which(req.mat.list[[1]][,8] > 0.003)], breaks = 500)

which(req.mat.list[[1]][,8] > 0.003)

total.reads.genes <- rowSums(test.mut[,c(4,5,6,7,8,9,10,11)])
total.reads.genes <- total.reads.genes[which(total.reads.genes > 0)]


plot(total.reads.genes, req.mat.list[[1]][,8], pch = 19)



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

mutant.data <- read.table("/projects/cegs/rahul/013.alignMutants/01.MutantsStraud.GSE39901/SraRunTable_StroudMutants.txt", sep = "\t", header = T, as.is = T)
methreads.fol <- "/lustre/scratch/projects/cegs/rahul/016.bshap/005.mutants.GSE39901/001.methReads_aragenes/"
methreads.files <- list.files(methreads.fol, pattern = "summary.txt")

required.frac.ind <- 11
init_cls[required.frac.ind]

required.frac.methreads <- c()
for(mut in methreads.files){
  t.id <- mutant.data$genotype_s[which(mutant.data$Run_s == unlist(strsplit(mut, "[.]"))[2])]
  t.req.frac <- rep(NA, length(araport.genes.stend))
  try(t.mut <- read.csv(file = file.path(methreads.fol, mut), header = F), silent = T)
  if(inherits(t.mut, 'try-error')) {next}
  req.sums <- rowSums(t.mut[,c(5,6,7,8,9,10,11)])
  ### taking the fraction of required column
  sample.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")[which(req.sums > 0)]
  common.ids <- sample.ids[which(sample.ids %in% araport.genes.stend)]
  t.req.inds <- which(sample.ids %in% common.ids)
  t.req.frac[which(araport.genes.stend %in% common.ids)] = t.mut[t.req.inds,required.frac.ind] / req.sums[t.req.inds]
  ### Now plotting the absolute numbers
  #sample.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")
  #common.ids <- sample.ids[which(sample.ids %in% araport.genes.stend)]
  #req.nums <- log10(t.mut[,required.frac.ind])
  #req.nums[req.nums == -Inf] = 0
  #t.req.frac[which(araport.genes.stend %in% common.ids)] = req.nums
  ###
  required.frac.methreads <- cbind(required.frac.methreads, t.req.frac)
  colnames(required.frac.methreads)[length(colnames(required.frac.methreads))] <- t.id
}

Heatmap(na.omit(required.frac.methreads), cluster_rows=T,cluster_columns =T, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_columns = "canberra", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = paste("%", init_cls[required.frac.ind]), color_bar = "continuous"))

####
## Checking for the columbia accession in various runs
###
init_cls = c("chr", "start", "end", "(0,0,0)","(1,0,0)","(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
araport.genes.stend <- paste(araport.genes$V1,araport.genes$V2,araport.genes$V3, sep = ",")[1:2500]

accs.data <- read.table("/projects/cegs/rahul/008.Col-0.Bis/sra_info.txt", sep = "\t", header = T, as.is = T)
methreads.fol <- "/projects/cegs/rahul/016.bshap/001.col-0/001.methReads_aragenes/"
methreads.files <- list.files(methreads.fol, pattern = "summary.txt")[2:5]
methreads.files <- c("meths.SRR771698.summary.txt", 'meths.SRR534193.summary.txt')

required.frac.ind <- 11
init_cls[required.frac.ind]


required.frac.methreads <- c()
for(mut in methreads.files){
  t.req.frac <- rep(NA, length(araport.genes.stend))
  t.id <- accs.data$Run_s[which(accs.data$Run_s == unlist(strsplit(mut, "[.]"))[2])]
  try(t.mut <- read.csv(file = file.path(methreads.fol, mut), header = F), silent = T)
  if(inherits(t.mut, 'try-error')) {next}
  req.sums <- rowSums(t.mut[,c(5,6,7,8,9,10,11)])
  ### taking the fraction of required column
  sample.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")[which(req.sums > 0)]
  common.ids <- sample.ids[which(sample.ids %in% araport.genes.stend)]
  t.req.inds <- which(sample.ids %in% common.ids)
  t.req.frac[which(araport.genes.stend %in% common.ids)] = t.mut[t.req.inds,required.frac.ind] / req.sums[t.req.inds]
  ### Now plotting the absolute numbers
  #sample.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")
  #common.ids <- sample.ids[which(sample.ids %in% araport.genes.stend)]
  #req.nums <- log10(t.mut[,required.frac.ind])
  #req.nums[req.nums == -Inf] = 0
  #t.req.frac[which(araport.genes.stend %in% common.ids)] = req.nums
  ###
  required.frac.methreads <- cbind(required.frac.methreads, t.req.frac)
  colnames(required.frac.methreads)[length(colnames(required.frac.methreads))] <- t.id
}


head(required.frac.methreads)

Heatmap(na.omit(required.frac.methreads[1:10000,]), cluster_rows=T,cluster_columns =T, col=brewer.pal(11, "Spectral")[11:1], show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = paste("%", init_cls[required.frac.ind]), color_bar = "continuous"))

clustering_distance_rows = "binary"


####
## Checking the read bshap for different tissues in root meristem. 
###
init_cls = c("chr", "start", "end", "(0,0,0)","(1,0,0)","(0,1,0)", "(0,0,1)", "(1,1,0)", "(1,0,1)", "(0,1,1)", "(1,1,1)")
root.data <- read.table("/projects/cegs/rahul/017.RootMeristem.Taiji2016/SraRunTable_rootTaiji.txt", sep = "\t", header = T, as.is = T)
methreads.fol <- "/projects/cegs/rahul/016.bshap/004.taiji.rootmeristem/003.selectedDMRs/"
setwd(methreads.fol)

req_list=read.table("dmr_chh_rms_results_collapsed.filtered.bed", as.is = T)
methreads.files <- list.files(methreads.fol, pattern = "[.]dmr.summary.txt")

req_list = read.table("dmr_chh_rms_results_collapsed.random.bed", as.is = T)
methreads.files <- list.files(methreads.fol, pattern = "[.]nodmr.summary.txt")


req_list = paste(req_list$V1,req_list$V2,req_list$V3, sep = ",")
required.frac.ind <- 11
init_cls[required.frac.ind]

required.frac.methreads <- c()
for(mut in methreads.files){
  t.req.frac <- rep(NA, length(req_list))
  t.id <- root.data$source_name_s[which(root.data$Run_s == unlist(strsplit(mut, "[.]"))[2])]
  try(t.mut <- read.csv(file = file.path(methreads.fol, mut), header = F), silent = T)
  if(inherits(t.mut, 'try-error')) {next}
  req.sums <- rowSums(t.mut[,c(5,6,7,8,9,10,11)])
  ### taking the fraction of required column
  sample.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")[which(req.sums > 0)]
  common.ids <- sample.ids[which(sample.ids %in% req_list)]
  t.req.inds <- which(sample.ids %in% common.ids)
  t.req.frac[which(req_list %in% common.ids)] = t.mut[t.req.inds,required.frac.ind] / req.sums[t.req.inds]
  ### Now plotting the absolute numbers
  #sample.ids <- paste(t.mut[,c(1)],t.mut[,c(2)],t.mut[,c(3)], sep = ",")
  #common.ids <- sample.ids[which(sample.ids %in% req_list)]
  #req.nums <- log10(t.mut[,required.frac.ind])
  #req.nums[req.nums == -Inf] = 0
  #t.req.frac[which(req_list %in% common.ids)] = req.nums
  ###
  required.frac.methreads <- cbind(required.frac.methreads, t.req.frac)
  colnames(required.frac.methreads)[length(colnames(required.frac.methreads))] <- t.id
}

nrow(required.frac.methreads)
head(required.frac.methreads)

pdf("temp.nodmr.pdf")
Heatmap(na.omit(required.frac.methreads), cluster_rows=T,cluster_columns =T, col=brewer.pal(11, "Spectral")[11:1], show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = paste("%", init_cls[required.frac.ind]), color_bar = "continuous"), clustering_method_columns = "ward.D")
dev.off()











