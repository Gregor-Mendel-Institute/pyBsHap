#### JSON parser for the meths on reads
library(jsonlite)
library(RColorBrewer)
library(rtracklayer)
library(data.table)
resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}
context.colors <- brewer.pal(3, "Set2")

## Histogram for methylation in reads
## 
getMethDensityPlot <- function(meths){
  meth.calls <- c(as.numeric(as.character(unlist(meths[['Chr1']]))), as.numeric(as.character(unlist(meths[['Chr2']]))), as.numeric(as.character(unlist(meths[['Chr3']]))), as.numeric(as.character(unlist(meths[['Chr4']]))), as.numeric(as.character(unlist(meths[['Chr5']]))))
  meth.calls <- meth.calls[!is.na(meth.calls)]
  meth.calls <- meth.calls[which(meth.calls != -1)]
  hist.meth <- hist(meth.calls, plot = F, breaks = 50)
  hist.meth$density <- hist.meth$counts/sum(hist.meth$counts)
  meth.density <- density(meth.calls, from = 0, to = 1)
  zones=matrix(c(1,1,2,3), ncol=2, byrow=T)
  colors.plot <- brewer.pal(5, "Paired")[c(1,2,1,2,1)]
  layout(zones, heights=c(1/5,4/5))
  par(mar=c(6,5,0,1))
  plot.new()
  title(main = "Methylation on BS-seq reads for Col_0", line = -6, cex.main = cex.plot * 1.4)
  title(main = "SRR771698", line = -8, cex.main = cex.plot)
  plot(meth.density$x * 100, (meth.density$y/max(meth.density$y))*max(hist.meth$density), type = "l", col = colors.plot[2], lwd = 4, ylim = c(0, 0.27), ylab = "Density", xlab = "% methylation", cex.axis = cex.plot, cex.lab = cex.plot, main = "")
  hist(meth.calls[which(meth.calls != 0)] * 100, breaks = 80, xlab = "% methylation", cex.axis = cex.plot, cex.lab = cex.plot, xlim = c(0, 100), col = colors.plot[1], main = "")
  abline(v = 50)
}



cex.plot <- 1.2
getMethDensityPlot(meths)
#pdf("Meth.reads.pdf")
#dev.off()

##################
## Plotting a manhattan plot

#araport.te.gene <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_TE_gene.201606.bed", as.is = T)
#araport.genes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_genes_201606.bed", as.is = T)
#araport.tes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_TEs_201606.bed", as.is = T)
araport.gff <- import("/vol/HOME/ARAPORT11/Araport11_GFF3_genic_regions.genes.TEs.201606.gff")

getMethylationAverages <- function(binmeths, ind = 1){
  finAvg = NA
  if (ind == 1){
    finAvg <- sum(binmeths$mC) / sum(binmeths$tC)
  } else if (ind == 2){ ## Get the number of methylated positions
    finAvg <- mean(binmeths$call)
  } else if (ind == 3){  ### Get the normal average
    finAvg <- mean(mapply(FUN = function(x, y){x / y}, binmeths$mC, binmeths$tC))
  }
  return(finAvg)
}

getWMA_default_plot <- function(check.gr, input_h5file, meths_fol, input_id, updown, title = "", getMethInd = 1, binlen = 300){
  yaxis.labs <- c("Weighted methylation average", "Average methylated positions", "Average methylation levels")
  #binlen = h5read(input_h5file, "binlen")
  meths.chrs <- h5read(input_h5file, "chrs")
  check.pos <- c(which(meths.chrs == as.character(seqnames(check.gr))), start(check.gr) - updown, end(check.gr) + updown)
  #check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% binlen), check.pos[3], binlen))
  check.region <- seq(check.pos[2] + 1 - (check.pos[2] %% binlen), check.pos[3], binlen)
  meth.files <- list.files(meths_fol, pattern = paste("^allc_", input_id, "_.*[12345].tsv", sep = ""))
  bs.mC <- fread(file.path(meths_fol, meth.files[check.pos[1]]), showProgress = T)
  colnames(bs.mC) <- c("chr", "pos", "strand", "context", "mC", "tC", "call")
  bs.mC <- subset(bs.mC, pos >= check.pos[2] - 100 & pos <= check.pos[3] + 100)
  weightedAvg.meths.cn <- numeric()
  weightedAvg.meths.cg <- numeric()
  weightedAvg.meths.chg <- numeric()
  weightedAvg.meths.chh <- numeric()
  Avg.cytosines.meths <- numeric()
  for (i in check.region){
    binmeths <- subset(bs.mC, pos >= i & pos < i + binlen)
    cg.ind <- which(grepl("CG", binmeths$context, perl = T))
    chg.ind <- which(grepl("C[ATC]G", binmeths$context, perl = T))
    chh.ind <- which(grepl("C[ATC][ATC]", binmeths$context, perl = T))
    Avg.cytosines.meths <- rbind(Avg.cytosines.meths, cbind(binmeths$pos, mapply(FUN = function(x, y){x / y}, binmeths$mC, binmeths$tC), 1))
    Avg.cytosines.meths[chg.ind, 3] <- 2
    Avg.cytosines.meths[chh.ind, 3] <- 3
    weightedAvg.meths.cg <- c(weightedAvg.meths.cg, getMethylationAverages(binmeths[cg.ind, ], ind = getMethInd))
    weightedAvg.meths.chg <- c(weightedAvg.meths.chg, getMethylationAverages(binmeths[chg.ind,], ind = getMethInd))
    weightedAvg.meths.chh <- c(weightedAvg.meths.chh, getMethylationAverages(binmeths[chh.ind,], ind = getMethInd))
    weightedAvg.meths.cn <- c(weightedAvg.meths.cn, getMethylationAverages(binmeths, ind = getMethInd))
  }
  plot(x = check.region, y = weightedAvg.meths.cn, ylim = c(0, 1), ylab = yaxis.labs[getMethInd], cex.lab = cex.plot, xaxt = "n", type = "l", lwd = 2, xlab = paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "), sub = paste(title, ",", input_h5file))
  points(x = check.region, y = weightedAvg.meths.cg, type = "l", lwd = 2, col = context.colors[1])
  points(x = check.region, y = weightedAvg.meths.chg, type = "l", lwd = 2, col = context.colors[2])
  points(x = check.region, y = weightedAvg.meths.chh, type = "l", lwd = 2, col = context.colors[3])
  #points(x = Avg.cytosines.meths[,1], Avg.cytosines.meths[,2], col = context.colors[Avg.cytosines.meths[,3]], pch = 19)
  legend("topright", legend = c("CG", "CHG", "CHH"), fill = context.colors)
  axis(1, at = check.region[1], labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = start(check.gr), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = check.region[length(check.region)], labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = end(check.gr), labels = "end", lwd.ticks = 5, line =0.5)
}

get_meth_colors = function(num.rows){
  if (num.rows == 10){
    #meth.colors <- c("grey", brewer.pal(num.rows, "Spectral")[num.rows:2])
    #return(c("grey", brewer.pal(num.rows-1, "Blues")))
    return(c("grey", brewer.pal(num.rows-1, "GnBu")))
  }
}

getMethWinds <- function(check.gr,input_h5file, updown){
  num.rows <- 10
  binlen = h5read(input_h5file, "binlen")
  meths.chrs <- h5read(input_h5file, "chrs")
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  check.pos <- c(which(meths.chrs == as.character(seqnames(check.gr))), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% binlen), check.pos[3], binlen))
  chrom.mat.cn = chrom.mat.cg = chrom.mat.chg = chrom.mat.chh = numeric()
  input_h5file_groups <- h5ls(input_h5file)
  for (x in check.region){
    H5close()
    bin_name = paste("b", meths.chrs[check.pos[1]], x, sep = "_")
    if (bin_name %in% input_h5file_groups$name){
      binmeths = h5read(input_h5file, bin_name)
      ## context inds cg = 2, chg = 3, chh = 4
      cn.meth = as.numeric(table(cut(binmeths[,1], breaks = meth.breaks, include.lowest = T, right = F)))
      cg.meth = as.numeric(table(cut(binmeths[,2], breaks = meth.breaks, include.lowest = T, right = F)))
      chg.meth = as.numeric(table(cut(binmeths[,3], breaks = meth.breaks, include.lowest = T, right = F)))
      chh.meth = as.numeric(table(cut(binmeths[,4], breaks = meth.breaks, include.lowest = T, right = F)))
    } else {
      cn.meth = cg.meth = chg.meth = chh.meth = rep(0,num.rows)
    }
    chrom.mat.cn <- cbind(chrom.mat.cn, cn.meth)
    chrom.mat.cg <- cbind(chrom.mat.cg, cg.meth)
    chrom.mat.chg <- cbind(chrom.mat.chg, chg.meth)
    chrom.mat.chh <- cbind(chrom.mat.chh, chh.meth)
  }
  meths.names <- levels(cut(c(0,0), breaks=round(meth.breaks, 2), include.lowest = T, right = F))
  rownames(chrom.mat.cn) = rownames(chrom.mat.cg) = rownames(chrom.mat.chg) = rownames(chrom.mat.chh) = meths.names
  colnames(chrom.mat.cn) = colnames(chrom.mat.cg) = colnames(chrom.mat.chg) = colnames(chrom.mat.chh) = check.region
  return(list(chrom.mat.cn, chrom.mat.cg, chrom.mat.chg, chrom.mat.chh))
}

meth.region.plot <- function(check.gr,input_h5file, title="", updown = 2000){
  num.rows <- 10
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  chrom.mat = getMethWinds(check.gr,input_h5file, updown)[[1]] ### Taking only the totals
  meth.colors <- get_meth_colors(num.rows = num.rows)
  input_bam <- h5read(input_h5file, "input_bam")
  binlen = h5read(input_h5file, "binlen")
  cex.plot <- 1.5
  par(mar = c(5.5, 4.5, 3, 4) )
  barplot(chrom.mat, space = 0, col = meth.colors, border = "#bdbdbd", las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3, xlab = paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "), sub = paste(title, ",", input_h5file))
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)) - as.integer(updown/binlen), labels = "end", lwd.ticks = 5, line =0.5)
  par(fig=c(0.8,0.93,0.8,0.83), new=T)
  par(mar=c(0,0,0,0))
  ylim <- length(meth.breaks) - 1
  plot(0, type='n', ann=F, axes=F, xaxs="i", yaxs="i", ylim=c(0,1), xlim=c(0,ylim))
  axis(1, at=c(1.5, ylim), labels=c(0, 1), tick=FALSE, line=-1.2, las=1, cex.axis = cex.plot * 0.8)
  mtext(text="% methylation", las=1, side=3, line=0, outer=FALSE, cex=cex.plot)
  for(z in seq(ylim)){
    rect(z-1, 0, z, 1, col=meth.colors[z], border='black', lwd=0.5)
  }
}

drawMethPlot <- function(check.gr, input_h5file, context, max_reads, updown, cex.plot = 1.5){
  num.rows <- 10
  binlen = h5read(input_h5file, "binlen")
  meths.chrs <- h5read(input_h5file, "chrs")
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  check.pos <- c(which(meths.chrs == as.character(seqnames(check.gr))), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% binlen), check.pos[3], binlen))
  meth.colors <- get_meth_colors(num.rows = num.rows)
  plot(-10, -10, xlim = c(0, length(check.region)), ylim = c(0, max_reads + 30), ylab = "Number of reads", xaxt = "n", xlab = "", cex.lab = cex.plot, frame.plot=FALSE)
  input_h5file_groups <- h5ls(input_h5file)
  for (xind in seq(length(check.region))){
    bin_name = paste("b", meths.chrs[check.pos[1]], check.region[xind], sep = "_")
    if (bin_name %in% input_h5file_groups$name){
      binmeths = h5read(input_h5file, bin_name)
      binmeths <- binmeths[order(binmeths[,2], binmeths[,3], binmeths[,4]),]
      if (xind %% 10 == 0){
        abline(v = xind, col = "gray60", lty = 2)
        }
      for (j in seq(nrow(binmeths))){
        methcol = meth.colors[as.numeric(cut(binmeths[j-1,context], breaks = meth.breaks, include.lowest = T, right = F))]
        rect(xleft = xind, xright = xind + 1, ybottom = j, ytop = j + 1, col = methcol, border = F)
      }
    }
  }
}

meth.region.plot.contexts <- function(check.gr,input_h5file, mtitle = "", updown = 2000,max_reads = 40){
  binlen = h5read(input_h5file, "binlen")
  meths.chrs <- h5read(input_h5file, "chrs")
  check.pos <- c(which(meths.chrs == as.character(seqnames(check.gr))), start(check.gr) - updown, end(check.gr) + updown)
  input_bam <- h5read(input_h5file, "input_bam")
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% binlen), check.pos[3], binlen))
  num.rows = 10
  cex.plot = 1.5
  meth.colors <- get_meth_colors(num.rows = num.rows)
  zones=matrix(c(1,1,2,2,3,3,4), ncol=1, byrow=T)
  layout(zones)
  par(mar=c(1,4.5,1,2))
  drawMethPlot(check.gr = check.gr, input_h5file, context = 2, max_reads = max_reads, updown = updown)
  mtext(text = "CG", side = 1, line = 1)
  par(mar=c(1,4.5,1,2))
  drawMethPlot(check.gr = check.gr, input_h5file, context = 3, max_reads = max_reads, updown = updown)
  mtext(text = "CHG", side = 1, line = 1)
  par(mar=c(1,4.5,1,2))
  drawMethPlot(check.gr = check.gr, input_h5file, context = 4, max_reads = max_reads, updown = updown)
  mtext(text = "CHH", side = 1, line = 1)
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(check.region), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(check.region) - as.integer(updown/binlen), labels = "end", lwd.ticks = 5, line =0.5)
  plot.new()
  mtext(text=paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "), cex = cex.plot, line = -3)
  mtext(text= paste(mtitle, ",", input_h5file), cex = 1, line = -5)
  legend("bottom", c("NA", 0,rep("", (10 - 3)), 1), fill=meth.colors, horiz = T, border = F, cex = cex.plot, bty = "n")
}

getPCAmeths <- function(check.gr, input_h5file,main,updown=5000,k=3){
  cex.plot <- 1.2
  input_bam <- h5read(input_h5file, "input_bam")
  binlen = h5read(input_h5file, "binlen")
  meths.chrs <- h5read(input_h5file, "chrs")
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  check.pos <- c(which(meths.chrs == as.character(seqnames(check.gr))), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% binlen), check.pos[3], binlen))
  input_h5file_groups <- h5ls(input_h5file)
  meths.contexts <- data.frame()
  for(x in check.region){
    H5close()
    bin_name = paste("b", meths.chrs[check.pos[1]], x, sep = "_")
    if (bin_name %in% input_h5file_groups$name){
      binmeths = h5read(input_h5file, bin_name)
      binmeths[binmeths == -1] <- NA
      meths.contexts <- rbind(meths.contexts, as.data.frame(na.omit(binmeths)))
    }
  }
  colnames(meths.contexts) = c("CN","CG","CHG","CHH")
  meths.contexts.filtered <- meths.contexts[which(meths.contexts$CN>0), c(2,3,4)]
  #var(meths.contexts.filtered,na.rm=T,use="pairwise.complete.obs")
  #pcont <- princomp(na.omit(meths.contexts.filtered), cor = TRUE, scores = TRUE)
  #biplot(pcont)
  #library("rgl")
  #plot3d(pcont$scores[,1:3])
  #plot(meths.contexts.filtered$CG, meths.contexts.filtered$CHG)
  #
  x_label <- paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", ")
  file_label <- paste(main, ",", input_h5file)
  draw(Heatmap(meths.contexts.filtered, row_title = file_label, column_title = x_label, cluster_rows=T,cluster_columns =F, col=brewer.pal(11, "Spectral")[11:1], clustering_distance_rows = "binary", show_row_names = FALSE, show_heatmap_legend=T, heatmap_legend_param = list(title = "% methylation", color_bar = "continuous")))
  #di <- dist(meths.contexts.filtered,method="binary")
  #tree <- hclust(di, method="average")
  #plot(tree, main="Reads clustered based on their methylation states",xlab=x_label, cex.lab = cex.plot, sub = file_label)
  #rect.hclust(tree, k = k, border="red")
}

## Try checking this plots for different regions
#### h5py FILES
library(rhdf5)
library(RColorBrewer)
library("ComplexHeatmap")
library("reshape")
ref_seq <- "/vol/HOME/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta"
output_fol <- "~/Templates/"
setwd(output_fol)
Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "/home/GMI/rahul.pisupati/py2_kernel/bin:/home/GMI/rahul.pisupati/anaconda2/bin", sep = ":"))

## Random files

sra_names = c("WT col-0")
sra_files = c("/projects/cegs/rahul/008.Col-0.Bis/03.sra.SRR501624/SRR501624_processed_reads_no_clonal.bam")

i = 1

#### Straud mutants, bam files
data_dir = "/projects/cegs/rahul/013.alignMutants/01.MutantsStraud.GSE39901/002.methylpy"
sra_table = read.csv("/projects/cegs/rahul/013.alignMutants/01.MutantsStraud.GSE39901/SraRunTable_StroudMutants.txt", as.is = T, sep = "\t")
sra_names = sra_table[,c('genotype_s')] 
sra_files = as.character()
h5files = as.character()
for (ef in sra_table[,c('Run_s')]) { sra_files = c(sra_files, file.path(data_dir, ef, paste(ef, "_processed_reads_no_clonal.bam", sep = "")))  }

fix(sra_table)
i = 1
i = 59 ## Met1
i = 30 ## drm1/2
i = 32 ## drm1/2 and cmt3
i = 45 ## CMT2
i = 50 ## nrpe1
i = 17 ## CMT3
i = 89 ## vim123

#### Root tissues, bam files
data_dir = "/projects/cegs/rahul/017.RootMeristem.Taiji2016/002.methylpy/"
sra_table = read.csv("/projects/cegs/rahul/017.RootMeristem.Taiji2016/SraRunTable_rootTaiji.txt", as.is = T, sep = "\t")
sra_names = sra_table[,c('source_name_s')] 
sra_files = as.character()
for (ef in sra_table[,c('Run_s')]) { sra_files = c(sra_files, paste(data_dir, "/", ef, "/", ef, "_processed_reads_no_clonal.bam", sep = "") )  }

i = 8
i = 6 ## CRC
i = 1 # Epidermis
i = 5 # Stele

### ___________________________________________
### Choose a region in the genome
ara.ind <- 161
ara.ind <- 106
ara.ind <- 107
ara.ind <- 14468  ### ROS1 gene
ara.ind <- 17145  ## DML-2 protein gene
ara.ind <- which(elementMetadata(araport.gff)$ID == "AT1G65960")
ara.ind <- which(as.character(seqnames(araport.gff)) == "ChrC")[4]

check.gr <- araport.gff[ara.ind]
#check.gr <- subset(araport.gff, ID == "AT4G37650")

check.gr <- GRanges(seqnames = c("ChrC"), ranges = IRanges(start = 15000, end = 20000), Name = "ChrC", type = "chroloplast genome:5kb")   ### Entire ChrC

## Checking DMRs between the root tissues.
check.gr <- GRanges(seqnames = c("Chr1"), ranges = IRanges(start = 423086, end = 424360), Name = "DMR1", type = "AT1TE01370")   ### TE1 in DMR
## DMR for sperm and vegetative cell
check.gr <- GRanges(seqnames = c("Chr1"), ranges = IRanges(start = 431466	, end = 431983), Name = "DMR1", type = "AT1TE01395")   ### TE1 in DMR
## RdDM TE
check.gr <- GRanges(seqnames = c("Chr1"), ranges = IRanges(start = 2473099	, end = 2473306), Name = "DMR1", type = "AT1TE08055")   ### TE1 in DMR
### CMT2 TE
check.gr <- GRanges(seqnames = c("Chr1"), ranges = IRanges(start = 9642845, end = 9643142), Name = "DMR1", type = "AT1TE31080")   ### CMT2 TE
check.gr <- GRanges(seqnames = c("Chr1"), ranges = IRanges(start = 13324758, end = 13329639), Name = "DMR1", type = "AT1TE43585")   ### CMT2 TE
### ___________________________________________
## Run all the commands below
sra_names[i]
input_file <- sra_files[i]

output_id <- strsplit(basename(input_file), "_")[[1]][1]
updown <- 3000
check_pos <- paste(as.character(seqnames(check.gr)), start(check.gr)-updown, end(check.gr)+updown, sep = ",")

pybshap.command <- paste("bshap getmeth -i", input_file, "-r", ref_seq, "-v -o", output_id, "-s",  check_pos)

system(pybshap.command)

input_h5file <- paste("meths.", output_id,  ".hdf5", sep = "")

#meth.region.plot(check.gr,input_h5file, updown = updown, title = sra_names[i])
#meth.region.plot.contexts(check.gr = check.gr, input_h5file = input_h5file, mtitle = sra_names[i], updown = updown, max_reads = 40)
#dev.off()

#getPCAmeths(check.gr, input_h5file, main=sra_names[i], updown = updown)

#getWMA_default_plot(check.gr = check.gr, input_h5file = input_h5file, updown = updown, meths_fol = dirname(as.character(sra_files[i])), input_id = output_id, title = sra_names[i], getMethInd = 1, binlen = 500)



mhl.command <- paste("bshap getmhl -i", input_file, "-r", ref_seq, "-x",  check_pos)
mhl_regions = system(mhl.command, intern = T)

#dev.off()
meth.region.plot(check.gr,input_h5file, updown = updown, title = sra_names[i])
par(resetPar())
par(mar = c(5.5, 4.5, 10, 4) )
par(new=T)
plot(seq( length(mhl_regions) -1), tail(as.numeric(sapply(mhl_regions, function(x){ unlist(strsplit(x, ","))[4]   } )), -1) , axes = F, ylim = c(0, 1), xlab = "", ylab = "", pch = 19)
axis(4)
mtext("MHL", line = 2, 4)



### looping
#
n = 12
i = 4
updown <- 2000

zones=matrix(seq(n), ncol = 3, byrow = T)
layout(zones)

input_file <- bs.bams[[i]]
names(bs.bams)[i]
output_id <- strsplit(basename(input_file), "_")[[1]][1]
#check.list <- sample(seq(length(araport.gff)), n)
pdf("plot_met1.pdf")

for (ara.ind in check.list){
  check.gr <- araport.gff[ara.ind]
  check_pos <- paste(as.character(seqnames(check.gr)), start(check.gr)-updown, end(check.gr)+updown, sep = ",")
  pybshap.command <- paste("bshap getmeth -i", input_file, "-r", ref_seq, "-v -o", output_id, "-s",  check_pos)
  system(pybshap.command)
  input_h5file <- paste("meths.", output_id,  ".hdf5", sep = "")
  getPCAmeths(check.gr, input_h5file,main=names(bs.bams)[i])
  #meth.region.plot(check.gr,input_h5file, title=names(bs.bams)[i], updown = 2000)
}

dev.off()

##___________
round_int <- function(x, by=3){by*(as.integer(x/by)+as.logical(x%%by))}
updown <- 2000

ara.ind = 32131

zones=matrix(seq(round_int(length(bs.bams))), ncol = 3, byrow = T)
layout(zones)

check.gr <- araport.gff[ara.ind]
check_pos <- paste(as.character(seqnames(check.gr)), start(check.gr)-updown, end(check.gr)+updown, sep = ",")

pdf("plot_gene.pdf")
for (i in seq(length(bs.bams))){
  input_file <- bs.bams[[i]]
  output_id <- strsplit(basename(input_file), "_")[[1]][1]
  pybshap.command <- paste("bshap getmeth -i", input_file, "-r", ref_seq, "-v -o", output_id, "-s",  check_pos)
  system(pybshap.command)
  input_h5file <- paste("meths.", output_id,  ".hdf5", sep = "")
  getPCAmeths(check.gr, input_h5file,main=names(bs.bams)[i])
  #meth.region.plot(check.gr,input_h5file, title=names(bs.bams)[i], updown = 2000)
}
dev.off()


##### ============================================
##
##  Now trying to do the correlations in each window between the contexts
##
##
##### ============================================

pval_thres = 0.05
getEstimateCorTest <- function(conx_cony_corr){
  t.out = NA
  if (!is.na(conx_cony_corr$p.value)){
    if(conx_cony_corr$p.value < pval_thres){
      t.out = conx_cony_corr$estimate
    }
  } 
  return(t.out)
}

getMethCorrWind <- function(check.gr, meths.all, updown = 2000){
  check.pos <- c(as.numeric(sub("Chr", "", as.character(seqnames(check.gr)), ignore.case = T)), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% meths.all$binlen), check.pos[3], meths.all$binlen))
  chrom.corr.cg.chg = chrom.corr.cg.chh = chrom.corr.chg.chh =  numeric()
  for (x in check.region){
    binmeths = meths.all[[meths.all$chrs[check.pos[1]]]][[x]]
    binmeths[binmeths == -1] <- NA
    ## context inds cg = 2, chg = 3, chh = 4
    t.cg.chg <- try(cor.test(binmeths[,2], binmeths[,3]), silent = T)
    t.cg.chh <- try(cor.test(binmeths[,2], binmeths[,4]), silent = T)
    t.chg.chh <- try(cor.test(binmeths[,3], binmeths[,4]), silent = T)
    t.e.cg.chg = t.e.cg.chh = t.e.chg.chh = NA
    if(!inherits(t.cg.chg, 'try-error')) {t.e.cg.chg = getEstimateCorTest(t.cg.chg)}
    if(!inherits(t.cg.chh, 'try-error')) {t.e.cg.chh = getEstimateCorTest(t.cg.chh)}
    if(!inherits(t.chg.chh, 'try-error')) {t.e.chg.chh = getEstimateCorTest(t.chg.chh)}
    chrom.corr.cg.chg <- c(chrom.corr.cg.chg, as.numeric(t.e.cg.chg))
    chrom.corr.cg.chh <- c(chrom.corr.cg.chh, as.numeric(t.e.cg.chh))
    chrom.corr.chg.chh <- c(chrom.corr.chg.chh, as.numeric(t.e.chg.chh))
  }
  meth.colors <- brewer.pal(4, "Dark2")
  #par(resetPar())
  cex.plot <- 1.5
  plot(seq(length(check.region)), chrom.corr.cg.chg, col = meth.colors[1], pch = 19, cex = cex.plot, ylim = c(-1, 1), ylab = "Correlation with different contexts", xaxt = "n", cex.lab = cex.plot, xlab = paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "), sub = paste("input:", meths.all$input_bam))
  points(seq(length(check.region)), chrom.corr.cg.chh, col = meth.colors[2], pch = 19, cex = cex.plot)
  points(seq(length(check.region)), chrom.corr.chg.chh, col = meth.colors[3], pch = 19, cex = cex.plot)
  legend("topright", fill = meth.colors, c("CG-CHG", "CG-CHH", "CHG-CHH"))
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/meths.all$binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(check.region), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(check.region) - as.integer(updown/meths.all$binlen), labels = "end", lwd.ticks = 5, line =0.5)
}


meth.region.plot(check.gr = check.gr,  meths.all, updown = 5000)
meth.region.plot.contexts(check.gr = check.gr, meths.all = meths.all, updown = 3000)

par(resetPar())
getMethCorrWind(check.gr, meths.all)


check.region <- as.character(seq(start(check.gr) + 1 - (start(check.gr) %% meths.all$binlen), end(check.gr), meths.all$binlen))
meths.all[[as.character(seqnames(check.gr))]][[check.region[6]]]

