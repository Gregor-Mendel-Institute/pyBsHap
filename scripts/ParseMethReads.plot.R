#### JSON parser for the meths on reads
library(jsonlite)
library(RColorBrewer)
library(rtracklayer)
resetPar <- function() {
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  op
}


colors.plot <- brewer.pal(5, "Paired")[c(1,2,1,2,1)]

json.file <- "/lustre/scratch/projects/cegs/rahul/006.SpermAndVegetativeCells/004.zilbermann/SRR516164/meths.SRR516164.bins.json"
json.file <- "/lustre/scratch/projects/cegs/rahul/006.SpermAndVegetativeCells/004.zilbermann/SRR516175/meths.SRR516175.json"
json.file <- "/lustre/scratch/projects/cegs/rahul/008.Col-0.Bis/02.methylpy/meths.SRR771698.json"
json.file <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/meths.SRR771698.CN.json"
json.file <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/temp.meths.CN.json"
json.file.cg <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/temp.meths.CG.json"
json.file.chg <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/temp.meths.CHG.json"
json.file.chh <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/temp.meths.CHH.json"

meths <- fromJSON(json.file)
meths.sperm <- fromJSON(json.file)
meths.cg <- fromJSON(json.file.cg)
meths.chg <- fromJSON(json.file.chg)
meths.chh <- fromJSON(json.file.chh)

json.file <- c("/lustre/scratch/projects/cegs/rahul/016.bshap/004.taiji.rootmeristem/", "SRR3311820")
meths.cn <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CN.json", sep = "")))
meths.cg <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CG.json", sep = "")))
meths.chg <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CHG.json", sep = "")))
meths.chh <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CHH.json", sep = "")))

## Histogram for methylation in reads
## 
meth.calls <- c(as.numeric(as.character(meths$Chr1)), as.numeric(as.character(meths$Chr2)), as.numeric(as.character(meths$Chr3)), as.numeric(as.character(meths$Chr4)), as.numeric(as.character(meths$Chr5)))
meth.calls <- meth.calls[!is.na(meth.calls)]
meth.calls <- meth.calls[which(meth.calls != -1)]
hist.meth <- hist(meth.calls, plot = F, breaks = 50)
hist.meth$density <- hist.meth$counts/sum(hist.meth$counts)
meth.density <- density(meth.calls, from = 0, to = 1)


cex.plot <- 1.2

#pdf("Meth.reads.pdf")
zones=matrix(c(1,1,2,3), ncol=2, byrow=T)
layout(zones, heights=c(1/5,4/5))
par(mar=c(6,5,0,1))
plot.new()
title(main = "Methylation on BS-seq reads for Col_0", line = -6, cex.main = cex.plot * 1.4)
title(main = "SRR771698", line = -8, cex.main = cex.plot)
plot(meth.density$x * 100, (meth.density$y/max(meth.density$y))*max(hist.meth$density), type = "l", col = colors.plot[2], lwd = 4, ylim = c(0, 0.27), ylab = "Density", xlab = "% methylation", cex.axis = cex.plot, cex.lab = cex.plot, main = "")
hist(meth.calls * 100, breaks = 80, xlab = "% methylation", cex.axis = cex.plot, cex.lab = cex.plot, ylim = c(0, 20000), col = colors.plot[1], main = "")
abline(v = 50)
#dev.off()

##################
## Plotting a manhattan plot

#araport.te.gene <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_TE_gene.201606.bed", as.is = T)
#araport.genes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_genes_201606.bed", as.is = T)
#araport.tes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_TEs_201606.bed", as.is = T)
araport.gff <- import.gff3("/vol/HOME/ARAPORT11/Araport11_GFF3_genic_regions.genes.TEs.201606.gff")

getMethWinds <- function(check.gr,meths.all, updown){
  num.rows <- 10
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  check.pos <- c(as.numeric(sub("Chr", "", as.character(seqnames(check.gr)), ignore.case = T)), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% meths.all$binlen), check.pos[3], meths.all$binlen))
  chrom.mat.cn = chrom.mat.cg = chrom.mat.chg = chrom.mat.chh = numeric()
  for (x in check.region){
    binmeths = meths.all[[meths.all$chrs[check.pos[1]]]][[x]]
    ## context inds cg = 2, chg = 3, chh = 4
    if (nrow(binmeths) > 0){
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

meth.region.plot <- function(check.gr,meths, updown = 2000){
  num.rows <- 10
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  chrom.mat = getMethWinds(check.gr,meths, updown)[[1]] ### Taking only the totals
  meth.colors <- c("grey", brewer.pal(num.rows, "Spectral")[(num.rows-1):1])
  par(resetPar())
  cex.plot <- 1.5
  barplot(chrom.mat, space = 0, col = meth.colors, border = F, las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3, xlab = paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "), sub = paste("input:", meths.all$input_bam))
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/meths$binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)) - as.integer(updown/meths$binlen), labels = "end", lwd.ticks = 5, line =0.5)
   par(fig=c(0.8,0.93,0.8,0.83), new=T)
   par(mar=c(0,0,0,0))
   ylim <- length(meth.breaks) - 1
   plot(0, type='n', ann=F, axes=F, xaxs="i", yaxs="i", ylim=c(0,1), xlim=c(0,ylim))
   axis(1, at=c(1.5, ylim), labels=c(0, 1), tick=FALSE, line=-1.2, las=1, cex.axis = cex.plot * 0.8)
   mtext(text="% methylation", las=1, side=3, line=0, outer=FALSE, cex=cex.plot)
   for(z in seq(ylim)){
     rect(z-1, 0, z, 1, col=meth.colors[z], border='black', lwd=0.5)
   }
  par(resetPar())
}

meth.region.plot.contexts <- function(check.gr,meths.all, updown = 2000){
  num.rows <- 10
  meth.breaks <- c(-1, seq(0, 1, length.out = num.rows))
  chrom.mat = getMethWinds(check.gr,meths.all, updown)
  meth.colors <- c("grey", brewer.pal(num.rows, "Spectral")[(num.rows-1):1])
  cex.plot <- 1.5
  zones=matrix(c(1,1,2,2,3,3,4), ncol=1, byrow=T)
  layout(zones)
  par(mar=c(1,4.5,1,2))
  barplot(chrom.mat[[2]], space = 0, col = meth.colors, border = F, las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3)
  mtext(text = "CG", side = 1, line = 1)
  par(mar=c(1,4.5,1,2))
  barplot(chrom.mat[[3]], space = 0, col = meth.colors, border = F, las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3)
  mtext(text = "CHG", side = 1, line = 1)
  par(mar=c(1,4.5,1,2))
  barplot(chrom.mat[[4]], space = 0, col = meth.colors, border = F, las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3)
  mtext(text = "CHH", side = 1, line = 1)
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/meths.all$binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat[[1]])), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat[[1]])) - as.integer(updown/meths.all$binlen), labels = "end", lwd.ticks = 5, line =0.5)
  plot.new()
  mtext(text=paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "), cex = cex.plot, line = -3)
  mtext(text=paste("input:", meths.all$input_bam), cex = 1, line = -5)
  legend("bottom", c("NA", 0,rep("", (num.rows - 3)), 1), fill=meth.colors, horiz = T, border = F, cex = cex.plot, bty = "n")
}

## Try checking this plots for different regions

json.file <- "/lustre/scratch/projects/cegs/rahul/016.bshap/004.taiji.rootmeristem/meths.SRR3311819.json"
meths.all <- fromJSON(json.file)

meth.region.plot(check.gr = GRanges("Chr1", IRanges(11705950,11724300)), meths = meths, updown = 2000)

ara.ind <- 161
ara.ind <- 106
ara.ind <- 14468
ara.ind <- 17145

meth.region.plot(check.gr = araport.gff[ara.ind],  meths.all, updown = 2000)
meth.region.plot.contexts(check.gr = araport.gff[ara.ind], meths.all = meths.all, updown = 2000)


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



json.file <- "/lustre/scratch/projects/cegs/rahul/016.bshap/004.taiji.rootmeristem/meths.SRR3311819.json"
meths.all <- fromJSON(json.file)

meth.region.plot(check.gr = GRanges("Chr1", IRanges(11705950,11724300)), meths = meths, updown = 2000)

ara.ind <- 161
ara.ind <- 106
ara.ind <- 14468  ### ROS1 gene
ara.ind <- 17145
check.gr <- araport.gff[ara.ind]

output_fol <- "~/Templates/"
ref_seq <- "/vol/HOME/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta"
input_folder <- "/lustre/scratch/projects/cegs/rahul/016.bshap/004.taiji.rootmeristem"
input_file <- "SRR3311820_processed_reads_no_clonal.bam"
output_id <- strsplit(input_file, "_")[[1]][1]
updown <- 5000
check_pos <- paste(as.character(seqnames(check.gr)), start(check.gr)-updown, end(check.gr)+updown, sep = ",")

##bshap getmeth -i SRR3311819_processed_reads_no_clonal.bam -r  -v -s Chr1,281350,295203  -o SRR3311819
pybshap.command <- paste("bshap getmeth -i", file.path(input_folder, input_file), "-r", ref_seq, "-v -o", output_id, "-s",  check_pos)
setwd(output_fol)
system(pybshap.command)
meths.all <- fromJSON(paste("meths.", output_id,  ".json", sep = ""))

meth.region.plot(check.gr = check.gr,  meths.all, updown = 5000)
meth.region.plot.contexts(check.gr = araport.gff[ara.ind], meths.all = meths.all, updown = 5000)

par(resetPar())
getMethCorrWind(check.gr, meths.all)


check.region <- as.character(seq(start(check.gr) + 1 - (start(check.gr) %% meths.all$binlen), end(check.gr), meths.all$binlen))
meths.all[[as.character(seqnames(check.gr))]][[check.region[6]]]



