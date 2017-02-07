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
json.file.cg <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/meths.SRR771698.CG.json"
json.file.chg <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/meths.SRR771698.CHG.json"
json.file.chh <- "/lustre/scratch/projects/cegs/rahul/016.bshap/001.col-0.SRR771698/meths.SRR771698.CHH.json"

meths <- fromJSON(json.file)
meths.sperm <- fromJSON(json.file)
meths.cg <- fromJSON(json.file.cg)
meths.chg <- fromJSON(json.file.chg)
meths.chh <- fromJSON(json.file.chh)

json.file <- c("/lustre/scratch/projects/cegs/rahul/016.bshap/004.taiji.rootmeristem/methReads/", "SRR3311826")
meths.cn <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CN.json", sep = "")))
meths.cg <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CG.json", sep = "")))
meths.chg <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CHG.json", sep = "")))
meths.chh <- fromJSON(file.path(json.file[1], paste("meths.", json.file[2], ".CHH.json", sep = "")))




## Histogram for methylation in reads
## 
meth.calls <- c(as.numeric(as.character(meths$Chr1)), as.numeric(as.character(meths$Chr2)), as.numeric(as.character(meths$Chr3)), as.numeric(as.character(meths$Chr4)), as.numeric(as.character(meths$Chr5)))
meth.calls <- meth.calls[!is.na(meth.calls)]
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


# cex.plot = 1.2
# plot(0,0, xlim = c(0, sum(meths$chrslen)), ylim = c(0,100), xaxt = "n", ylab = "% methylation", xlab = "Chromosomes", type = "n", cex.lab = cex.plot, xaxs="i")
# chr.index = 1
# for (i in 1:length(meths$chrslen)){
#   axis(1, at = (chr.index + (meths$chrslen[i]/2)), labels = meths$chrs[i], cex.axis = cex.plot, tick = F)
#   chr.index = chr.index + meths$chrslen[i]
#   abline(v = chr.index, col = "grey")
# }


### Checking specific regions
check.region <- names(meths[[meths$chrs[i]]])[1:100]
#check.region <- as.character(seq(1, meths$chrslen[i], meths$binlen))
check.pos <- c(1, 54676, 57576) ## TE1, 1kb up and down
check.pos <- c(1, 33365,37871) ### Gene AT1G01060
check.pos <- c(1, 33365,37871) ## ge
check.pos <- c(1, 513699,515422) ## gene  AT1G02475
check.pos <- c(1,7064494,7075561) ### TEgene, 1kb up and down
check.pos <- c(1,3779235,3787293) ## TE gene, 1kb up and down, AT1G11265

#araport.te.gene <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_TE_gene.201606.bed", as.is = T)
#araport.genes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_genes_201606.bed", as.is = T)
#araport.tes <- read.table("/vol/HOME/ARAPORT11/Araport11_GFF3_TEs_201606.bed", as.is = T)
araport.gff <- import.gff3("/vol/HOME/ARAPORT11/Araport11_GFF3_genic_regions.genes.TEs.201606.gff")


meth.region.plot <- function(check.gr,meths, updown = 1000){
  num.rows <- 9
  meth.breaks <- seq(0, 1, length.out = num.rows)
  check.pos <- c(as.numeric(sub("Chr", "", as.character(seqnames(check.gr)), ignore.case = T)), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% meths$binlen), check.pos[3], meths$binlen))
  chrom.mat <- sapply(check.region, function(x){ if(length(meths[[meths$chrs[check.pos[1]]]][[x]]) > 0) {y = cut(meths[[meths$chrs[check.pos[1]]]][[x]], breaks = meth.breaks, include.lowest = T); as.numeric(table(y))}else {rep(0,length(meth.breaks)-1)}})
  #chrom.mat <- chrom.mat[,order(as.numeric(colnames(chrom.mat)))]
  rownames(chrom.mat) = levels(cut(c(0,0), breaks=round(meth.breaks, 2), include.lowest = T))
  meth.colors <-brewer.pal(num.rows, "Spectral")[num.rows:1]
  par(resetPar())
  cex.plot <- 1.5
  barplot(chrom.mat, space = 0, col = meth.colors, border = F, las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3, xlab = paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "))
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/meths$binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)) - as.integer(updown/meths$binlen), labels = "end", lwd.ticks = 5, line =0.5)
   par(fig=c(0.8,0.93,0.8,0.83), new=T)
   par(mar=c(0,0,0,0))
   ylim <- length(meth.breaks) - 1
   plot(0, type='n', ann=F, axes=F, xaxs="i", yaxs="i", ylim=c(0,1), xlim=c(0,ylim))
   axis(1, at=c(0, ylim), labels=c(0, 1), tick=FALSE, line=-1.2, las=1, cex.axis = cex.plot * 0.8)
   mtext(text="% methylation", las=1, side=3, line=0, outer=FALSE, cex=cex.plot)
   for(z in seq(ylim)){
     rect(z-1, 0, z, 1, col=meth.colors[z], border='black', lwd=0.5)
   }
  par(resetPar())
}

meth.region.plot.nolegend <- function(check.gr,meths, updown = 1000){
  meth.breaks <- seq(0, 1, length.out = 10)
  check.pos <- c(as.numeric(sub("Chr", "", as.character(seqnames(check.gr)), ignore.case = T)), start(check.gr) - updown, end(check.gr) + updown)
  check.region <- as.character(seq(check.pos[2] + 1 - (check.pos[2] %% meths$binlen), check.pos[3], meths$binlen))
  chrom.mat <- sapply(check.region, function(x){ if(length(meths[[meths$chrs[check.pos[1]]]][[x]]) > 0) {y = cut(meths[[meths$chrs[check.pos[1]]]][[x]], breaks = meth.breaks, include.lowest = T); as.numeric(table(y))}else {rep(0,length(meth.breaks)-1)}})
  rownames(chrom.mat) = levels(cut(c(0,0), breaks=round(meth.breaks, 2), include.lowest = T))
  meth.colors <- brewer.pal(nrow(chrom.mat), "Spectral")[nrow(chrom.mat):1]
  cex.plot <- 1.5
  barplot(chrom.mat, space = 0, col = meth.colors, border = F, las = 2, xaxt = "n", cex.lab = cex.plot, cex.axis = cex.plot, ylab = "Number of reads", cex.main = cex.plot * 1.3, xlab = paste(elementMetadata(check.gr)$type, elementMetadata(check.gr)$Name, elementMetadata(check.gr)$locus_type, sep = ", "))
  axis(1, at = 0, labels = paste("-", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = as.integer(updown/meths$binlen), labels = "start", lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)), labels = paste("+", updown/1000, "Kb", sep = ""), lwd.ticks = 5, line =0.5)
  axis(1, at = length(colnames(chrom.mat)) - as.integer(updown/meths$binlen), labels = "end", lwd.ticks = 5, line =0.5)
}


#check.pos <- as.numeric(getupdown.araport(ara.ind, aradf = araport.tes))
ara.ind <- 161
ara.ind <- 106
ara.ind <- 14468
meth.region.plot(check.gr = araport.gff[ara.ind], meths = meths, updown = 2000)
meth.region.plot(check.gr = araport.gff[ara.ind], meths = meths.sperm, updown = 2000)
meth.region.plot(check.gr = araport.gff[ara.ind], meths = meths.cg, updown = 2000)
araport.gff[ara.ind]

start(araport.gff[ara.ind]) - 2000



meth.region.plot(check.gr = GRanges("Chr1", IRanges(11705950,11724300)), meths = meths, updown = 2000)

ara.ind <- 289

meth.region.plot(check.gr = araport.gff[ara.ind], meths = meths, updown = 2000)

zones=matrix(c(1,2,3), ncol=1, byrow=T)
layout(zones)
meth.region.plot.nolegend(check.gr = araport.gff[ara.ind], meths = meths.cg, updown = 2000)
meth.region.plot.nolegend(check.gr = araport.gff[ara.ind], meths = meths.chg, updown = 2000)
meth.region.plot.nolegend(check.gr = araport.gff[ara.ind], meths = meths.chh, updown = 2000)


