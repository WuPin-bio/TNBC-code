#script for figure4
#.....................................reconstruct tree...........................
#step 1: generate input
library(magrittr)
options(stringsAsFactors=F)
rm(list = ls())
setwd('./')
indir <- './'
df_sam <- read.csv('sample_tumor.csv', row.names=1) %>% subset(Tumor > 40)
#df_sam$new_id <- paste0(rownames(df_sam), '-', df_sam$Tumor)
patient <- c('P01','P02','P03','P04','P05','P06','P07','P08','P09','P10')
all_file <- dir(indir, full.names=T) %>% grep('_', ., value=T)
str(all_file)
for (p in patient) {
  print(paste0('Running for ', p))
  df_f <- data.frame(files = grep(paste0(p, '_'), all_file, value=T))
  df_f$sample <- sub('\\..*', '', basename(df_f$files))
  df_f[nrow(df_f) + 1, ] <- c('Normal.hg38.100k.k100.nobad.varbin.data.txt', 'Normal')

  len <- (max(nchar(df_f$sample)) + 5 - nchar(df_f$sample))
  df_f$n_blank <- sapply(len, function(x) paste(rep(' ', x), collapse=''))
  # sample filtering
  df_f <- subset(df_f, sample %in% c(rownames(df_sam), 'Normal'))
  #  df_f$new_id <- df_sam[df_f$sample, 'new_id']
  # header
  n_var <- system(paste("grep -vP '^24|chrom'", df_f$files[1], '| wc -l'), intern=T)
  n_sam <- nrow(df_f)
  lines <- paste0(n_sam, '      ', n_var)
  for ( i in rev(1:nrow(df_f)) ) {
    cnv <- data.table::fread(df_f$files[i], data.table=F, showProgress=F) %>%
		subset(chrom != 24)
    # classify: 1 ( cnv < 1.5); 2 ( 1.5 <= cnv < 2.5 ); 3 ( cnv >= 2.5)
    cnv$seg.mean.LOWESS <- 2 * cnv$seg.mean.LOWES
    cnv[cnv$seg.mean.LOWESS < 1.5, 'seg.mean.LOWESS'] <- 1
    cnv[cnv$seg.mean.LOWESS >= 1.5 & cnv$seg.mean.LOWESS < 2.5, 'seg.mean.LOWESS'] <- 2
    cnv[cnv$seg.mean.LOWESS >= 2.5, 'seg.mean.LOWESS'] <- 3
    # rename sample
	sam <- df_f$sample[i]
    val <- paste(cnv$seg.mean.LOWESS, collapse = '')
    lines <- c(lines, paste0(sam, df_f$n_blank[i], val))
  }
   file.remove(paste0(p, '_forPhylip.txt'))
   writeLines(lines, con = paste0(p, '_forPhylip.txt'))
  rm(df_f, n_var, n_sam, len, lines)
}

#.....................................run GISTIC to obtain focal-level SCNA...........................
#step1: varbin worflow
#!/usr/bin/env python
import sys
import time
from operator import itemgetter
import bisect

def main():

	bincount = 100000
	MAP = open("./6-output-hg38.chrom.mappable.bowtie.k100.txt", "r")
	GOOD = open("./5-output-hg38.goodzones.umap.k100.bed", "r")
	CHROMLEN = open("./4-output-hg38.chrom.sizes.txt", "r")
	OUTFILE = open("./hg38.bin.boundaries.100k.bowtie.100.txt", "w")

	#OUTFILE.write("chrom\tbin.start.chrompos\tbin.start.abspos\tbin.end.chrompos\tbin.length\tmappable.positions\n")

	chromlen = dict()
	for x in CHROMLEN:
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisChromlen = int(arow[1])
		thisAbspos = int(arow[2])
		chromlen[thisChrom] = [thisChromlen, thisAbspos]

	chromarray = []
	chroms = dict()
	totalLength = int(0)
	for x in MAP:
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisLength = int(arow[1])
		if thisChrom == "chrM":
			continue
		totalLength += thisLength
		chroms[thisChrom] = thisLength

	bincountUsed = 0
	for k, v in chroms.items():
		chromBincount = float(bincount) * (float(v) / float(totalLength))
		i = int(chromBincount)
		bincountUsed += i 
		r = chromBincount - i
		chromarray.append([k, i, r])

	a = []
	for i in chromarray:
		bisect.insort(a, (-i[2], i))

	chromarray = []
	for j in a:
		chromarray.append(j[1]) ##  This j[1] index has to match the index of i in the bisect.insort 2nd parameter.

	a = []
	
	print(chromarray)	

	remain = bincount - bincountUsed
	for i in range(remain):
		chromarray[i][1] += 1

	chroms2 = dict()
	for i in range(len(chromarray)):
		chromlength = chroms[chromarray[i][0]]
		chrombins = chromarray[i][1]
		binlength = float(chromlength) / float(chrombins)
		chroms2[chromarray[i][0]] = [chrombins, binlength]

	
	chromorder = sorted(chroms2.keys())
	print(chromorder)

	print 
	print("Starting to get bin boundaries")
	print
	goodEOF = False
	for chrom in chromorder:
		print(chrom)
		firstBin = True
		#  position GOOD file
		x = GOOD.readline()
		arow = x.rstrip().split("\t")
		thisChrom = arow[0]
		thisStart = int(arow[1])
		thisEnd = int(arow[2])
		print("new chrom", arow)
		while thisChrom != chrom:
			x = GOOD.readline()
			arow = x.rstrip().split("\t")
			if len(x) == 0:
				goodEOF = True
				break
			thisChrom = arow[0]
			thisStart = int(arow[1])
			thisEnd = int(arow[2])
			print("position chrom", arow)
		if goodEOF:
			break
		print("after position")
		currentStart = thisStart 
		chromBincount = chroms2[chrom][0]
		chromBinlength = chroms2[chrom][1]
		chromExcess = chromBinlength - int(chromBinlength)
		currentExcess = 0.0
		thisBincount = 0
		binStart = 0
		binEnd = 0
		currentLength = 0
		print(chromBincount, chromBinlength)
		while thisBincount < chromBincount:
			thisBincount += 1
			print("thisBincount", thisBincount)
			thisBinlength = int(chromBinlength)
			currentExcess += chromExcess
			print("currentExcess", currentExcess)
			if currentExcess >= 1.0:
				currentExcess -= 1.0
				thisBinlength += 1
			print("thisBinlength", thisBinlength)
			binStart = currentStart
			currentLength = 0
			if (binStart + thisBinlength) < thisEnd:
				print("got bin from current GOOD")
				binEnd = binStart + thisBinlength
				currentStart = binStart + thisBinlength
				currentLength = thisBinlength
			else:
				print("getting more GOOD")
				currentLength += (thisEnd - currentStart)
				print(thisEnd, currentStart, currentLength, thisBinlength)
				while currentLength < thisBinlength:
					x = GOOD.readline()
					print(x)
					if len(x) == 0:
						goodEOF = True
						break
					arow = x.rstrip().split("\t")
					thisChrom = arow[0]
					thisStart = int(arow[1])
					thisEnd = int(arow[2])
					print("adding length", arow)
					if thisChrom != chrom:
						print("ERROR: Past end of chrom.", thisChrom, chrom)
						break

					if (thisEnd - thisStart) < (thisBinlength - currentLength):
						currentLength += (thisEnd - thisStart)
					else:
						currentStart = thisStart + (thisBinlength - currentLength)
						currentLength = thisBinlength
					print("currentLength", currentLength)
				binEnd = currentStart

			if firstBin:
				binStart = 0
			if thisBincount == chromBincount:
				binEnd = chromlen[chrom][0]

			binStartAbspos = chromlen[chrom][1] + binStart

			OUTFILE.write(chrom)
			OUTFILE.write("\t")
			OUTFILE.write(str(binStart))
			OUTFILE.write("\t")
			OUTFILE.write(str(binStartAbspos))
			OUTFILE.write("\t")
			OUTFILE.write(str(binEnd))
			OUTFILE.write("\t")
			OUTFILE.write(str(binEnd - binStart))
			OUTFILE.write("\t")
			OUTFILE.write(str(currentLength))
			OUTFILE.write("\n")

			firstBin = False

			if goodEOF:
				break


	OUTFILE.close()
	MAP.close()
	GOOD.close()


if __name__ == "__main__":
	main()


#!/usr/bin/env python
import sys
import pysam
def main():

	infilename = sys.argv[1]
	outfilename = sys.argv[2]
	statfilename = sys.argv[3]

	chrominfo = fileToDictionary("./4-output-hg38.chrom.sizes.txt", 0)
	bins = fileToArray("./hg38.bin.boundaries.100k.bowtie.100.sorted.txt", 0)
	#bins = fileToArray(sys.argv[4],0)
	#INFILE = open(infilename, "r")
	bf = pysam.AlignmentFile(infilename,"rb",threads=1)
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	print (len(binCounts))
	print (len(bins))

	counter = 0
	totalReads = 0
	prevChrompos = ""
    
    
	for i in range(len(bins)):
		contig ,start,stop= bins[i][0],int(bins[i][1]),int(bins[i][3])
		region_count= bf.count(contig=contig, start=start, stop=stop)
		binCounts[i]=region_count
		totalReads += region_count
		counter += region_count
		
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")

	binCounts.sort()
	
	STATFILE.write("TotalReads\tMedianBinCount\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(binCounts[int(len(bins)/2)]))
	STATFILE.write("\n")

	#INFILE.close()
	OUTFILE.close()
	STATFILE.close()


def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")

	rd = dict()
#	input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
		if id in rd:
			#rd[id].append(arow)
			print ("duplicate knowngene id = " + id)
			print ("arow =   " + str(arow) )
			print ("rd[id] = " + str(rd[id]))
		else:
			rd[id] = arow
		
	input.close()
	return(rd)


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []

	for i in range(skipFirst):
		input.readline()

	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
		
	input.close()
	return(ra)


if __name__ == "__main__":
	main()

library("DNAcopy", lib.loc="./R/library/")

lowess.gc <- function(jtkx, jtky) {
        jtklow <- lowess(jtkx, log(jtky), f=0.05)
        jtkz <- approx(jtklow$x, jtklow$y, jtkx)
        return(exp(log(jtky) - jtkz$y))
}

cbs.segment01 <- function(indir, outdir, bad.bins, varbin.gc, varbin.data, sample.name, alt.sample.name, alpha, nperm, undo.SD, min.width) {
	gc <- read.table(varbin.gc, header=T)
	bad <- read.table(bad.bins, header=F)

	chrom.numeric <- substring(gc$bin.chrom, 4)
	chrom.numeric[which(gc$bin.chrom == "chrX")] <- "23"
	chrom.numeric[which(gc$bin.chrom == "chrY")] <- "24"
	chrom.numeric <- as.numeric(chrom.numeric)

	thisRatio <- read.table(paste(indir, varbin.data, sep="/"), header=F) 
	names(thisRatio) <- c("chrom", "chrompos", "abspos", "bincount", "ratio")
	thisRatio$chrom <- chrom.numeric
	a <- thisRatio$bincount + 1
	thisRatio$ratio <- a / mean(a)
	thisRatio$gc.content <- gc$gc.content
	thisRatio$lowratio <- lowess.gc(thisRatio$gc.content, thisRatio$ratio)

	a <- quantile(gc$bin.length, 0.985)
	thisRatioNobad <- thisRatio[which(bad[, 1] == 0), ]
	
	set.seed(25) 
	CNA.object <- CNA(log(thisRatio$lowratio, base=2), thisRatio$chrom, thisRatio$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatio), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatio$seg.mean.LOWESS <- m[, 1]

	chr <- thisRatio$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatio$abspos[which(chr != chr.shift) + 1], thisRatio$abspos[nrow(thisRatio)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
	y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
	y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")

	png(paste(outdir, "/", sample.name, ".wg.png", sep=""), height=400, width=600)
	par(mar=c(5.1,4.1,4.1,4.1))
	plot(x=thisRatio$abspos, y=thisRatio$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", yaxt="n", ylab="Ratio", col="#CCCCCC", cex=0.5)
	axis(1, at=x.at, labels=x.labels)
	axis(2, at=y.at, labels=y.labels)
	lines(x=thisRatio$abspos, y=thisRatio$lowratio, col="#CCCCCC", cex=0.5)
	points(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	lines(x=thisRatio$abspos, y=thisRatio$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	write.table(thisRatio, sep="\t", file=paste(outdir, "/", sample.name, ".hg38.100k.k100.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg38.100k.k100.varbin.short.txt", sep=""), quote=F, row.names=F) 



	set.seed(25) 
	CNA.object <- CNA(log(thisRatioNobad$lowratio, base=2), thisRatioNobad$chrom, thisRatioNobad$chrompos, data.type="logratio", sampleid=sample.name) 
	smoothed.CNA.object <- smooth.CNA(CNA.object) 
	segment.smoothed.CNA.object <- segment(smoothed.CNA.object, alpha=alpha, nperm=nperm, undo.splits="sdundo", undo.SD=undo.SD, min.width=min.width) 
	thisShort <- segment.smoothed.CNA.object[[2]]

	m <- matrix(data=0, nrow=nrow(thisRatioNobad), ncol=1)	
	prevEnd <- 0
	for (i in 1:nrow(thisShort)) {
		thisStart <- prevEnd + 1
		thisEnd <- prevEnd + thisShort$num.mark[i]
		m[thisStart:thisEnd, 1] <- 2^thisShort$seg.mean[i]
		prevEnd = thisEnd
	}
	
	thisRatioNobad$seg.mean.LOWESS <- m[, 1]

	chr <- thisRatioNobad$chrom
	chr.shift <- c(chr[-1], chr[length(chr)])
	vlines <- c(1, thisRatioNobad$abspos[which(chr != chr.shift) + 1], thisRatioNobad$abspos[nrow(thisRatioNobad)])
	hlines <- c(0.5, 1.0, 1.5, 2.0)
	chr.text <- c(1:22, "X", "Y")
	chr.text <- chr.text[c(1:14, 16, 18, 20, 22:24)]
	vlines.shift <- c(vlines[-1], 4*10^9)
	chr.at <- vlines + (vlines.shift - vlines) / 2
	chr.at <- chr.at[c(1:14, 16, 18, 20, 22:24)]
	x.at <- c(0, 0.5, 1, 1.5, 2, 2.5, 3) * 10^9
	x.labels <- c("0", "0.5", "1.0", "1.5", "2.0", "2.5", "3.0")
	y.at <- c(0.005, 0.020, 0.100, 0.500, 2.000)
	y.labels <- c("0.005", "0.020", "0.100", "0.500", "2.000")

	png(paste(outdir, "/", sample.name, ".wg.nobad.png", sep=""), height=400, width=600)
	par(mar=c(5.1,4.1,4.1,4.1))
	plot(x=thisRatioNobad$abspos, y=thisRatioNobad$lowratio, log="y", main=paste(sample.name, alt.sample.name), xaxt="n", xlab="Genome Position Gb", yaxt="n", ylab="Ratio", col="#CCCCCC", cex=0.5)
	axis(1, at=x.at, labels=x.labels)
	axis(2, at=y.at, labels=y.labels)
	lines(x=thisRatioNobad$abspos, y=thisRatioNobad$lowratio, col="#CCCCCC", cex=0.5)
	points(x=thisRatioNobad$abspos, y=thisRatioNobad$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	lines(x=thisRatioNobad$abspos, y=thisRatioNobad$seg.mean.LOWESS, col="#0000AA", cex=0.5)
	abline(h=hlines)
	abline(v=vlines)
	mtext(chr.text, at = chr.at)
	dev.off()
	
	write.table(thisRatioNobad, sep="\t", file=paste(outdir, "/", sample.name, ".hg38.100k.k100.nobad.varbin.data.txt", sep=""), quote=F, row.names=F) 
	write.table(thisShort, sep="\t", file=paste(outdir, "/", sample.name, ".hg38.100k.k100.nobad.varbin.short.txt", sep=""), quote=F, row.names=F) 

}
cbs.segment01(indir=".", outdir=".", bad.bins="Supplementaryhg19.50k.k50.bad.bins.txt", varbin.gc="filepath/hg38.varbin.gc.content.100k.bowtie.k100.txt", varbin.data="P01_A_10.pysam.count", sample.name="P01_A_10", alt.sample.name="", alpha=0.05, nperm=1000, undo.SD=1.0, min.width=5)




#step2:run GISTIC
#source activate rna_env
path=./
source /hwfsxx1/ST_HN/PUB/conda/install/bin/activate rna_env

/hwfsxx1/ST_HN/PUB/conda/install/envs/rna_env/bin/gistic2 \
-b ${path} \
-seg ${path}/data_cna_hg38-modify.seg  \
-refgene ./hg38.UCSC.add_miR.160920.refgene.mat \
-conf 0.99 \
-genegistic 1 \
-savegene 1 \
-broad 1


#.....................................focal-level SCNA with phylogenetic signals...........................
library(pheatmap)
setwd(',.')
#sample type feature1 feature2 ...
# s1    high    10       12
data<-read.table('Amp-parcor.txt',header=T,row.names = 1,check.names = F)
data<-data.frame(data)
annotation<-read.table('Amp-parPadjust.txt',header=T,row.names = 1,check.names = F)
annotation<-data.frame(annotation)
N_pairs<-data$N_pairs
type<-data$type
annotation_row<-data.frame(sample_DNA=rownames(data),N_pairs=N_pairs,type=type)
rownames(annotation_row)<-annotation_row$sample_DNA
annotation_row$sample_DNA<-NULL
df<-data[,3:length(data)]
bk<-c(seq(-1,-0.1,by=0.1),seq(0,1,by=0.1))
pheatmap(df,
         scale = "none",main="Partial Correlation Coefficient",
         display_numbers = matrix(ifelse(annotation<0.1, ifelse(annotation<0.05,ifelse(annotation<0.01,"***","**"),"*"), ""), nrow(data)),
         cluster_rows = T,show_colnames = T,show_rownames = T,
         cluster_cols = F,border_color = NA,annotation_row = annotation_row,
         color = c(colorRampPalette(colors = c("#137AB5","white"))(length(bk)/2),colorRampPalette(colors =c("white","#CC0724"))(length(bk)/2)),legend_breaks=seq(-1,1,2),breaks=bk
)


#...........................survival...........................
#survival analysis
#install.packages(c("survival", "survminer"))
library("survival")
library("survminer")

setwd('./')
data<-read.table('data.txt',header=T,row.names = 1,check.names = F)
#data.txt
#patient_id   PFI    PFI_time    geneA_mutation_status    geneB_mutation_status
#P1            1       100              0                         1
#P2            0        50              1                         1

gene_list<-colnames(data)[3:ncol(data)]

pdf('PFI-survival.pdf')
i=2

for (gene in gene_list){
print(gene)
i=i+1
fit <- survfit(Surv(PFI_time,PFI) ~ data[,i], data = data)
print(fit)
summary(fit)$table

p<-ggsurvplot(fit,data=data,
           pval = TRUE, conf.int = TRUE,
           risk.table = TRUE, 
           risk.table.col = "strata", 
           linetype = "strata", 
           surv.median.line = "hv",
           ggtheme = theme_bw(), 
           palette = c("#E7B800", "#2E9FDF"),
           title=paste(gene,'PFI'))

print(p)
}
dev.off()