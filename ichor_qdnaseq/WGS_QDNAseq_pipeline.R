library("Biobase")
library("QDNAseq")
library("DNAcopy") 
library("CGHcall")
#library("ACE")


##### initialization #####

args = commandArgs(trailingOnly=TRUE)

path_to_bam_file = args[1]
inputPath = normalizePath(dirname(path_to_bam_file))

NAMEEE_intermediary1 = gsub("\\", "/", normalizePath(path_to_bam_file), fixed = TRUE) # for window path
NAMEEE_intermediary2 = sub("*/*.bam", "", NAMEEE_intermediary1)

NAMEEE = sub(".*/", "", NAMEEE_intermediary2)

output_relative_Path = args[2]
outputPath = normalizePath(output_relative_Path)

bin_size = args[3]


inputPath
NAMEEE
outputPath
bin_size


##### QDNAseq #####

bin_size = as.numeric(bin_size)

bins <- getBinAnnotations(binSize=bin_size)

paste0(inputPath,"/",NAMEEE,".bam")

readCounts <- binReadCounts(bins, bamfiles=paste0(inputPath,"/",NAMEEE,".bam"))

# exportBins(readCounts, file="readcount_unfiltered.tsv", 
#            format="tsv", type=c("copynumber"), logTransform = FALSE)
 
readCountsFiltered <- applyFilters(readCounts)
readCountsFiltered <- estimateCorrection(readCountsFiltered)

readCountsFiltered <- applyFilters(readCounts, chromosomes=NA, residual=TRUE, blacklist=TRUE)
#readCountsFiltered <- applyFilters(readCounts, chromosomes=c("Y", "MT"),
 #                                  residual=TRUE, blacklist=TRUE)

 
# exportBins(readCountsFiltered, file="D038R4_readcount_filtered_10kb.tsv", 
#            format="tsv", type=c("copynumber"), logTransform = FALSE)

copyNumbers <- correctBins(readCountsFiltered)
 
copyNumbersNormalized <- normalizeBins(copyNumbers)
 
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)


##### Downstream Analysis #####

#1. calculate CN by segments
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun = 'none',
                                    alpha = 0.05, nperm = 10000, p.method = "hybrid",
                                    min.width=5, kmax=25, nmin=200, eta=0.05,
                                    trim = 0.025, undo.splits = "sdundo", undo.SD=1)
 
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)
 
exportBins(copyNumbersSegmented, file=paste0(outputPath,"/",NAMEEE,".bam.CN"), format="tsv", type=c("copynumber")) 
exportBins(copyNumbersSegmented, file=paste0(outputPath,"/",NAMEEE,".bam.seg"), format="tsv", type=c("segments"))

X = read.table(paste0(outputPath,"/",NAMEEE,".bam.CN"), sep = "\t", header = TRUE)
Y = read.table(paste0(outputPath,"/",NAMEEE,".bam.seg"), sep = "\t", header = TRUE)
 
X = cbind(X, Y[,5])
X = X[,-4]
X = X[,-1]
#print(X)

colnames(X) <- c("chromosome", "start", "ratio", "ratio_median")
X[,1] <- gsub('X','23',X[,1])
X[,1] = as.numeric(as.character(X[,1]))

#2. create the input file of scarHRD, .ratio
write.table(X, file = paste0(outputPath,"/",NAMEEE,".bam.ratio.tsv"), sep = "\t", row.names = FALSE)

#3. make the calls of CNA
copyNumbersCalled <- callBins(copyNumbersSegmented)
exportBins(copyNumbersCalled, file=paste0(outputPath,"/",NAMEEE,".bam.cna.seg"), format="seg")   # exporta los cn en los segmentos detectados alterados. Este es el input para gistic, aunque no tiene el mismo tamaño de seg por muestra
exportBins(copyNumbersCalled, file=paste0(outputPath,"/",NAMEEE,".bam.cna.seg.called"), format="tsv", type=c("calls"))
#exportBins(copyNumbersCalled, file=paste0(outputPath,"/",NAMEEE,".bam.seg.forgistic"), format="vcf", type=c("segments"))
