library(methylumi)
library(wateRmelon)
library(FDb.InfiniumMethylation.hg19)
library(pcaMethods)

setwd("/data/omicscore/Belsky-Belsky-20200413/")

#sampleList <- gsub("_Grn.idat", "", list.files("./rawData/", pattern="*_Grn.idat"))
rawData <- readEPIC(idatPath="./rawData/")

detectionPvals <- pvals(rawData)
numberOfBadProbes <- apply(detectionPvals, 2, function(x) { length(which(x >= 0.05)) })
save(detectionPvals, file="./results/detectionPvals.Rdata")


# Load in the phenotype data:
sampleData.df <- read.table("./rawData/Calerie_EPIC_samplesheet_All runs.csv", sep=",", header=T, skip=7)

#Fix some mis-labeled sample time points:
sampleData.df$Sample_Name <- gsub("\\sMonth", "", sampleData.df$Sample_Name)

# Remove the technical replicates
samplesToRemove <- paste(
	sampleData.df$Sentrix_ID[grep("_rep", sampleData.df$Sample_Name)],
	sampleData.df$Sentrix_Position[grep("_rep", sampleData.df$Sample_Name)],
	sep="_"
)
sampleData.df <- sampleData.df[-grep("_rep", sampleData.df$Sample_Name),]

sampleData.df$Visit <- gsub(".+_", "", sampleData.df$Sample_Name)
sampleData.df$CalorieID <- gsub("_.+", "", sampleData.df$Sample_Name)

rownames(sampleData.df) <- paste(sampleData.df$Sentrix_ID, sampleData.df$Sentrix_Position, sep="_")
rawData <- rawData[,rownames(sampleData.df)]

pData(rawData)$Sample_Name <- sampleData.df$Sample_Name
pData(rawData)$CalorieID <- sampleData.df$CalorieID
pData(rawData)$Visit <- sampleData.df$Visit
pData(rawData)$Run <- sampleData.df$Run
pData(rawData)$Sample_Plate <- sampleData.df$Sample_Plate

normData <- normalizeMethyLumiSet(rawData)

save(normData, file="./results/normData.Rdata")


## Get QC probe names
qcProbeNames <- row.names(methylated(QCdata(normData)))

## generate a matrix of beta values from the QC probes
qc.mat <- methylated(QCdata(normData)) / ( methylated(QCdata(normData)) + unmethylated(QCdata(normData)))

## We only want QC probes that have 'Norm' in the name
qc.mat <- qc.mat[grep("Norm", qcProbeNames),]

## Replace NAs with NAs
qc.mat[which(is.na(qc.mat))] <- NA

pc <- pca(t(qc.mat), nPcs=50, method="ppca")
save(pc, file="./results/PCA_ControlProbes.Rdata")

sd.45 <- sDev(pc)
sd.45 <- sd.45**2 / (sum(sd.45**2))
cutoff <- as.numeric(which(cumsum(sd.45) >= 0.9)[1])

pca.df <- as.data.frame(scores(pc)[,1:cutoff])
pca.df$barcode <- rownames(pca.df)
pca.df <- pca.df[,c(25,1:24)]

write.table(pca.df, file="./results/CALERIE_Control_PCs.txt", sep="\t", row.names=F, quote=F)

detectionPvals <- pvals(rawData)

betas <- betas(normData)

#Cleanup bad probes:
betas[which(detectionPvals >= 0.05)] <- NA

apply(betas, 1, function(x) { length(which(is.na(x))) }) -> bp

betas <- betas[-which(bp > 0.05 * ncol(betas)),]

geneticProbes <- betas[grep("^rs", rownames(betas)),]

geneticCorrelation.mat <- cor(geneticProbes, use="complete.obs")


summary(sapply(unique(pData(normData)$CalorieID), function(x) {
	barcodes <- as.character(pData(normData)$barcode[which(pData(normData)$CalorieID == x)])
	if( length(barcodes) == 1 ) {
		1
	} else {
		sub.mat <- geneticCorrelation.mat[barcodes, barcodes]
		min(sub.mat)
	}
}))

betas <- betas[grep("^cg", rownames(betas)),]
phenoData.df <- pData(normData)

save(betas, phenoData.df, file="./results/CALERIE_EPIC_data.Rdata")
