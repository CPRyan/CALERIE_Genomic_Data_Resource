library(methylumi)
library(wateRmelon)
library(FDb.InfiniumMethylation.hg19)
library(pcaMethods)
pacman::p_load(methylumi, wateRmelon, FDb.InfiniumMethylation.hg19, pcaMethods, tidyverse, readr)

sample_file <-read_csv(here::here("Data/IDATs", "adipose_sample_annotation_file.csv")) %>% 
  filter(!is.na(barcode))
rawData <- readEPIC(idatPath= here::here("Data/IDATs/idatFiles"))

# Rename the barcode without the folder name
split_vec <- sapply(strsplit(pData(rawData)$barcode, "/"), `[`, 2)
sampleNames(rawData) <-split_vec

# Filter to samples in our data only
rawData <-rawData[,sampleNames(rawData) %in% sample_file$Basename]
pData(rawData)$DEID <- sample_file$DEID


detectionPvals <- pvals(rawData)
numberOfBadProbes <- apply(detectionPvals, 2, function(x) { length(which(x >= 0.05)) })
save(detectionPvals, file= here::here("Output/Data", "Corcoran_Pipeline/detectionPvals.Rdata"))

normData <- normalizeMethyLumiSet(rawData)

save(normData, file= here::here("Output/Data", "Corcoran_Pipeline/normData.Rdata"))


## Get QC probe names
qcProbeNames <- row.names(methylated(QCdata(normData)))

## generate a matrix of beta values from the QC probes
qc.mat <- methylated(QCdata(normData)) / ( methylated(QCdata(normData)) + unmethylated(QCdata(normData)))

## We only want QC probes that have 'Norm' in the name
qc.mat <- qc.mat[grep("Norm", qcProbeNames),]

## Replace NAs with NAs
qc.mat[which(is.na(qc.mat))] <- NA

pc <- pca(t(qc.mat), nPcs=50, method="ppca")
save(pc, file=here::here("Output/Data", "Corcoran_Pipeline/PCA_ControlProbes.Rdata"))

sd.45 <- sDev(pc)
sd.45 <- sd.45**2 / (sum(sd.45**2))
cutoff <- as.numeric(which(cumsum(sd.45) >= 0.9)[1])

pca.df <- as.data.frame(scores(pc)[,1:cutoff])
pca.df$barcode <- rownames(pca.df)
# pca.df <- pca.df[,c(25,1:24)]
pca.df <-pca.df[,c(2,1)]

write.table(pca.df, file=here::here("Output/Data", "Corcoran_Pipeline/CALERIE_Control_PCs.txt"), sep="\t", row.names=F, quote=F)

detectionPvals <- pvals(rawData)

betas <- betas(normData)

#Cleanup bad probes:
betas[which(detectionPvals >= 0.05)] <- NA

apply(betas, 1, function(x) { length(which(is.na(x))) }) -> bp

betas <- betas[-which(bp > 0.05 * ncol(betas)),]

geneticProbes <- betas[grep("^rs", rownames(betas)),]

geneticCorrelation.mat <- cor(geneticProbes, use="complete.obs")



summary(sapply(unique(pData(normData)$DEID), function(x) {
  barcodes <- as.character(pData(normData)$barcode[which(pData(normData)$DEID == x)])
  if( length(barcodes) == 1 ) {
    1
  } else {
    sub.mat <- geneticCorrelation.mat[barcodes, barcodes]
    min(sub.mat)
  }
}))

# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -0.4381  0.1517  0.4118  0.5486  1.0000  1.0000 

betas <- betas[grep("^cg", rownames(betas)),]

sample_file <-sample_file %>% dplyr::select(-c(OriginalOrder:ExternalSampleID, ProjectID, Plate.ID, ExternalSampleID.1, DNA.Concentration.1, ID.Wanagat:`...19`, -sex))


save(betas, sample_file, file=here::here("Output/Data", "Corcoran_Pipeline/CALERIE_Adipose_Corcoran_data.Rdata"))

write_csv(betas %>% as_tibble(rownames = "probe"), file=here::here("Output/Data", "Corcoran_Pipeline/CALERIE_Adipose_Corcoran_data.csv"))
write_csv(sample_file, file=here::here("Output/Data", "Corcoran_Pipeline/CALERIE_Adipose_Corcoran_sample_file.csv"))

