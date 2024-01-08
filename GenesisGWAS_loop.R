library(GENESIS)
library(SNPRelate)
library(GWASTools)
library(data.table)
library(ggplot2)
library(ggsci)

for (chr in 1:23) { 
gdsfmt::showfile.gds(closeall=TRUE)
PathToData <- paste0("~/Analysis/GWASData/XAuto/Split-imputed/MADCaPMerged_XAuto_", chr)
GDSautopath <- snpgdsBED2GDS(bed.fn = paste0(PathToData, ".bed"), bim.fn = paste0(PathToData, ".bim"),
                             fam.fn = paste0(PathToData, ".fam"),  out.gdsfn = paste0(PathToData, ".gds"))
GDSauto <- GenotypeData(GdsGenotypeReader(GDSautopath))
#reload GDS file
# gdsfmt::showfile.gds(closeall=TRUE)

#Running PC-AiR
mypcair <- pcair(GDSauto, num.cores = 9L, autosome.only = F) #, snp.include = prunedSNP
saveRDS(mypcair, paste0(PathToData, "_mypcair.rds"))
# mypcair <- readRDS(paste0(PathToData, "_mypcair.rds"))

# Running PC-Relate
GDSauto <- GenotypeBlockIterator(GDSauto, snpBlock=500) #, snpInclude = prunedSNP
mypcrelate <- pcrelate(GDSauto, pcs = mypcair$vectors[,1:2], training.set = mypcair$unrels, BPPARAM = BiocParallel::MulticoreParam(9))
saveRDS(mypcrelate, paste0(PathToData, "_mypcrel.rds"))

#LMM
#loading phenotype and covariate file
# mypcair <- readRDS(paste0(PathToData, "_mypcair.rds"))
# mypcrelate <- readRDS(paste0(PathToData, "_mypcrel.rds"))
pheno <- fread("~/Analysis/GWASData/XAuto/PrCaPheno.pheno")
pheno <- pheno[pheno$IID %in% mypcair$sample.id,] 
annot <- data.frame(scanID = mypcair$sample.id, pc1 = mypcair$vectors[,1], pheno = pheno[,3]); 
annot <- ScanAnnotationDataFrame(annot)
myGRM <- pcrelateToMatrix(mypcrelate)[annot$scanID, annot$scanID]

nullmod <- fitNullModel(annot, outcome = "PrCa", covars = "pc1", cov.mat = myGRM, family = "gaussian")
genoIterator <- GenotypeBlockIterator(GDSauto)
assoc <- assocTestSingle(genoIterator, null.model = nullmod, BPPARAM = BiocParallel::MulticoreParam(9))
#assoc <- assoc[order(assoc$Score.pval),]
fwrite(assoc, paste0(PathToData, "_LMM.txt"), sep = "\t")
rm(list = ls()); gc()
}

rm(list = ls()); gc()
assoc <- fread("/ssd/LachanceLab/GWASData/POLMM_input/MADCaP_West_Agg_Merged_AllChr_LMM.LZoom")
ggplot(assoc, aes(x=V2, y=-log10(V5), color=V1)) + geom_point(size = 0.5) + labs(x = "Chromosome",  y = "-log10(p)") + 
  theme_minimal() + scale_color_igv() + 
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank(),
                          panel.grid.minor.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_text(size=12), axis.title = element_text(size=18))
ggsave("/ssd/LachanceLab/Plots/MADCaP_West/Manhtt_X_Military.png", units = "cm", height = 7, width = 22, dpi = 1000, bg = "white")
