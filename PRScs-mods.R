library(data.table)
library(dplyr)

sumstats <- fread("~/Analysis/GWASData/XAuto/GENESIS_2/MADCaPMerged_XAuto_LMM.txt")
sumstats$chr <- as.integer(sumstats$chr)
conv <- fread("~/Analysis/Misc/MADCaP_imputed_snp_rsid_int.txt")
sumstats <-  sumstats %>% inner_join(conv, by=c('chr'='chr', 'pos'='pos'))
sumstats$rsid[sumstats$rsid == "."] <- NA; sumstats <- sumstats[!is.na(sumstats$rsid),]
split_cols <- data.frame(do.call('rbind', strsplit(as.character(sumstats$variant.id),':',fixed=TRUE)))
sumstats_1 <- data.frame(SNP = sumstats$rsid, A1 = split_cols$X4, A2 = split_cols$X3, BETA = sumstats$Est, P = sumstats$Score.pval)
# sumstats <- bigreadr::fread2("~/Analysis/GWASData/XAuto/GENESIS_2/MADCaPMerged_XAuto_GWASv2.txt")
# sumstats_1 <- data.frame(SNP = sumstats$rsid, A1 = sumstats$a1, A2 = sumstats$a0, BETA = sumstats$beta, P = sumstats$p)
fwrite(sumstats_1, "~/Analysis/PRS/PRScsx/MADCaPMerged_XAuto_GWAS_PRScsx.txt", sep = "\t")

#convert gwide bim
bimfile <- fread(paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed.bim"))
bimfile <-  bimfile %>% inner_join(conv, by=c('V1'='chr', 'V4'='pos'))
bimfile$rsid[bimfile$rsid == "."] <- NA; bimfile <- bimfile[!is.na(bimfile$rsid),]
bimfile <- data.frame(bimfile$V1, bimfile$rsid, bimfile$V3, bimfile$V4, bimfile$V5, bimfile$V6)
fwrite(bimfile, paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed_rsid.bim"), sep = "\t")



#convert bim to rsid cols
conv <- fread("~/Analysis/Misc/MADCaP_imputed_snp_rsid_int.txt")
for (i in 1:22)
  {
bimfile <- fread(paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed_", i, ".bim"))

#conv <- conv %>% mutate(chr = str_replace(chr, "chr", ""))
#conv$chr <- as.integer(conv$chr)
#fwrite(conv, "/ssd/LachanceLab/GWASData/Misc/MADCaP_imputed_snp_rsid_int.txt", sep = "\t")

#sumstats <-  sumstats %>% bind_rows(semi_join(conv, sumstats, by = c("chr", "pos")))
bimfile <-  bimfile %>% inner_join(conv, by=c('V1'='chr', 'V4'='pos'))
bimfile$rsid[bimfile$rsid == "."] <- NA; bimfile <- bimfile[!is.na(bimfile$rsid),]
bimfile <- data.frame(bimfile$V1, bimfile$rsid, bimfile$V3, bimfile$V4, bimfile$V5, bimfile$V6)
fwrite(bimfile, paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed_", i, "_rsid.bim"), sep = "\t")
}


#After calculaion, prepare SNP iDs for PRS calculation
# effect_sizes <- fread("/Users/joelvaz/Analysis/PRS/PRScsx/MADCaPMerged_UKBB_effect_sizes_combined.txt", header = F)
# #conv <- fread("~/Analysis/Misc/MADCaP_imputed_snp_rsid_int.txt")
# effect_sizes$V7 <-  paste0("chr", effect_sizes$V1, ":", effect_sizes$V3, ":", effect_sizes$V5, ":", effect_sizes$V4)
# fwrite(effect_sizes, "/Users/joelvaz/Analysis/PRS/PRScsx/MADCaPMerged_UKBB_effect_sizes_combined_SNPIDs.txt", sep = "\t")


#After calculaion, prepare SNP iDs for PRS calculation
effect_sizes <- fread("/Users/joelvaz/Analysis/PRS/PRScsx/MADCaPMerged_UKBB_effect_sizes_combined.txt", header = F)
bimfile <- fread(paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed.bim"))
#get the conversion done from above ------
effect_sizes <- effect_sizes %>% inner_join(bimfile, by=c('V2'='rsid'))
effect_sizes_split <- data.frame(do.call('rbind', strsplit(as.character(effect_sizes$V2.y),':',fixed=TRUE)))
effect_sizes <- data.frame(effect_sizes$V1.x, effect_sizes$V2.y, effect_sizes_split$X2, effect_sizes_split$X4,
                           effect_sizes_split$X3, effect_sizes$V6.x)
effect_sizes = effect_sizes[!duplicated(effect_sizes$effect_sizes.V2.y),]
fwrite(effect_sizes, "/Users/joelvaz/Analysis/PRS/PRScsx/MADCaPMerged_UKBB_effect_sizes_combined_SNPIDs.txt", sep = "\t")









