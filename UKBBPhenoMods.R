library(data.table)
library(reshape)

pheno = data.frame(t(fread("/ssd/LachanceLab/UKBB/Ethnic1Cancer9.txt")))[-1,-3]
colnames(pheno) <- c("ID", "Ethnicity", "Cancer0", "Cancer1", "Cancer2", "Cancer3",
                     "Cancer4", "Cancer5", "Cancer6", "Cancer7", "Cancer8")
sex_fam = data.frame(fread("/ssd/LachanceLab/UKBB/PLINK/UKB_BED_AFTER_QC_removed182IDs_chr22.fam"))

#only males
pheno <- pheno[pheno$ID %in% sex_fam[sex_fam$V5 == 1,]$V1,]
#only africans
pheno_Afr <- pheno[pheno$Ethnicity == "African",]; pheno_Afr <- pheno_Afr[!is.na(pheno_Afr$Ethnicity),]
#only cancer prostate (c6)
pheno_Afr_case <- pheno[pheno$Cancer0 == "C61",] #change pheno to pheno_Afr if africa applicable
pheno_Afr_case <- pheno_Afr_case[!is.na(pheno_Afr_case$Ethnicity),] #remove NAs

#just cases
Afr_PCa_IDs <- data.frame(FID = pheno_Afr_case$ID, IID = pheno_Afr_case$ID, PrCa = 1)

#add controls
pheno_Afr_ctrl <- pheno[sample(which(is.na(pheno$Cancer0)),7296), ] #change pheno to pheno_Afr if africa applicable, and number of samples too
pheno_Afr_ctrl <- data.frame(FID = pheno_Afr_ctrl[,1], IID = pheno_Afr_ctrl[,1], PrCa = 0)
#combine everything
Afr_PCa_IDs <- rbind(Afr_PCa_IDs, pheno_Afr_ctrl)

fwrite(Afr_PCa_IDs, "/ssd/LachanceLab/UKBB/PLINK/UKBB_PrCa_7296.pheno", sep = "\t")

#all pheno
pheno_Afr_all <- data.frame(FID = pheno$ID, IID = pheno$ID, PrCa = ifelse(pheno$Cancer0 == "C61", 1, 0))
pheno_Afr_all$PrCa[is.na(pheno_Afr_all$PrCa)] <- 0

fwrite(pheno_Afr_all, "/ssd/LachanceLab/UKBB/PLINK/UKBB_Afr_PrCa.pheno", sep = "\t")

summstats = fread("/ssd/LachanceLab/GWASData/METAL/META_X_Auto_LDPred2_mod23.txt")
updatedref = fread("/ssd/LachanceLab/GWASData/METAL/updated_ref_alt.txt")
summstats[,"a0"] <- updatedref$Ref
summstats[,"a1"] <- updatedref$Alt
fwrite(summstats, "/ssd/LachanceLab/GWASData/METAL/META_X_Auto_LDPred2_mod23.txt", sep = "\t")

#-------------------------------------------------------------------------------------------------------------
#LDPred2 mods
#snp to rsid
info <- readRDS(runonce::download_file(
  "https://ndownloader.figshare.com/files/25503788",
  fname = "map_hm3_ldpred2.rds"))

info$snpid <- paste0("chr", info$chr, ":", info$pos, ":", info$a0, ":", info$a1)
