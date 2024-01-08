# install.packages("doSNOW", repos = "http://cran.us.r-project.org", dependencies = TRUE)
library(remotes)
#prepare workspace and load bigsnpr
library(bigsnpr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(parallel)
library(data.table)
library(magrittr)
library(dplyr)
library(fmsb)

# clusterEvalQ(cl, .libPaths("~/R_4.2"))

print("Installed packages successfully")

numCores <- detectCores() - 1
print(numCores)
#read in phenotype and covariates
phenotype <- fread("~/Analysis/UKBB/UKBB_PrCa_7296.pheno")
covariate <- fread("~/Analysis/UKBB/UKB_BED_AFTER_QC_removed182IDs_final_subset20000_Age.cov")
# generate required table
pheno <- merge(phenotype, covariate) 

info <- fread("~/Analysis/Misc/hapmap3plus.txt")

#Load summary statistic file
# Read in the summary statistic file
sumstats <- bigreadr::fread2("~/Analysis/GWASData/XAuto/GENESIS_2/MADCaPMerged_XAuto_GWASv2.txt")
#previous sumstats modifications

# LDpred 2 require the header to follow the exact naming
#######("chr", "pos", "rsid", "a1", "a0", "n_eff", "beta_se", "p", "OR", "INFO", "MAF")

# Transform the OR into log(OR)
## sumstats$beta <- log(sumstats$OR)

sumstats <- sumstats[sumstats$rsid%in% info$rsid,]


# Get maximum amount of cores
NCORES <- nb_cores()
# Open a temporary file
tmp <- tempfile(tmpdir = "tmp-data"); #tmp1 <- tempfile(tmpdir = "tmp-data1")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE); #on.exit(file.remove(paste0(tmp1, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
fam.order <- NULL
# preprocess the bed file (only need to do once for each data set)
snp_readBed("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed.bed")
# now attach the genotype object
obj.bigSNP <- snp_attach("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed.rds")
# extract the SNP information from the genotype
map <- obj.bigSNP$map[-3]
names(map) <- c("chr", "rsid", "pos", "a1", "a0")
# perform SNP matching
info_snp <- snp_match(sumstats, map, match.min.prop = 0.2)
# Assign the genotype to a variable for easier downstream analysis
genotype <- obj.bigSNP$genotypes
#genotype <- snp_fastImputeSimple(genotype, ncores = numCores) 
print(genotype)
# Rename the data structures
CHR <- map$chr
POS <- map$pos
# get the CM information from 1000 Genome
# will download the 1000G file to the current directory (".")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "../1000G/", ncores = NCORES)
# calculate LD
for (chr in 1:22) { #change chr
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(info_snp$chr == chr)
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = NCORES,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 1) { #change chr
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}
# We assume the fam order is the same across different chromosomes
fam.order <- as.data.table(obj.bigSNP$fam)
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))
#Perform LD score regression
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
#removing ld NAs from ld
print(paste0("ld before NAs ", length(ld)))
 ld_NA <- which(is.na(ld))
# ld <- ld[!is.na(ld)]
# print(paste0("ld after NAs ", length(ld)))
#removing ld NAs from corr
 corr <- as_SFBM(corr[-ld_NA, -ld_NA], tmp1) #added extra
 print(paste0("cor after NAs ", length(corr)))
#remove NAs from info_snp
 info_snp <- info_snp[-ld_NA, ]
 df_beta_subset <- df_beta[-ld_NA, ]
# print(ld)
# print(nrow(df_beta_subset))
# print((df_beta_subset$beta / df_beta_subset$beta_se)^2)
print(paste0("min LD ", min(ld)))
print(paste0("max chi ", min(df_beta$beta / df_beta$beta_se)^2))
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta$beta / df_beta$beta_se)^2,
                    sample_size = df_beta$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]


#Calculate the null R2 (binary trait)
# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
y <- pheno[fam.order, on = c("FID", "IID")]
#y <- na.omit(y)
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste0(paste("PC", 1:10, sep = "", collapse = "+"),"+Age") %>%
  paste0("PrCa~", .) %>%
  as.formula %>%
  glm(., data = y, family=binomial) %>%
  summary
null.r2 <- fmsb::NagelkerkeR2(null.model)

#auto model
# Get adjusted beta from the auto model
print(h2_est)
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = 1),
  ncores = NCORES,
  shrink_corr = 0.95,
  allow_jump_sign = F,
  verbose = FALSE,
  use_MLE = FALSE
)
beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)
print(paste0("beta_auto length", length(beta_auto)))

#obtain model PRS
if(is.null(obj.bigSNP)){
  obj.bigSNP <- snp_attach("~/scratch/UKBB_PrCa/UKBB_Afr_PrCa_imputed.rds")
}
genotype <- obj.bigSNP$genotypes
# calculate PRS for all samples
ind.test <- 1:nrow(genotype)
#genotype <- snp_fastImputeSimple(genotype, ncores = numCores) 
pred_auto <-
  big_prodMat(genotype,
              beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`,
              ncores = NCORES)
# scale the PRS generated from AUTO
pred_scaled <- apply(pred_auto, 2, sd)
final_beta_auto <-
  rowMeans(beta_auto[,
                     abs(pred_scaled -
                           median(pred_scaled)) <
                       3 * mad(pred_scaled)])
pred_auto <-
  big_prodVec(genotype,
              final_beta_auto,
              ind.row = ind.test,
              ind.col = info_snp$`_NUM_ID_`)

#Get the final performance of the LDpred models
reg.formula <- paste0(paste("PC", 1:10, sep = "", collapse = "+"),"+Age") %>%
  paste0("PrCa~PRS+", .) %>%
  as.formula
reg.dat <- y
reg.dat$PRS <- pred_auto
reg.dat <- na.omit(reg.dat)
auto.model <- lm(reg.formula, dat=reg.dat) %>% summary
(result <- data.table(
  auto = auto.model$r.squared - null.r2,
  null = null.r2
))

fwrite(reg.dat, "~/scratch/UKBB_PrCa/UKB_all_PrCa_PRS.txt", sep = "\t")
fwrite(auto.model, "~/scratch/UKBB_PrCa/UKB_all_PrCa_PRS_autoperf.txt", sep = "\t")