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

numCores <- detectCores()
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
# Open a temporary file
tmp <- tempfile(tmpdir = "~/Analysis/PRS/tmp-data")
tmp1 <- tempfile(tmpdir = "~/Analysis/PRS/tmp-data1")
on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
on.exit(file.remove(paste0(tmp1, ".sbk")), add = TRUE)
# Initialize variables for storing the LD score and LD matrix
corr <- NULL
ld <- NULL
# We want to know the ordering of samples in the bed file 
info_snp <- NULL
fam.order <- NULL
corr_list = list()

for (chr in c(6,8)) { #1:22
  print(paste0("working on chromosome", chr))
  gc()
  # preprocess the bed file (only need to do once for each data set)
  # Assuming the file naming is EUR_chr#.bed
  snp_readBed(paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed_",chr,".bed")) #, backingfile = paste0("/Volumes/projects/TEMP_joel/PRS/tmp-data", chr)
  # now attach the genotype object
  obj.bigSNP <- snp_attach(paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed_",chr,".rds"))
  # extract the SNP information from the genotype
  map <- obj.bigSNP$map[-3]
  names(map) <- c("chr", "rsid", "pos", "a1", "a0")
  # perform SNP matching
  # paste0("Length of extended sumstats is ", nrow(sumstats[sumstats$chr==chr,-which(names(sumstats) %in% "rsid")]))
  tmp_snp <- snp_match(sumstats[sumstats$chr==chr, ], map, match.min.prop = 0)
  info_snp <- rbind(info_snp, tmp_snp)
  # Assign the genotype to a variable for easier downstream analysis
  genotype <- obj.bigSNP$genotypes
  # Rename the data structures
  CHR <- map$chr
  POS <- map$pos
  # get the CM information from 1000 Genome
  # will download the 1000G file to the current directory (".")
  POS2 <- snp_asGeneticPos(CHR, POS, dir = "~/Analysis/GWASData/")
  # calculate LD
  # Extract SNPs that are included in the chromosome
  ind.chr <- which(tmp_snp$chr == chr)
  ind.chr2 <- tmp_snp$`_NUM_ID_`[ind.chr]
  # Calculate the LD
  corr0 <- snp_cor(
    genotype,
    ind.col = ind.chr2,
    ncores = numCores,
    infos.pos = POS2[ind.chr2],
    size = 3 / 1000
  )
  if (chr == 6) { #change to 1
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp)
    # print(ld)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
    # print(ld)
  }
  # We assume the fam order is the same across different chromosomes
  if(is.null(fam.order)){
    fam.order <- as.data.table(obj.bigSNP$fam)
  }
}
# Rename fam order
setnames(fam.order,
         c("family.ID", "sample.ID"),
         c("FID", "IID"))

#Perform LD score regression
df_beta <- info_snp[,c("beta", "beta_se", "n_eff", "_NUM_ID_")]
ld_NA <- which(is.na(ld))
print(paste0("ld length before NAs = ", length(ld)))
ld <- ld[!is.na(ld)]
df_beta_subset <- df_beta[-ld_NA, ]
print(paste0("ld length after NAs = ", length(ld)))
print(paste0("no. of rows in beta = ", nrow(df_beta_subset)))
ldsc <- snp_ldsc(   ld, 
                    length(ld), 
                    chi2 = (df_beta_subset$beta / df_beta_subset$beta_se)^2,
                    sample_size = df_beta_subset$n_eff, 
                    blocks = NULL)
h2_est <- ldsc[["h2"]]
print(paste0("heritability = ", h2_est))

#removing ld NAs from corr
corr <- as_SFBM(corr[-ld_NA, -ld_NA], tmp1) #added extra
#remove NAs from info_snp
info_snp <- info_snp[-ld_NA, ]

#Calculate the null R2 (binary trait)
# Reformat the phenotype file such that y is of the same order as the 
# sample ordering in the genotype file
y <- pheno[fam.order, on = c("FID", "IID")]
#y <- na.omit(y)
# Calculate the null R2
# use glm for binary trait 
# (will also need the fmsb package to calculate the pseudo R2)
null.model <- paste0(paste("PC", 1:10, sep = "", collapse = "+")) %>%
  paste0("PrCa~", .) %>%
  as.formula %>%
  glm(., data = y, family=binomial) 
null.r2 <- fmsb::NagelkerkeR2(null.model)

#auto model
# Get adjusted beta from the auto model
multi_auto <- snp_ldpred2_auto(
  corr,
  df_beta_subset,
  h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.9, length.out = 1),
  ncores = 9,
  shrink_corr = 0.95,
  allow_jump_sign = F,
  verbose = F,
  use_MLE = F
)
beta_auto <- sapply(multi_auto, function(auto)
  auto$beta_est)
print(length(beta_auto))

#obtain model PRS
pred_auto <- NULL
for(chr in c(6,8)){ #1:22
  obj.bigSNP <- snp_attach(paste0("~/Analysis/UKBB/UKBB_Afr_PrCa_imputed_",chr,".rds"))
  genotype <- obj.bigSNP$genotypes
  # calculate PRS for all samples
  ind.test <- 1:nrow(genotype)
  # Extract SNPs in this chromosome
  chr.idx <- which(info_snp$chr == chr)
  ind.chr <- info_snp$`_NUM_ID_`[chr.idx]
  beta_auto_subset <- matrix(beta_auto[chr.idx,])
  print(nrow(ind.chr))
  #imputation
  #genotype <- snp_fastImputeSimple(genotype, ncores = numCores) #added new
  print(length(genotype))
  tmp <-
    big_prodVec(genotype,
                beta_auto_subset,
                ind.row = ind.test,
                ind.col = ind.chr)
  if(is.null(pred_auto)){
    pred_auto <- tmp
  }else{
    pred_auto <- pred_auto + tmp
  }
}

#Get the final performance of the LDpred models
reg.formula <- paste0(paste("PC", 1:10, sep = "", collapse = "+")) %>%
  paste0("PrCa~PRS+", .) %>%
  as.formula
reg.dat <- y
reg.dat$PRS <- pred_auto
reg.dat <- na.omit(reg.dat)
auto.model <- lm(reg.formula, dat=reg.dat) %>% summary
(result <- data.table(
  auto = auto.model$r.squared - null.r2$R2,
  null = null.r2$N
))

fwrite(reg.dat, "~/Analysis/PRS/UKB_Afr_PrCa_PRS.txt", sep = "\t")
fwrite(auto.model, "~/Analysis/PRS/UKB_Afr_PrCa_autoperf.txt", sep = "\t")