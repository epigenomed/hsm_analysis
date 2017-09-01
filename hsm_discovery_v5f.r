## R
## HSM analysis
## March 2016
## v5disc 1 sample each family with blood and other covariates

###------------------------------------------------------------------------------------###     

setwd("/home/hsm_analysis/HSM_results")

###------------------------------------------------------------------------------------###  
# packages

library(Cairo)
library(lme4)
library(car)
library(plyr)

###------------------------------------------------------------------------------------###  
# global

argv <- commandArgs(TRUE)

###------------------------------------------------------------------------------------###
#All Epitwin File
TwoSet.df <- read.table("../data/combined_set1and2_f2.txt.gz", header = T) 
Incl.cols <- read.table("../data/InclCols_from_comb_header.txt")
TwoSetRED.df <- TwoSet.df[TwoSet.df$column %in% Incl.cols[, 1], ]  

 
# Batch
Batch.v <- c()
Batch.v = recode(TwoSetRED.df$Batch, "'Batch10'=1; 'Batch3'=2; 'Batch4'=3; 'Batch5'=4; 
'Batch6'=5; 'Batch8'=6; 'Batch9'=7; 'KCLDepression'=8; 'Longitudinal'=9; 'Pain'=10;
'PainIllumina'=11; 'PainQCPilot'=12; 'PainRA'=13; 'T2D'=14", as.factor.result=TRUE)

# Zygosity
Zygosity <-c()
Zygosity[TwoSetRED.df$ACTUAL_ZYGOSITY == "MZ"] <- as.numeric(0)
Zygosity[TwoSetRED.df$ACTUAL_ZYGOSITY == "DZ"] <- as.numeric(1)
Zygosity_MZ0DZ1 <- as.factor(Zygosity)

# Sex
SEX <- c()
SEX[TwoSetRED.df$SEX == "F"] <- as.numeric(0)
SEX[TwoSetRED.df$SEX == "M"] <- as.numeric(1)
SEX_F0 <- as.factor(SEX)

## Get None NA columns
# Include these Columns
#[4]"KCLfam"             
#[7]"ageDNAextraction"
#[8]"Batch"            
#[11]"SEX"
#[12]"ACTUAL_ZYGOSITY"    
#[13]"CLOSEST_SMOKE_DNA"
#[14]"eosinophils"
#[15]"lymphocytes"      
#[16]"monocytes"
#[17]"neutrophils"


# Need Blood Info
Blood.df <- read.table("../data/FinalEpiKeyMerge_bloodCellCounts_BD_LR_2015_Feb08.txt.gz", header = T, stringsAsFactors = F)

TwoSetRdB.df <- merge(TwoSetRED.df[, c(2:13, 15)], Blood.df[, c(3, 14:19)], by = "BGIid", all.x = T)
TwoSetRdBs.df <- TwoSetRdB.df[order(TwoSetRdB.df$column), ]

Nb1.v <- c(4, 7, 8, 11:17) 
Nb1.r <- which(complete.cases(TwoSetRdBs.df[, Nb1.v]) == TRUE)


###------------------------------------------------------------------------------------###

# v5disc Selection
# Replace Nb1.r from here:

# Female Only
Fem.r <- which(TwoSetRdBs.df$SEX == "F")
NInit.r <- which(!duplicated(TwoSetRdBs.df$KCLfam))
Nfem.r <- intersect(Nb1.r, Fem.r)
NDisc.r <- intersect(Nfem.r, NInit.r)

###------------------------------------------------------------------------------------###

#Define HSM function

##HSMfx1 <- function(z) {
##  					fnull.v <- lmer(as.numeric(x)~1 + as.numeric(TwoSetRED.df$ageDNAextraction[Nd2.r]) + as.factor(TwoSetRED.df$CLOSEST_SMOKE_DNA[Nd2.r]) + Batch.v[Nd2.r] + as.numeric(TwoSetRdBs.df$lymphocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$monocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$neutrophils)[Nd2.r] + as.numeric(TwoSetRdBs.df$eosinophils)[Nd2.r] + (1|as.factor(TwoSetRdBs.df$KCLfam[Nd2.r])) + (1|as.factor(Zygosity_MZ0DZ1[Nd2.r])),REML=FALSE)
##                      fit.v   <- lmer(as.numeric(x)~1 + as.numeric(TSGenoS.df[, SNPnow.v][Nd2.r]) + as.numeric(TwoSetRED.df$ageDNAextraction[Nd2.r]) + as.factor(TwoSetRED.df$CLOSEST_SMOKE_DNA[Nd2.r]) + Batch.v[Nd2.r] + as.numeric(TwoSetRdBs.df$lymphocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$monocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$neutrophils)[Nd2.r] + as.numeric(TwoSetRdBs.df$eosinophils)[Nd2.r] + (1|as.factor(TwoSetRdBs.df$KCLfam[Nd2.r])) + (1|as.factor(Zygosity_MZ0DZ1[Nd2.r])),REML=FALSE)
##                      av.v = anova(fit.v, fnull.v)
##                      s = summary(fit.v)
##                      return(c(coef(s)[2,'Estimate'], coef(s)[2,'Std. Error'], av.v$Pr[2]))
##                      }

#F Only
HSMav5disc_noS <- function(z) {
                              fnull.v <- lm(qqnorm(z[Nd2.r], plot.it = FALSE)$x ~1 + as.numeric(TwoSetRdBs.df$ageDNAextraction[Nd2.r]) + as.numeric(TwoSetRdBs.df$lymphocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$monocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$neutrophils)[Nd2.r] + as.numeric(TwoSetRdBs.df$eosinophils)[Nd2.r] + as.factor(TwoSetRdBs.df$CLOSEST_SMOKE_DNA[Nd2.r]) + Batch.v[Nd2.r]))
                              fit.v   <- lm(qqnorm(z[Nd2.r], plot.it = FALSE)$x ~1 + as.numeric(TSGenoS.df[, SNPnow.v][Nd2.r]) + as.numeric((TwoSetRdBs.df$ageDNAextraction[Nd2.r])) + as.numeric(TwoSetRdBs.df$lymphocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$monocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$neutrophils)[Nd2.r] + as.numeric(TwoSetRdBs.df$eosinophils)[Nd2.r]+ as.factor(TwoSetRdBs.df$CLOSEST_SMOKE_DNA[Nd2.r]) + Batch.v[Nd2.r]))
                              av.v    <- anova(fit.v, fnull.v) 
                              oldskS.v <- summary(fit.v)$coeff[2,1]
                              return(c(av.v$Pr[2], oldskS.v))
                              # below need to t
                              }
                              
HSMav6_noS <- function(z) {
                          fnull.v <- lmer(qqnorm(z[Nd2.r], plot.it = FALSE)$x ~1 + as.factor(TwoSetRdBs.df$CLOSEST_SMOKE_DNA[Nd2.r]) + as.factor(Batch.v[Nd2.r]) + as.numeric(TwoSetRdBs.df$ageDNAextraction[Nd2.r]) + as.numeric(TwoSetRdBs.df$lymphocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$monocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$neutrophils)[Nd2.r] + as.numeric(TwoSetRdBs.df$eosinophils)[Nd2.r] + (1|as.factor(Zygosity_MZ0DZ1[Nd2.r])) + (1|as.factor(TwoSetRdBs.df$KCLfam[Nd2.r])), REML=FALSE)
                          fit.v   <- lmer(qqnorm(z[Nd2.r], plot.it = FALSE)$x ~1 + as.numeric(TSGenoS.df[, SNPnow.v][Nd2.r]) + as.factor(TwoSetRdBs.df$CLOSEST_SMOKE_DNA[Nd2.r]) + as.factor(Batch.v[Nd2.r]) + as.numeric(TwoSetRdBs.df$ageDNAextraction[Nd2.r]) + as.numeric(TwoSetRdBs.df$lymphocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$monocytes)[Nd2.r] + as.numeric(TwoSetRdBs.df$neutrophils)[Nd2.r] + as.numeric(TwoSetRdBs.df$eosinophils)[Nd2.r] + (1|as.factor(Zygosity_MZ0DZ1[Nd2.r])) + (1|as.factor(TwoSetRdBs.df$KCLfam[Nd2.r])), REML=FALSE)
                          av.v    <- anova(fit.v, fnull.v)
                          oldsk.v <- lm(TSGenoS.df[, SNPnow.v][Nd2.r] ~ qqnorm(z[Nd2.r], plot.it = FALSE)$x)
                          oldskS.v <- summary(oldsk.v)$coeff[2,1]
                          return(c(av.v$Pr[2], oldskS.v))
                          # below need to t
                          }
                                                    
HSMreplic7 <- function(z) {
              fnull.v <- lmer(qqnorm(z[SUxNd2.r], plot.it = FALSE)$x ~1 + as.numeric(TwoSetRdBs.df$ageDNAextraction[SUxNd2.r]) + as.factor(Batch.v[SUxNd2.r]) + (1|as.factor(Zygosity_MZ0DZ1[SUxNd2.r])) + (1|as.factor(TwoSetRdBs.df$KCLfam[SUxNd2.r])) + as.factor(SEX_F0[SUxNd2.r]), REML=FALSE)
              fit.v   <- lmer(qqnorm(z[SUxNd2.r], plot.it = FALSE)$x ~1 + as.numeric(TSGenoS.df[, SNPnow.v][SUxNd2.r]) + as.numeric(TwoSetRdBs.df$ageDNAextraction[SUxNd2.r]) + as.factor(Batch.v[SUxNd2.r]) + (1|as.factor(Zygosity_MZ0DZ1[SUxNd2.r])) + (1|as.factor(TwoSetRdBs.df$KCLfam[SUxNd2.r])) + as.factor(SEX_F0[SUxNd2.r]), REML=FALSE)
              av.v    <- anova(fit.v, fnull.v)
              oldsk.v <- lm(TSGenoS.df[, SNPnow.v][SUxNd2.r] ~ qqnorm(z[SUxNd2.r], plot.it = FALSE)$x)
              oldskS.v <- summary(oldsk.v)$coeff[2,1]
              return(c(av.v$Pr[2], oldskS.v))
              # below need to t
              }
                 
HSMall <- function(z) {
                      oldsk.v <- lm((qqnorm(z[SUxNd2.r], plot.it = FALSE)$x) ~ (as.numeric(TSGenoS.df[, SNPnow.v][SUxNd2.r])))
                      oldskP.v <- summary(oldsk.v)$coeff[2,4]
                      oldskS.v <- summary(oldsk.v)$coeff[2,1]
                      return(c(oldskP.v, oldskS.v))
                      # below need to t
                      } 
                      
                      
###------------------------------------------------------------------------------------###
##
## By Chrom  

chr.v <- as.numeric(argv[1])

  SNPs.name <- paste("../SNP_genotypes/chr", chr.v, ".raw", sep = "")

  MeDIP.path <- paste("../MeDIP_LD_blocks/chr", chr.v, sep = "")
  
  File.name <- paste("../SNP_LD_blocks/chr", chr.v, "/chr", chr.v, "_GWAS_SNP_LD_bedfiles_all.txt", sep = "")
  
  # PhenoInfo
  Pheno.name <- paste("../SNP_LD_blocks/chr", chr.v,  "/chr", chr.v, "_GWAS_SNP_pheno.bed", sep ="") 
  Pheno.info <- read.table(Pheno.name, header = F, stringsAsFactors = F, sep = "\t")

  # Get All Chrom Genotype Info
  SNPinfo.df <- read.table(SNPs.name, header = T)
  colnames(SNPinfo.df)[2] <- c("KCLid")

  TSGeno.df <- merge(TwoSetRED.df, SNPinfo.df[, c(2, 7:(length(colnames(SNPinfo.df))))], by = "KCLid", all.x = T)
  TSGenoS.df <- TSGeno.df[order(TSGeno.df$column), ]

  # Check MeDIP header info
  # MeDIP.head <- paste("/home/DTR/EPITWIN/Analysis/combined_MeDIP/reduced_rpm_bed/chr", chr.v, "_Combined_RED_rpm_bedheader.txt", sep = "")
  # MeDIP.header <- read.table(MeDIP.head)
  # as.numeric(MeDIP.header[1, 4:4353]) - as.numeric(TSGenoS.df$BGI) #Should All Be Zero

  Files.v <- read.table(File.name, header = F, stringsAsFactors = F) 
  Files.l <- length(Files.v[, 1])
  TopHit.df <- data.frame(matrix(nrow = (Files.l), ncol = 6))

  # LD info - get correct LD block that matches to SNP
  # chr*_GWAS_SNP_LD_info.txt
  # col4 = SNP, col5 = LD equivalent SNP
  LDinfo.name <- paste("../SNP_LD_blocks/chr", chr.v,  "/chr", chr.v, "_GWAS_SNP_LD_info.txt", sep ="")
  LDinfo.df <-  read.table(LDinfo.name, sep = "\t", header = F)
 
  # SNP info
  SNPinfo.head <- colnames(TSGenoS.df)
  SNPinfoR.head <- unlist(strsplit(SNPinfo.head[16:length(SNPinfo.head)], "_"))
  oddie.v <- seq(1, (length(SNPinfoR.head)), by = 2)
  SNPinfoR.in <- SNPinfoR.head[oddie.v] 
  
  # Complete SNP count
  CmplSNP.df <- data.frame(matrix(nrow = Files.l, ncol = 2))
   
  # Mean LD block info
  LDbkMean.df <- data.frame(matrix(nrow = Files.l, ncol = 2))
   
  for(zz in 1:(Files.l))

    {
    chrSNP.v <- strsplit (Files.v[zz, 1], ".bed")[[1]][1]
    SNP.v <- strsplit(chrSNP.v, "_")[[1]][2] 
    Gene.v <- Pheno.info[Pheno.info$V4 == SNP.v, ][, 13]
    Pheno.v <-  Pheno.info[Pheno.info$V4 == SNP.v, ][, 8]
    SubT.v <- paste(Gene.v, ";", Pheno.v)
  
    #Get right SNP
    SNPnow.v <- (which(SNPinfoR.in == SNP.v)) + 15
  
    #Count Number of Missing Genotypes
    CmplSNP.df[zz, 1] <- SNP.v
    #Get Only non NA genotype rows
    ComlSNP.r <- which(!(is.na(TSGenoS.df[, SNPnow.v])))
    
    #v5disc version
    Nd2.r <- intersect(ComlSNP.r, NDisc.r) 
    CmplSNP.df[zz, 2] <- length(Nd2.r)
    
    #Check if monomorphic
    if(length(table(TSGenoS.df[, SNPnow.v][Nd2.r])) != 3){next}
    
    #Get correct LD block    
    MeSNPLD.name <- LDinfo.df$V5[which(LDinfo.df$V4 == unlist(strsplit(unlist(strsplit(Files.v[zz, 1], ".bed")), paste("chr", chr.v, "_", sep ="")))[2])]
    MeSNP.get <- which(Files.v[, 1] == paste("chr", chr.v, "_", MeSNPLD.name, ".bed", sep =""))[1]
        
    MeDIP.name  <- paste(MeDIP.path, "/", Files.v[MeSNP.get, 1], "_LD_rpm.bed.gz", sep = "")
    #Check if empty
    if(file.info(MeDIP.name)[1, 1] < 100){next}
     
    MeDIP.bed <- read.table(MeDIP.name, sep = "\t", header = F)  
   
    #v5disc version F ONLY  
    rHSM.v <- t(apply(MeDIP.bed[, 4:4353], 1, HSMav5disc_noS))
    
    # plus want null p value
    rHSMgpt.v <- paste(chrSNP.v, "_methyl",  sep = "")
    rHSMout.v <- paste("v5disc/chr", chr.v, "/", chrSNP.v, "_hsm_results_v5disc.txt", sep = "")
    
    rHSMout.df <- cbind(MeDIP.bed[, 1:3], rHSM.v)
    #colnames(rHSMout.df0[4] <- c("hsm_p")

    write.table(rHSMout.df, rHSMout.v, sep = '\t', quote = F, row.names = F, col.names = F)
          
    }
  
  CmplSNP.name <- paste("v5disc/chr", chr.v, "/chr", chr.v, "_complete_intersect_number_vAnova5disc.txt", sep = "")
  
  write.table(CmplSNP.df, CmplSNP.name, sep = '\t', quote = F, row.names = F, col.names = F)


###------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------###

