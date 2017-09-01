## R
## ReGraph hsm LD Block Results 
## With additional Colours for Genetic Effects

# *NEED to have BEDtools loaded in bash [bedtools/2.17.0]

###------------------------------------------------------------------------------------###     

setwd("/home/hsm_analysis/HSM_results")

###------------------------------------------------------------------------------------###  
# packages

library(Cairo)
library("plyr")
library("RColorBrewer") 
library(ggplot2)


###------------------------------------------------------------------------------------###  
# global

argv <- commandArgs(TRUE)

###------------------------------------------------------------------------------------###  

## Colour Define
## Red = hsm original
## Green = CNV
## Blue = Indel
## Dk Orange = STR
## Lt Orange = eSTR
## Black = Blacklist Regions
## Brown = tagging SNP 

mypalette<-brewer.pal(12,"Paired")


## colz.v
colz.v <- c()
colz.v[1] <- mypalette[6] 			#novariant	red
colz.v[2] <- c("black") 			#black		black
colz.v[3] <- mypalette[10]			#multiple	purple
colz.v[4] <- mypalette[2]			#indel		dk blue
colz.v[5] <- mypalette[4]	 		#cnv		green
colz.v[6] <- mypalette[8]			#str		dk orange
colz.v[7] <- mypalette[7]			#estr		lt orange
#colz.v[8] <- mypalette[12]			#snp name	brown

##--------------------------------------------------------------------------------------###

##IntersectBed Wrapper

#bedtools wrapper
bedTools.2in <- function(functionstring = "intersectBed", bed1, bed2, opt.string = "-wo -f 0.9")
{
  #create temp files
  a.file = tempfile()
  b.file = tempfile()
  out    = tempfile()
  options(scipen =99) # not to use scientific notation when writing out
 
  #write bed formatted dataframes to tempfile
  write.table(bed1, file = a.file, quote = F, sep = "\t", col.names = F, row.names = F)
  write.table(bed2, file = b.file, quote = F, sep = "\t", col.names = F, row.names = F)
 
  # create the command string and call the command using system()
  command = paste(functionstring, "-a", a.file, "-b", b.file, opt.string, ">", out, sep=" ")
  cat(command, "\n")
  try(system(command))
 
  res = read.table(out, header = F)
  unlink(a.file); unlink(b.file); unlink(out)
  return(res)
} 

##---------------------------------------------------------------------------------------###

# Variant database # Defined in reGraph_hsm_v1.r
# ("black")("multiple")("indel")("cnv")("str")("estr")
# hg19_methyl_wind_CombinedVariant.bed


##---------------------------------------------------------------------------------------###
##---------------------------------------------------------------------------------------###

## BACK TO GW HSM
#All Epitwin File

TwoSet.df <- read.table("../data/combined_set1and2_f2.txt.gz", header = T, stringsAsFactors = F) 
Incl.cols <- read.table("../data/InclCols_from_comb_header.txt")
TwoSetRED.df <- TwoSet.df[TwoSet.df$column %in% Incl.cols[, 1], ]  

# Need Blood Info
Blood.df <- read.table("../data/FinalEpiKeyMerge_bloodCellCounts_BD_LR_2015_Feb08.txt.gz", header = T, stringsAsFactors = F)

TwoSetRdB.df <- merge(TwoSetRED.df[, c(2:13, 15)], Blood.df[, c(3, 14:19)], by = "BGIid", all.x = T)
TwoSetRdBs.df <- TwoSetRdB.df[order(TwoSetRdB.df$column), ]

Nb1.v <- c(4, 7, 8, 11:17) 
Nb1.r <- which(complete.cases(TwoSetRdBs.df[, Nb1.v]) == TRUE)

           
##---------------------------------------------------------------------------------------###
##
## By Chrom  

chr.v <- as.numeric(argv[1])
  
  # Get Variant Bedfile for Chromosome
  ComVar.name <- paste("../methylome_windows/combinedVariants/chr", chr.v, "_methyl_wind_CombinedVariant.bed.gz", sep = "")
  CombVar.df <- read.table(ComVar.name, header = F, stringsAsFactors = F)
  
  
  ## SNP etc from old
  SNPs.name <- paste("../SNP_genotypes/chr", chr.v, ".raw", sep = "")
  File.name <- paste("../SNP_LD_blocks/chr", chr.v, "/chr", chr.v, "_GWAS_SNP_LD_bedfiles_all.txt", sep = "")

  #PhenoInfo
  Pheno.name <- paste("../SNP_LD_blocks/chr", chr.v,  "/chr", chr.v, "_GWAS_SNP_pheno.bed", sep ="") 
  Pheno.info <- read.table(Pheno.name, header = F, stringsAsFactors = F, sep = "\t")

  # Get All Chrom Genotype Info
  SNPinfo.df <- read.table(SNPs.name, header = T)
  colnames(SNPinfo.df)[2] <- c("KCLid")

  TSGeno.df <- merge(TwoSetRED.df, SNPinfo.df[, c(2, 7:(length(colnames(SNPinfo.df))))], by = "KCLid", all.x = T)
  TSGenoS.df <- TSGeno.df[order(TSGeno.df$column), ]

  Files.v <- read.table(File.name, header = F, stringsAsFactors = F) 
  Files.l <- length(Files.v[, 1])
 
  # SNP info
  SNPinfo.head <- colnames(TSGenoS.df)
  SNPinfoR.head <- unlist(strsplit(SNPinfo.head[16:length(SNPinfo.head)], "_"))
  oddie.v <- seq(1, (length(SNPinfoR.head)), by = 2)
  SNPinfoR.in <- SNPinfoR.head[oddie.v] 
     
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
     ComlSNP.r <- which(!(is.na(TSGenoS.df[, SNPnow.v])))
    
    # Which to Include so no NA for anova
    Nd2.r <- intersect(ComlSNP.r, Nb1.r) 
    
    #Check if monomorphic
    if(length(table(TSGenoS.df[, SNPnow.v][Nd2.r])) != 3){next}
    
    #rHSMgph.v <- paste("vAll/chr", chr.v, "/", chrSNP.v, "_methyl_vAll.pdf", sep = "") 
    rHSMgph2.v <- paste("vAll/chr", chr.v, "/reGraphAll/", chrSNP.v, "_methyl_vAll_gg.pdf", sep = "") 
    rHSMgph3.v <- paste("vAll/chr", chr.v, "/reGraphAll/", chrSNP.v, "_methyl_vAll_gg2.pdf", sep = "") 
  
    rHSMgpt.v <- paste(chrSNP.v, "_methyl",  sep = "")
    
    rHSMin.v <- paste("vAll/chr", chr.v, "/", chrSNP.v, "_hsm_results_vAll.txt", sep = "")

    ##---***---***---***----------------------------------------------------------------------##

    # Get the Results File
    rHSMin.df <- read.table(rHSMin.v, header = F, stringsAsFactors = F) 
      
    ##---***---***---***----------------------------------------------------------------------##

    ## Graph  
 
    ##----------------------------------------------------------------------------------------##

    # Intersect CombVar.df [is chr specific file]
    # bedTools.2in("bedIntersect", bed1, bed2 -wo -f 0.9) #== "intersectBed -a bed1 -b bed2 -wo - f 0.9"
    hsm.df <- bedTools.2in("intersectBed", rHSMin.df, CombVar.df)[, c(1:4, 10)]

    hsm.df$"V10" <- as.character(hsm.df$"V10")
    colnames(hsm.df) <- c("chr", "start", "stop", "p", "variant")
    
    # Ceiling at p-300 and replace 0
    hsm.df$"p"[hsm.df$"p" < 1E-300] <- (1E-300)

    ##----------------------------------------------------------------------------------------##
    # Define Variant Continuous Groupings #
    
    blackset.v <- which(hsm.df$variant == "black") 
    if(length(blackset.v) != 0){blackset.list <- unname(tapply(blackset.v, cumsum(c(1, diff(blackset.v)) != 1), range))} 

    multset.v <- which(hsm.df$variant == "multiple") 
    if(length(multset.v) != 0){multset.list <- unname(tapply(multset.v, cumsum(c(1, diff(multset.v)) != 1), range))} 

    InDelset.v <- which(hsm.df$variant == "indel") 
    if(length(InDelset.v) != 0){InDelset.list <- unname(tapply(InDelset.v, cumsum(c(1, diff(InDelset.v)) != 1), range))} 

    CNVset.v <- which(hsm.df$variant == "cnv") 
    if(length(CNVset.v) != 0){CNVset.list <- unname(tapply(CNVset.v, cumsum(c(1, diff(CNVset.v)) != 1), range))} 

    STRset.v <- which(hsm.df$variant == "str") 
    if(length(STRset.v) != 0){STRset.list <- unname(tapply(STRset.v, cumsum(c(1, diff(STRset.v)) != 1), range))} 
 
    eSTRset.v <- which(hsm.df$variant == "estr") 
    if(length(eSTRset.v) != 0){eSTRset.list <- unname(tapply(eSTRset.v, cumsum(c(1, diff(eSTRset.v)) != 1), range))} 

   ##----------------------------------------------------------------------------------------##
    # Graphing #
    
    pic <- c()
    
    pic <- ggplot(data = hsm.df, aes(y = -log10(p), x = ((start + 250)/1000000))) 
 
    pic <- pic + geom_area(fill = colz.v[1])
    pic <- pic + labs(list(x = "position (Mb)", y = expression(-log[10](italic(p)))))
    pic <- pic + ggtitle(bquote(atop(.(SNP.v), atop(italic(.(SubT.v[1])), "")))) 
    
    ## STR subset 
    if(length(STRset.v) != 0){
      pp <- c() 
      for(pp in 1:length(STRset.list)){
        begin.v <- hsm.df$start[STRset.list[[pp]][1]]
        end.v <- hsm.df$start[STRset.list[[pp]][2]]
        pic <- pic + geom_area(data = subset(hsm.df, start >=  begin.v & start <= end.v), fill = colz.v[6])
        }  
      }

    ## eSTR subset
    if(length(eSTRset.v) != 0){
      pp <- c() 
      for(pp in 1:length(eSTRset.list)){
        begin.v <- hsm.df$start[eSTRset.list[[pp]][1]]
        end.v <- hsm.df$start[eSTRset.list[[pp]][2]]
        pic <- pic + geom_area(data = subset(hsm.df, start >=  begin.v & start <= end.v), fill = colz.v[7])
        }  
      }
   
    ## CNV subset 
    if(length(CNVset.v) != 0){
      pp <- c() 
      for(pp in 1:length(CNVset.list)){
        begin.v <- hsm.df$start[CNVset.list[[pp]][1]]
        end.v <- hsm.df$start[CNVset.list[[pp]][2]]
        pic <- pic + geom_area(data = subset(hsm.df, start >=  begin.v & start <= end.v), fill = colz.v[5])  
        }  
      }
   
    ## InDel subset 
    if(length(InDelset.v) != 0){
      pp <- c() 
      for(pp in 1:length(InDelset.list)){
        begin.v <- hsm.df$start[InDelset.list[[pp]][1]]
        end.v <- hsm.df$start[InDelset.list[[pp]][2]]
        pic <- pic + geom_area(data = subset(hsm.df, start >=  begin.v & start <= end.v), fill = colz.v[4])  
        }  
      }
 
    ## mult subset 

    if(length(multset.v) != 0){
      pp <- c() 
      for(pp in 1:length(multset.list)){
        begin.v <- hsm.df$start[multset.list[[pp]][1]]
        end.v <- hsm.df$start[multset.list[[pp]][2]]
        pic <- pic + geom_area(data = subset(hsm.df, start >=  begin.v & start <= end.v), fill = colz.v[3])  
        }  
      }
 
    ## black subset 
    if(length(blackset.v) != 0){
      pp <- c() 
      for(pp in 1:length(blackset.list)){
        begin.v <- hsm.df$start[blackset.list[[pp]][1]]
        end.v <- hsm.df$start[blackset.list[[pp]][2]]
        pic <- pic + geom_area(data = subset(hsm.df, start >=  begin.v & start <= end.v), fill = colz.v[2])  
        }  
      }

    
	#SNP location
	SNPloc.v <- Pheno.info$V3[Pheno.info$V4 == SNP.v][1]
    pic <- pic + annotate("segment", x = (SNPloc.v/(1E+6)), xend = (SNPloc.v/(1E+6)), y = 0, yend = (max(-log10(hsm.df$p))/10), colour = "black", linetype = 2)
    pic <- pic + annotate("text", x = SNPloc.v/(1E+6), y = max(-log10(hsm.df$p))/8, label = SNP.v, size = 4, colour = "black")
    pic <- pic + annotate("point", x = SNPloc.v/(1E+6), y = 0, colour = "black")
 
    # CairoPDF(file = rHSMgph2.v, width = 16, height = 9)
    # Standard width = 9, height = 6) 
    CairoPDF(file = rHSMgph2.v, width = 9, height = 6)
    print(pic)
    dev.off()
    
    # theme_bw(base_size=20)
     
    # More Landscape
    CairoPDF(file = rHSMgph3.v, width = 16, height = 9) 
    print(pic)
    dev.off()
    
  }
  
##---***---***---***--------------------------------------------------------------------##
##---***---***---***--------------------------------------------------------------------##

