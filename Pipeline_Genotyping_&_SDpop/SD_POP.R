library(dplyr)
library(qqman)
library(ggplot2)

#####################################
##### 1. Chargement GENE TABLE ######
#####################################


GENE_TABLE_path <- "/Users/cibio/Desktop/Laurus/SDPOP/braker_fixed_ids_gene_gtf"
GENE_TABLE <- read.table(GENE_TABLE_path,header = FALSE,sep = "\t",fill = TRUE,comment.char = "",quote = "",col.names = c("scaffold", "feature", "start", "end", "gene_id"))

# 1. Numeroter les contigs
scaff_order <- c(
  unique(GENE_TABLE$scaffold[grepl("^OZ", GENE_TABLE$scaffold)]),
  unique(GENE_TABLE$scaffold[!grepl("^OZ", GENE_TABLE$scaffold)]))
GENE_TABLE$CHR_NUM <- match(GENE_TABLE$scaffold, scaff_order)

# 2. Clean table
GENE_TABLE$feature <- NULL
GENE_TABLE <- GENE_TABLE[, c("scaffold", "CHR_NUM", "gene_id","start", "end")]

#########################################
##### 2. Chargement SD_POP results ######
##########################################

SD_POP_nobilis_POR_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_3/SDPOP_ZW_L_nobilis_POR_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_POR_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_3/SDPOP_XY_L_nobilis_POR_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_POR_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_3/SDPOP_nosexchr_L_nobilis_POR_Lazorica_norm_filtered_EXONS.recode.out.contig"

SD_POP_nobilis_TUN_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_6/SDPOP_ZW_L_nobilis_TUN_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_TUN_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_6/SDPOP_XY_L_nobilis_TUN_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_TUN_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_6/SDPOP_nosexchr_L_nobilis_TUN_Lazorica_norm_filtered_EXONS.recode.out.contig"

SD_POP_nobilis_GRE_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_7/SDPOP_ZW_L_nobilis_GRE_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_GRE_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_7/SDPOP_XY_L_nobilis_GRE_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_GRE_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_7/SDPOP_nosexchr_L_nobilis_GRE_Lazorica_norm_filtered_EXONS.recode.out.contig"

SD_POP_nobilis_all_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_4/SDPOP_ZW_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_all_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_4/SDPOP_XY_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_all_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_4/SDPOP_nosexchr_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig"

SD_POP_nobilis_IT_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_5/SDPOP_ZW_L_nobilis_ITA_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_IT_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_5/SDPOP_XY_L_nobilis_ITA_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_IT_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_5/SDPOP_nosexchr_L_nobilis_ITA_Lazorica_norm_filtered_EXONS.recode.out.contig"

SD_POP_azorica_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_azorica_L_azorica/SD_POP/SDPOP_ZW_L_azorica_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_azorica_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_azorica_L_azorica/SD_POP/SDPOP_XY_L_azorica_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_azorica_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_azorica_L_azorica/SD_POP/SDPOP_nosexchr_L_azorica_Lazorica_norm_filtered_EXONS.recode.out.contig"

SD_POP_novocanariensis_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_novocanariensis_L_azorica/SDPOP_resex/SDPOP_ZW_L_novocanariensis_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_novocanariensis_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_novocanariensis_L_azorica/SDPOP_resex/SDPOP_XY_L_novocanariensis_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_novocanariensis_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_novocanariensis_L_azorica/SDPOP_resex/SDPOP_nosexchr_L_novocanariensis_Lazorica_norm_filtered_EXONS.recode.out.contig"

###
SD_POP_nobilis_all_12CHR_ZW_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica_12CHR/SDPOP/SDPOP_ZW_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_all_12CHR_XY_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica_12CHR/SDPOP/SDPOP_XY_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig"
SD_POP_nobilis_all_12CHR_nosexchr_path <- "/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica_12CHR/SDPOP/SDPOP_nosexchr_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig"


read_sdpop <- function(path, colname="gene"){
  df <- read.table(path, header=TRUE, sep="\t", fill=TRUE, comment.char="", quote="")
  df[,1] <- sub("^>", "", df[,1])
  names(df)[1] <- colname
  df
}



#NOBILIS_ALL
#SD_POP_ZW_nobilis_all  <- read_sdpop(SD_POP_nobilis_all_ZW_path)
SD_POP_XY_nobilis_all  <- read_sdpop(SD_POP_nobilis_all_XY_path)
#SD_POP_no_nobilis_all  <- read_sdpop(SD_POP_nobilis_all_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_nobilis_all  <- merge(SD_POP_ZW_nobilis_all, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_nobilis_all  <- merge(SD_POP_XY_nobilis_all, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_nobilis_all  <- merge(SD_POP_no_nobilis_all, GENE_TABLE, by.x=1, by.y=3)

#NOBILIS_POR
#SD_POP_ZW_nobilis_POR  <- read_sdpop(SD_POP_nobilis_POR_ZW_path)
SD_POP_XY_nobilis_POR  <- read_sdpop(SD_POP_nobilis_POR_XY_path)
#SD_POP_no_nobilis_POR  <- read_sdpop(SD_POP_nobilis_POR_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_nobilis_POR  <- merge(SD_POP_ZW_nobilis_POR, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_nobilis_POR  <- merge(SD_POP_XY_nobilis_POR, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_nobilis_POR  <- merge(SD_POP_no_nobilis_POR, GENE_TABLE, by.x=1, by.y=3)

#NOBILIS_ITA
#SD_POP_ZW_nobilis_IT  <- read_sdpop(SD_POP_nobilis_IT_ZW_path)
SD_POP_XY_nobilis_IT  <- read_sdpop(SD_POP_nobilis_IT_XY_path)
#SD_POP_no_nobilis_IT  <- read_sdpop(SD_POP_nobilis_IT_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_nobilis_IT  <- merge(SD_POP_ZW_nobilis_IT, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_nobilis_IT  <- merge(SD_POP_XY_nobilis_IT, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_nobilis_IT  <- merge(SD_POP_no_nobilis_IT, GENE_TABLE, by.x=1, by.y=3)

#NOBILIS_TUN
#SD_POP_ZW_nobilis_TUN  <- read_sdpop(SD_POP_nobilis_TUN_ZW_path)
SD_POP_XY_nobilis_TUN  <- read_sdpop(SD_POP_nobilis_TUN_XY_path)
#SD_POP_no_nobilis_TUN  <- read_sdpop(SD_POP_nobilis_TUN_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_nobilis_TUN  <- merge(SD_POP_ZW_nobilis_TUN, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_nobilis_TUN  <- merge(SD_POP_XY_nobilis_TUN, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_nobilis_TUN  <- merge(SD_POP_no_nobilis_TUN, GENE_TABLE, by.x=1, by.y=3)

#NOBILIS_GRE
#SD_POP_ZW_nobilis_GRE  <- read_sdpop(SD_POP_nobilis_GRE_ZW_path)
SD_POP_XY_nobilis_GRE  <- read_sdpop(SD_POP_nobilis_GRE_XY_path)
#SD_POP_no_nobilis_GRE  <- read_sdpop(SD_POP_nobilis_GRE_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_nobilis_GRE  <- merge(SD_POP_ZW_nobilis_GRE, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_nobilis_GRE  <- merge(SD_POP_XY_nobilis_GRE, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_nobilis_GRE  <- merge(SD_POP_no_nobilis_GRE, GENE_TABLE, by.x=1, by.y=3)


#AZORICA
#SD_POP_ZW_azorica  <- read_sdpop(SD_POP_azorica_ZW_path)
SD_POP_XY_azorica  <- read_sdpop(SD_POP_azorica_XY_path)
#SD_POP_no_azorica  <- read_sdpop(SD_POP_azorica_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_azorica  <- merge(SD_POP_ZW_azorica, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_azorica  <- merge(SD_POP_XY_azorica, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_azorica  <- merge(SD_POP_no_azorica, GENE_TABLE, by.x=1, by.y=3)

#NOVOCANARIENSIS
#SD_POP_ZW_novocanariensis  <- read_sdpop(SD_POP_novocanariensis_ZW_path)
SD_POP_XY_novocanariensis  <- read_sdpop(SD_POP_novocanariensis_XY_path)
#SD_POP_no_novocanariensis  <- read_sdpop(SD_POP_novocanariensis_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_novocanariensis  <- merge(SD_POP_ZW_novocanariensis, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_novocanariensis  <- merge(SD_POP_XY_novocanariensis, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_novocanariensis  <- merge(SD_POP_no_novocanariensis, GENE_TABLE, by.x=1, by.y=3)

#NOBILIS 12 CHR
#SD_POP_ZW_nobilis_all_12CHR  <- read_sdpop(SD_POP_nobilis_all_12CHR_ZW_path)
SD_POP_XY_nobilis_all_12CHR  <- read_sdpop(SD_POP_nobilis_all_12CHR_XY_path)
#SD_POP_no_nobilis_all_12CHR  <- read_sdpop(SD_POP_nobilis_all_12CHR_nosexchr_path)
#On merge le gene table pour avoir les infos des contigs
#SD_POP_ZW_nobilis_all_12CHR  <- merge(SD_POP_ZW_nobilis_all_12CHR, GENE_TABLE, by.x=1, by.y=3)
SD_POP_XY_nobilis_all_12CHR  <- merge(SD_POP_XY_nobilis_all_12CHR, GENE_TABLE, by.x=1, by.y=3)
#SD_POP_no_nobilis_all_12CHR  <- merge(SD_POP_no_nobilis_all_12CHR, GENE_TABLE, by.x=1, by.y=3)



######################################################
##### 3. Manhatan plot XY & AUTOSOMAL GENOME WIDE ####
######################################################

#Manhatan plot NOBILIS
par(mfrow = c(5, 2))
manhattan(
  subset(SD_POP_XY_nobilis_POR, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     
  logp = FALSE,     
  main = "Manhattan plot – XY nobilis POR",
  ylab = "xy_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_POR, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",    
  logp = FALSE,   
  main = "Manhattan plot – autosomal nobilis POR",
  ylab = "autosmal_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

#Manhatan plot NOBILIS IT
par(mfrow = c(2, 1))
manhattan(
  subset(SD_POP_XY_nobilis_IT, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     
  logp = FALSE,     
  main = "Manhattan plot – XY nobilis ITA",
  ylab = "xy_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_IT, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",    
  logp = FALSE,   
  main = "Manhattan plot – autosomal nobilis ITA",
  ylab = "autosmal_prior"
)

#Manhatan plot NOBILIS TUN
#par(mfrow = c(2, 1))
manhattan(
  subset(SD_POP_XY_nobilis_TUN, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     
  logp = FALSE,     
  main = "Manhattan plot – XY nobilis TUN",
  ylab = "xy_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_TUN, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",    
  logp = FALSE,   
  main = "Manhattan plot – autosomal nobilis TUN",
  ylab = "autosmal_prior"
)

#Manhatan plot NOBILIS GRE
#par(mfrow = c(2, 1))
manhattan(
  subset(SD_POP_XY_nobilis_GRE, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     
  logp = FALSE,     
  main = "Manhattan plot – XY nobilis GRE",
  ylab = "xy_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_GRE, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",    
  logp = FALSE,   
  main = "Manhattan plot – autosomal nobilis GRE",
  ylab = "autosmal_prior"
)


# Plot nob | azo | novo


#Manhatan plot NOBILIS ALL
par(mfrow = c(3, 2))
manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     
  logp = FALSE,     
  main = "Manhattan plot – XY nobilis all",
  ylab = "xy_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",    
  logp = FALSE,   
  main = "Manhattan plot – autosomal nobilis all",
  ylab = "autosmal_prior"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
#Manhatan plot AZORICA
manhattan(
  subset(SD_POP_XY_azorica, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     # ou "scaffold.x", "start.x", ou rien
  logp = FALSE,     # max n'est PAS un -log10(p), donc logp=FALSE
  main = "Manhattan plot – XY azorica",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_azorica, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",     # ou "scaffold.x", "start.x", ou rien
  logp = FALSE,     # max n'est PAS un -log10(p), donc logp=FALSE
  main = "Manhattan plot – autosomal azorica",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

#Manhatan plot NOVOCANARIENSIS
manhattan(
  subset(SD_POP_XY_novocanariensis, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",     # ou "scaffold.x", "start.x", ou rien
  logp = FALSE,     # max n'est PAS un -log10(p), donc logp=FALSE
  main = "Manhattan plot – XY novocanariensis",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_novocanariensis, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",     # ou "scaffold.x", "start.x", ou rien
  logp = FALSE,     # max n'est PAS un -log10(p), donc logp=FALSE
  main = "Manhattan plot – autosomal novocanariensis",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
par(mfrow = c(1, 1))


###################################################
##### 4. Manhatan plot XY & AUTOSOMAL CHR6 ####
###################################################
par(mfrow = c(1, 2))
par(mfrow = c(5, 2))
#NOBILIS
manhattan(
  subset(SD_POP_XY_nobilis_all_12CHR, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – XY nobilis POR (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_all_12CHR, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – autosomal nobilis POR (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

#NOBILIS_IT
manhattan(
  subset(SD_POP_XY_nobilis_IT, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – XY nobilis ITA (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_IT, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – autosomal nobilis ITA (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)


par(mfrow = c(3, 2))
#NOBILIS_ALL
manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – XY nobilis all (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – autosomal nobilis all (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

#NOVOCANARIENSIS
manhattan(
  subset(SD_POP_XY_novocanariensis, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – xy novocanariensis (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_novocanariensis, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – autosomal novocanariensis (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

#AZORICA
manhattan(
  subset(SD_POP_XY_azorica, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – xy azorica (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_azorica, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – autosomal azorica (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)


#########################
#### 5. Manhatan plot paralogs CHR 6 #####
#########################


par(mfrow = c(5, 1))

manhattan(
  subset(SD_POP_XY_nobilis, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralog nobilis POR (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_IT, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralog nobilis IT (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralog nobilis all (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_azorica, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralog azorica (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_novocanariensis, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralog nobilis all (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)



##################################
#### 6. paralogs GENOME WIDE #####
##################################

par(mfrow = c(5, 1))

manhattan(
  subset(SD_POP_XY_nobilis_POR, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralogs nobilis POR",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_IT, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralogs nobilis ITA",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_TUN, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralogs nobilis TUN",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_GRE, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralogs nobilis GRE",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)



manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralogs nobilis all",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)


manhattan(
  subset(SD_POP_XY_azorica, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralog azorica ",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)


manhattan(
  subset(SD_POP_XY_novocanariensis, CHR_NUM %in% 1:12),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_paralog",
  snp = "gene",
  logp = FALSE,
  main = "Manhattan plot – paralogs novocanariensis",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)


par(mfrow = c(1, 1))



###############################
##### 7. Filtre sur CHR6. #####
###############################
SD_POP_XY_azorica_CHR6 <- SD_POP_XY_azorica %>%
  filter(CHR_NUM == 6)

SD_POP_XY_nobilis_all_CHR6 <- SD_POP_XY_nobilis_all %>%
  filter(CHR_NUM == 6)

SD_POP_XY_novocanariensis_CHR6 <- SD_POP_XY_novocanariensis %>%
  filter(CHR_NUM == 6)

SD_POP_XY_azorica_CHR6_sex_genes <- SD_POP_XY_azorica_CHR6 %>%
  filter(noprior_xy > 0.8) 

SD_POP_XY_nobilis_all_CHR6_sex_genes <- SD_POP_XY_nobilis_all_CHR6 %>%
  filter(noprior_xy > 0.8) 

SD_POP_XY_novocanariensis_CHR6_sex_genes <- SD_POP_XY_novocanariensis_CHR6 %>%
  filter(noprior_xy > 0.8) 

SD_POP_XY_commun_CHR6_sex_genes <- SD_POP_XY_azorica_CHR6_sex_genes %>%
  inner_join(SD_POP_XY_nobilis_all_CHR6_sex_genes, by = "gene", suffix = c("_azorica", "_nobilis")) %>%
  inner_join(SD_POP_XY_novocanariensis_CHR6_sex_genes, by = "gene")

par(mfrow = c(1, 1))

#######################################
##### 8. Commun sexed linked gene #####
#######################################

setwd("/Users/cibio/Desktop/Laurus/SDPOP/results_15_04")
azo <- SD_POP_XY_azorica_CHR6 %>% rename_with(~paste0(.x, ".azo"), -1)
novo <- SD_POP_XY_novocanariensis_CHR6 %>% rename_with(~paste0(.x, ".novo"), -1)
nob <- SD_POP_XY_nobilis_all_CHR6 %>% rename_with(~paste0(.x, ".nob"), -1)

sexed_linked_genes_azo_novo_nob_CHR6 <-
  azo %>%
  full_join(novo, by = names(azo)[1]) %>%
  full_join(nob,  by = names(azo)[1]) %>%
  select(gene, scaffold.azo,CHR_NUM.azo,start.azo,end.azo, N_sites.azo ,mean_coverage.azo ,noprior_autosomal.azo, noprior_haploid.azo, noprior_paralog.azo,noprior_xhemizygote.azo,noprior_xy.azo,max.azo,N_sites.novo,mean_coverage.novo,noprior_autosomal.novo,noprior_haploid.novo,noprior_paralog.novo,noprior_xhemizygote.novo,noprior_xy.novo,max.novo,N_sites.nob,mean_coverage.nob,noprior_autosomal.nob,noprior_haploid.nob, noprior_paralog.nob,noprior_xhemizygote.nob,noprior_xy.nob,max.nob) %>%
  filter(
    noprior_xy.azo > 0.8 |
    noprior_xy.novo > 0.8 |
    noprior_xy.nob  > 0.8)

#write.csv(sexed_linked_genes_azo_novo_nob_CHR6, "sexed_linked_genes_azo_novo_nob_CHR6.csv", row.names = TRUE)


azo <- SD_POP_XY_azorica %>% rename_with(~paste0(.x, ".azo"), -1)
novo <- SD_POP_XY_novocanariensis %>% rename_with(~paste0(.x, ".novo"), -1)
nob <- SD_POP_XY_nobilis_all %>% rename_with(~paste0(.x, ".nob"), -1)

sexed_linked_genes_azo_novo_nob <-
  azo %>%
  full_join(novo, by = names(azo)[1]) %>%
  full_join(nob,  by = names(azo)[1]) %>%
  select(gene, scaffold.azo,CHR_NUM.azo,start.azo,end.azo, N_sites.azo ,mean_coverage.azo ,noprior_autosomal.azo, noprior_haploid.azo, noprior_paralog.azo,noprior_xhemizygote.azo,noprior_xy.azo,max.azo,N_sites.novo,mean_coverage.novo,noprior_autosomal.novo,noprior_haploid.novo,noprior_paralog.novo,noprior_xhemizygote.novo,noprior_xy.novo,max.novo,N_sites.nob,mean_coverage.nob,noprior_autosomal.nob,noprior_haploid.nob, noprior_paralog.nob,noprior_xhemizygote.nob,noprior_xy.nob,max.nob) %>%
  filter(
    noprior_xy.azo > 0.8 |
      noprior_xy.novo > 0.8 |
      noprior_xy.nob  > 0.8)

#write.csv(sexed_linked_genes_azo_novo_nob, "sexed_linked_genes_azo_novo_nob.csv", row.names = TRUE)

sexed_linked_genes_azo_novo <-
  azo %>%
  full_join(novo, by = names(azo)[1]) %>%
  full_join(nob,  by = names(azo)[1]) %>%
  select(gene,
         scaffold.azo, CHR_NUM.azo, start.azo, end.azo,
         N_sites.azo, mean_coverage.azo,
         noprior_xy.azo,
         N_sites.novo, mean_coverage.novo,
         noprior_xy.novo,
         N_sites.nob, mean_coverage.nob,
         noprior_xy.nob) %>%
  filter(
      noprior_xy.azo > 0.8 &
      noprior_xy.novo > 0.8 &
      noprior_xy.nob < 0.8
  )

common_genes_azo_novo <- sexed_linked_genes_azo_novo_nob$gene

##################################################################
## 9. SNP diagnostic (pour verifier le sexage des individus) #####
##################################################################

library(vcfR)

vcf_from_az <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP/SNP_diagnostique/diagnostic_snp_from_azorica.vcf")
gt_from_az <- extract.gt(vcf_from_az, element = "GT")

vcf_from_nob_pt <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP/SNP_diagnostique/diagnostic_snp_from_nobilis_port.vcf")
gt_from_nob_pt <- extract.gt(vcf_from_nob_pt, element = "GT")

vcf_novo_from_az <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/L_novocanariensis_L_azorica/SDPOP/SNP_diagnostic/diagnostic_snp_novocanariensis_from_azorica.vcf")
gt_novo_from_az <- extract.gt(vcf_novo_from_az, element = "GT")


count_genotypes <- function(x) {
  c(
    hom_ref = sum(x == "0/0", na.rm = TRUE),
    het     = sum(x %in% c("0/1", "1/0"), na.rm = TRUE),
    hom_alt = sum(x == "1/1", na.rm = TRUE),
    N.A      = sum(is.na(x) | x == "./.")
  )
}

setwd("/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SNP_diagnostic")
counts_from_nob_pt <- apply(gt_from_nob_pt, 2, count_genotypes)
counts_from_nob_pt <- t(counts_from_nob_pt)
#write.csv(counts_from_nob_pt, "counts_from_nob_pt.csv", row.names = TRUE)

counts_from_az <- apply(gt_from_az, 2, count_genotypes)
counts_from_az <- t(counts_from_az)
#write.csv(counts_from_az, "counts_from_az.csv", row.names = TRUE)

counts_novo_from_az <- apply(gt_novo_from_az, 2, count_genotypes)
counts_novo_from_az <- t(counts_novo_from_az)
#write.csv(counts_novo_from_az, "counts_novo_from_az.csv", row.names = TRUE)


###################################
#### 10. Divergence gametologs ####
###################################


WXYZ_nobilis_all <-  read.table("/Users/cibio/Desktop/Laurus/SDPOP/L_nobilis_all_L_azorica/SDPOP_4/wXYz_allgenes_L_nobilis_all_Lazorica_norm_filtered_EXONS.recode.out.contig",comment.char = "", header = TRUE, sep = "", stringsAsFactors = FALSE)
WXYZ_nobilis_all$X.contig_name <- sub("_X", "", WXYZ_nobilis_all$X.contig_name)
WXYZ_nobilis_all$X.contig_name <- sub("^>", "", WXYZ_nobilis_all$X.contig_name) 
WXYZ_nobilis_all <- WXYZ_nobilis_all %>%
  rename(gene = X.contig_name)

WXYZ_novo <-  read.table("/Users/cibio/Desktop/Laurus/SDPOP/L_novocanariensis_L_azorica/SDPOP_resex/wXYz_allgenes_L_novocanariensis_Lazorica_norm_filtered_EXONS.recode.out.contig",comment.char = "", header = TRUE, sep = "", stringsAsFactors = FALSE)
WXYZ_novo$X.contig_name <- sub("_X", "", WXYZ_novo$X.contig_name)
WXYZ_novo$X.contig_name <- sub("^>", "", WXYZ_novo$X.contig_name)
WXYZ_novo <- WXYZ_novo %>%
  rename(gene = X.contig_name)

WXYZ_azo <-  read.table("/Users/cibio/Desktop/Laurus/SDPOP/L_azorica_L_azorica/SD_POP/wXYz_allgenes_L_azorica_Lazorica_norm_filtered_EXONS.recode.out.contig",comment.char = "", header = TRUE, sep = "", stringsAsFactors = FALSE)
WXYZ_azo$X.contig_name <- sub("_X", "", WXYZ_azo$X.contig_name)
WXYZ_azo$X.contig_name <- sub("^>", "", WXYZ_azo$X.contig_name)
WXYZ_azo <- WXYZ_azo %>%
  rename(gene = X.contig_name)


SD_POP_XY_nobilis_all_CHR6_sex_genes <- SD_POP_XY_nobilis_all_CHR6_sex_genes %>%
  left_join(WXYZ_nobilis_all %>% select(gene, divergence),
            by = "gene")

SD_POP_XY_novocanariensis_CHR6_sex_genes <- SD_POP_XY_novocanariensis_CHR6_sex_genes %>%
  left_join(WXYZ_novo %>% select(gene, divergence),
            by = "gene")

SD_POP_XY_azorica_CHR6_sex_genes <- SD_POP_XY_azorica_CHR6_sex_genes %>%
  left_join(WXYZ_azo %>% select(gene, divergence),
            by = "gene")


gene_topology_ML <- c(
  g16665 = "C", g16678 = "D", g16681 = "B", g16684 = "D", g16687 = "D", g16692="F", g16729 = "B", g16731 = "D", g16732 = "B",
  g16736 = "B", g16747 = "B", g16750 = "B", g16774 = "B", g16806 = "B", g16702 = "D", g16717 = "A", g16721 = "A", 
  g16742 = "A", g16746 = "A", g16765 = "D", g16766 = "A", g16767 = "D", g16789 = "A", g16791 = "D", g16793 = "G", 
  g16804 = "A", g16813 = "A", g16814 = "A"
)

gene_topology_Bayesian <- c(
  g16665 = "C", g16678 = "B", g16681 = "B", g16684 = "B", g16687 = "B", g16692="D", g16729 = "B", g16731 = "B", g16732 = "B",
  g16736 = "B", g16747 = "B", g16750 = "B", g16774 = "B", g16806 = "B", g16702 = "A", g16717 = "A", g16721 = "A", 
  g16742 = "D", g16746 = "A", g16765 = "A", g16766 = "D", g16767 = "D", g16789 = "D", g16791 = "D", g16793 = "G", 
  g16804 = "A", g16813 = "D", g16814 = "A"
)

gene_topology <- c(
  g16665 = "C", g16678 = "B", g16681 = "B", g16684 = "B", g16687 = "B", g16692="F", g16729 = "B", g16731 = "B", g16732 = "B",
  g16736 = "B", g16747 = "B", g16750 = "B", g16774 = "B", g16806 = "B", g16702 = "A", g16717 = "A", g16721 = "A", 
  g16742 = "A", g16746 = "A", g16765 = "A", g16766 = "A", g16767 = "D", g16789 = "A", g16791 = "D", g16793 = "G", 
  g16804 = "A", g16813 = "A", g16814 = "A"
)




df_genes <- rbind(
  transform(SD_POP_XY_nobilis_all_CHR6_sex_genes, species = "L.nobilis"),
  transform(SD_POP_XY_novocanariensis_CHR6_sex_genes, species = "L.novocanariensis"),
  transform(SD_POP_XY_azorica_CHR6_sex_genes, species = "L.azorica")
)

df_genes$gene_topology <- gene_topology_Bayesian[df_genes$gene]
df_genes$gene_topology[is.na(df_genes$gene_topology)] <- "E"

df_genes$species <- factor(df_genes$species,
                           levels = c("L.azorica", "L.novocanariensis", "L.nobilis"))

################################
## 10.1 PLOT Ds par espece  ####
################################
ggplot(df_genes, aes(x = start, y = divergence)) +
  geom_point(aes(color = species),
             size = 2,
             alpha = 0.6) +
  # tendance
  geom_smooth(aes(color = species),
              method = "loess",
              se = FALSE,
              linewidth = 1) +
  facet_wrap(~ species, ncol = 1) +
  scale_color_manual(values = c(
    "L.novocanariensis" = "#01558D",
    "L.azorica" = "#AC0C2F",
    "L.nobilis" = "#E6AF2E"
  )) +
  scale_x_continuous(
    breaks = seq(0, 120e6, by = 30e6),
    labels = seq(0, 120, by = 30)
  ) +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "Genomic position (CHR6, Mb)",
    y = "Divergence"
  )

##############################
## 10.2 Plot Ds grouped   ####
##############################

ggplot(df_genes, aes(x = start, y = divergence)) +
  geom_point(aes(color = species),
             size = 1.4,
             alpha = 0.7) +
  # tendance
  geom_smooth(aes(color = species),
              method = "loess",
              se = FALSE,
              linewidth = 1) +
  scale_color_manual(values = c(
    "L.novocanariensis" = "#01558D",
    "L.azorica" = "#AC0C2F",
    "L.nobilis" = "#E6AF2E"
  )) +
  scale_x_continuous(
    breaks = seq(0, 120e6, by = 30e6),
    labels = seq(0, 120, by = 30)
  ) +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "Genomic position (CHR6, Mb)",
    y = "Divergence"
  )


par(mfrow= c(1,1))


################################################
### 10.3 PLOT Ds + TOPOLOGY GAMETOLOGUES.    ###
################################################
par(mfrow= c(1,1))

###ETAPE 1: Dessin des backrounds

#Background zone commune: (azorica | novocanariensis | nobilis)
rect_commun_3 <- data.frame(
  xmin = 45e6, xmax = 94.5e6, ymin = -Inf, ymax = Inf)

#Background zone commune: (azorica | novocanariensis)
rect_commun_2 <- data.frame(
  species = factor(c("L.azorica", "L.novocanariensis"), levels = levels(df_genes$species)),
  xmin = 34.39e6, xmax = 45e6, ymin = -Inf, ymax = Inf)

# Background zones uniques : azorica
rect_uniq_azo <- data.frame( species = factor(c("L.azorica", "L.azorica"), levels = levels(df_genes$species)),
  xmin = c(23.78e6, 94.5e6), xmax = c(34.39e6, 103.78e6), ymin = -Inf, ymax = Inf)

# Background zones uniques : azorica
rect_uniq_azo_2 <- data.frame( species = factor(c("L.azorica"), levels = levels(df_genes$species)),
                             xmin = 1.93e6, xmax = 23.78e6, ymin = -Inf, ymax = Inf)

# Background zones uniques : novocanariensis
rect_uniq_novo <- data.frame( species = factor("L.novocanariensis", levels = levels(df_genes$species)),
                             xmin = 16.8e6, xmax =34.39e6 , ymin = -Inf, ymax = Inf)

# Background zones uniques : nobilis
#rect_uniq_nob <- data.frame( species = factor( "L.nobilis", levels = levels(df_genes$species)),
                              #xmin = 28e6, xmax =44e6 , ymin = -Inf, ymax = Inf)

lines_azorica <- data.frame(
  species = factor("L.azorica", levels = levels(df_genes$species)),
  xintercept = c(23.78e6, 103.78e6)
)


ggplot(df_genes, aes(x = start, y = divergence)) +
  # zone arrière-plan
  geom_rect(
    data = rect_commun_3, inherit.aes = FALSE, 
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
    fill = "darkred", alpha = 0.15
  ) +
  geom_rect( 
    data = rect_commun_2, inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "red", alpha = 0.08
  ) +
  geom_rect( 
    data = rect_uniq_azo, inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "blue", alpha = 0.15
  ) +
  geom_rect( 
    data = rect_uniq_azo_2, inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "blue", alpha = 0.05
  ) +
  geom_rect( 
    data = rect_uniq_novo, inherit.aes = FALSE,
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = "blue", alpha = 0.05
  ) +
  #geom_rect( 
    #data = rect_uniq_nob, inherit.aes = FALSE,
    #aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    #fill = "blue", alpha = 0.05
  #) +
  geom_vline(data = lines_azorica, aes(xintercept = xintercept),
    linetype = "dashed", color = "black", linewidth = 0.8
  ) +
  geom_point(aes(fill = gene_topology, shape = gene_topology), size = 3, alpha = 0.9,stroke = 0.5) +
  geom_smooth(color = "grey40", method = "loess", se = FALSE, linewidth = 0.7) +
  #labs(y = expression(italic(d)[S])) +
  facet_wrap(~ species, ncol = 1) +
  scale_fill_manual(values = c(
    A = "darkred",
    B = "red",
    F = "orange",
    C = "blue",
    D = "grey70",
    E = "grey8",
    G = "grey70"
  )) +
  scale_shape_manual(values = c(
    A = 21,
    B = 21,
    F = 21,
    C = 21,
    D = 21,
    E = 1,
    G = 21
  )) +
  scale_x_continuous(
    breaks = seq(0, 120e6, by = 5e6),
    labels = seq(0, 120, by = 5)
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "italic"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA)
  ) +
  labs(
    x = "Position chromosomique (CHR6, Mb)",
    y = expression(italic(d)[S]),
    fill = "Topology",
    shape = "Topology"
  )

  
####################
## Venn Diagram. ###
####################

library(ggVennDiagram)
library(ggplot2)

# listes de gènes
genes_list <- list(
  azorica = unique(df_genes$gene[df_genes$species == "L.azorica"]),
  novocanariensis = unique(df_genes$gene[df_genes$species == "L.novocanariensis"]),
  nobilis = unique(df_genes$gene[df_genes$species == "L.nobilis"])
)

# plot
library(ggVennDiagram)
library(ggplot2)

genes_list <- list(
  L.azorica = unique(df_genes$gene[df_genes$species == "L.azorica"]),
  L.novocanariensis = unique(df_genes$gene[df_genes$species == "L.novocanariensis"]),
  L.nobilis = unique(df_genes$gene[df_genes$species == "L.nobilis"])
)


ggVennDiagram(genes_list, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "white") +
  theme( legend.position = "none", theme_minimal(base_size = 14)) +
  theme(panel.grid = element_blank(), legend.position = "none")

############################
#### AGE SEX CHR       #####
############################

L_azorica_strat_1_gene= list("g16692","g16716","g16717","g16702","g16721","g16729","g16731","g16732","g16736","g16738","g16742",
                        "g16745","g16746","g16747","g16750", "g16765","g16766","g16767","g16774", "g16788", "g16789","g16791",
                        "g16793","g16795", "g16802","g16804","g16806","g16809","g168013","g16814","g16816","g16817","g16854")
                        
L_novocanariensis_strat_1_gene=list("g16692","g16717","g16702","g16721","g16723","g16729","g16731","g16732","g16736","g16737","g16742","g16746",
                               "g16747","g16750", "g16765","g16766","g16767","g16774","g16789","g16791","g16792","g16793","g16794","g16804",
                               "g16806","g168013","g16814")

L_nobilis_strat_1_gene=list("g16992","g16717","g16721","g16742","g16746","g16766","g16789","g16804","g16813","g16814")
L_azorica_strat_2_gene= list("g16678","g16681","g16684","g16687")
L_novocanariensis_strat_2_gene= list("g16678","g16681","g16684","g16687","g16689")
L_azorica_strat_3_gene= list("g16574","g16590","g16591","g16661","g16665","g16816","g166817","g16854")
#L_azorica_strat_4_gene= list("g16816","g166817","g16854")

#Avec filtre Gametologues > 400nt

#"g16854" 
L_azorica_strat_3_gene_filt= list("g16574","g16590","g16591","g16661","g16665","g16816","g166817")
#g176738
L_azorica_strat_1_gene_filt= list("g16692","g16716","g16717","g16702","g16721","g16729","g16731","g16732","g16736","g16742",
                             "g16745","g16746","g16747","g16750", "g16765","g16766","g16767","g16774", "g16788", "g16789","g16791",
                             "g16793","g16795", "g16802","g16804","g16806","g16809","g168013","g16814","g16816","g16817","g16854")
#g16677
L_novocanariensis_strat_2_gene_filt= list("g16678","g16681","g16684","g16687","g16689")
#g16792
L_novocanariensis_strat_1_gene_filt=list("g16692","g16717","g16702","g16721","g16723","g16729","g16731","g16732","g16736","g16737","g16742","g16746",
                                    "g16747","g16750", "g16765","g16766","g16767","g16774","g16789","g16791","g16793","g16794","g16804",
                                    "g16806","g168013","g16814")


###########

summarise_div <- function(data, genes, species, strate) {
  data %>%
    filter(gene %in% genes) %>%
    summarise(
      mean_divergence   = mean(divergence, na.rm = TRUE),
      median_divergence = median(divergence, na.rm = TRUE),
      top_25_percent    = quantile(divergence, 0.75, na.rm = TRUE),
      top_10_percent    = quantile(divergence, 0.90, na.rm = TRUE),
      top_5_percent     = quantile(divergence, 0.95, na.rm = TRUE)
    ) %>%
    mutate(
      species = species,
      strate  = strate,
      n_genes = length(genes)
    )
}


L_azorica_strat_1 <- summarise_div(WXYZ_azo,L_azorica_strat_1_gene_filt,"azorica",1)

L_azorica_strat_2 <- summarise_div(WXYZ_azo,L_azorica_strat_2_gene,"azorica",2)

L_azorica_strat_3 <- summarise_div(WXYZ_azo,L_azorica_strat_3_gene_filt,"azorica",3)

#L_azorica_strat_4 <- summarise_div(WXYZ_azo,L_azorica_strat_4_gene,"azorica",4)

L_novocanariensis_strat_1 <- summarise_div(WXYZ_novo,L_novocanariensis_strat_1_gene_filt,"novocanariensis",1)

L_novocanariensis_strat_2 <- summarise_div(WXYZ_novo,L_novocanariensis_strat_2_gene_filt,"novocanariensis",2)

L_nobilis_strat_1 <- summarise_div(WXYZ_nobilis_all, L_nobilis_strat_1_gene, "nobilis", 1)


divergence_all <- bind_rows(
  L_azorica_strat_1,
  L_azorica_strat_2,
  L_azorica_strat_3,
  L_novocanariensis_strat_1,
  L_novocanariensis_strat_2,
  L_nobilis_strat_1
)

GT <- 10
u  <- 7 * 10^-9

divergence_all <- divergence_all %>%
  mutate(
    age_years_mean   = (mean_divergence   * GT) / (2 * u),
    age_years_median = (median_divergence * GT) / (2 * u),
    age_years_top25  = (top_25_percent    * GT) / (2 * u),
    age_years_top10  = (top_10_percent    * GT) / (2 * u),
    age_years_top5   = (top_5_percent     * GT) / (2 * u)
  )


#write.csv(divergence_all, "/Users/cibio/Desktop/Laurus/SDPOP/PLOT_FINAL/Panel-Divergence-Topology/divergence_all.csv", row.names = FALSE)

write.csv(divergence_all, "/Users/cibio/Desktop/Laurus/SDPOP/PLOT_FINAL/Panel-Divergence-Topology/divergence_all_filt.csv", row.names = FALSE)


#####################
#### CIRCOS PLOT ####
#####################

library(mgcv)
library(circlize)

### -----------------------
### 1. Ordre chromosomes
### -----------------------

chr_order <- as.character(c(6,7,8,9,10,11,12,1,2,3,4,5))

### -----------------------
### 2. Données
### -----------------------

data_list <- list(
  SD_POP_XY_azorica,
  SD_POP_XY_novocanariensis,
  SD_POP_XY_nobilis_all
)

data_list <- lapply(data_list, function(df){
  df <- df[df$CHR_NUM %in% chr_order, ]
  df$CHR_NUM <- factor(df$CHR_NUM, levels = chr_order)
  df[order(df$CHR_NUM, df$start), ]
})

species_colors <- c(
  "novocanariensis" = "#AC0C2F",
  "azorica" = "#01558D",
  "nobilis" = "#E6AF2E"
)
### -----------------------
### 3. Init cercle
### -----------------------

chr_limits <- do.call(rbind, lapply(data_list, function(df){
  aggregate(start ~ CHR_NUM, df, max)
}))
chr_limits <- aggregate(start ~ CHR_NUM, chr_limits, max)
circos.clear()
circos.par(
  start.degree = 90,
  gap.after = c(rep(2, length(chr_order)-1), 8)
)
circos.initialize(
  factors = chr_limits$CHR_NUM,
  xlim = cbind(0, chr_limits$start)
)

### -----------------------
### 4. Labels chromosomes
### -----------------------

circos.trackPlotRegion(
  ylim = c(0,1),
  track.height = 0.05,
  bg.border = NA,
  panel.fun = function(x, y){
    circos.text(
      CELL_META$xcenter, 0.5,
      paste0("CHR", CELL_META$sector.index),
      facing = "bending.inside",
      niceFacing = TRUE,
      cex = 0.7
    )
  }
)

### -----------------------
### 5. Fonction track
### -----------------------

add_track <- function(df, col_species, add_xaxis = FALSE){
  
  circos.trackPlotRegion(
    factors = df$CHR_NUM,
    ylim = c(0,1),
    track.height = 0.15,
    bg.border = "black",
    
    panel.fun = function(x, y){
      chr = CELL_META$sector.index
      sub_df = df[df$CHR_NUM == chr, ]
      
      # points
      circos.points(sub_df$start, sub_df$noprior_autosomal,
                    col = adjustcolor("grey", alpha.f = 0.3),pch = 16, cex = 0.3)
      
      circos.points(sub_df$start, sub_df$noprior_xy,
                    col = adjustcolor(col_species, alpha.f = 0.8), pch = 16, cex = 0.3)
      
      # seuil
      circos.lines(CELL_META$xlim, c(0.8,0.8), lty = 2)
      
      # Y axis (CHR6)
      if(as.character(chr) == "6"){
        circos.yaxis(
          side = "left",
          at = c(0,0.2,0.4,0.6,0.8,1),
          labels.cex = 0.3,
        )
      }
      
      # X axis
      if(add_xaxis){
        
        xlim = CELL_META$xlim
        ticks <- seq(0, xlim[2], by = 3e7)
        
        circos.axis(
          h = "bottom",
          major.at = ticks,
          labels = round(ticks/1e6),
          labels.cex = 0.3,
          direction = "inside", 
          major.tick.length = 0.04
        )
      }
      
      # regression 
      reg_data<- data.frame(
        x = sub_df$start,
        auto = sub_df$noprior_autosomal,
        xy = sub_df$noprior_xy
      )
      # GAM autosomal
      #fit_auto <- gam(auto ~ s(x, bs = "cs"), data = reg_data)
      # GAM XY
      #fit_xy <- gam(xy ~ s(x, bs = "cs"), data = reg_data)
      # grid de prédiction
      
      # LOESS autosomal
      fit_auto <- loess(auto ~ x, data = reg_data, span = 0.3)
      
      # LOESS XY
      fit_xy <- loess(xy ~ x, data = reg_data, span = 0.3)
      xs <- seq(min(reg_data$x), max(reg_data$x), length.out = 200)
        
      pred_auto <- predict(fit_auto, newdata = data.frame(x = xs))
      pred_xy   <- predict(fit_xy, newdata = data.frame(x = xs))
        
      # 🎯 lignes dans circlize
      circos.lines(xs, pred_auto, col = "darkgrey", lwd = 1.5)
      circos.lines(xs, pred_xy, col = col_species, lwd = 1.5)
    }
  )
}

### -----------------------
### 6. Ajouter tracks
### -----------------------

species_names <- c("azorica","novocanariensis", "nobilis")

for(i in seq_along(data_list)){
  add_track(
    data_list[[i]],
    col_species = species_colors[species_names[i]],
    add_xaxis = (i == 3)
  )
}



##################################### 
####### ZOOM CHR 6 ##################
#####################################

SD_POP_XY_novocanariensis_CHR6$species <- "novocanariensis"
SD_POP_XY_azorica_CHR6$species <- "azorica"
SD_POP_XY_nobilis_all_CHR6$species <- "nobilis"

df_all <- rbind(
  SD_POP_XY_novocanariensis_CHR6,
  SD_POP_XY_azorica_CHR6,
  SD_POP_XY_nobilis_all_CHR6
)
df_all$species <- factor(df_all$species,
                         levels = c("azorica", "novocanariensis","nobilis"))

species_colors <- c(
  "novocanariensis" = "#AC0C2F",
  "azorica" = "#01558D",
  "nobilis" = "#E6AF2E"
)

ggplot(df_all) +
  # autosomal (gris fixe)
  geom_point(aes(x = start, y = noprior_autosomal),
             color = adjustcolor("grey", 0.3),
             size = 1) +
  # XY (coloré par espèce)
  geom_point(aes(x = start, y = noprior_xy, color = species),
             size = 1) +
  geom_smooth(aes(x = start, y = noprior_autosomal),
              method = "loess", span = 0.3,
              se = FALSE, color = "darkgrey", linewidth = 1) +
  geom_smooth(aes(x = start, y = noprior_xy, color = species),
              method = "loess", span = 0.3,
              se = FALSE, linewidth = 1) +
  geom_hline(yintercept = 0.8,
             linetype = "dashed",
             color = "black",
             linewidth = 0.8) +
  facet_wrap(~ species, ncol = 1, labeller = label_value) +
  scale_color_manual(values = species_colors) +
  # X axis: 0–120 Mb every 30
  scale_x_continuous(
    breaks = seq(0, 120e6, by = 30e6),
    labels = seq(0, 120, by = 30)
  ) +
  # Y axis: 0–1
  scale_y_continuous(
    breaks = seq(0, 1, by = 0.2),
    limits = c(0, 1)
  ) +
  labs(
    x = "Position chromosomique (CHR6, Mb)",
    y = "Probabilité postérieure"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_line(color = "black")
  )


#####################
##### COUNT GENE ####
#####################


counts_azo <- data.frame(
  category = c("autosomal", "haploid", "paralog", "xhemizygote", "xy"),
  n_genes = c(
    sum(SD_POP_XY_azorica$noprior_autosomal > 0.8),
    sum(SD_POP_XY_azorica$noprior_haploid > 0.8),
    sum(SD_POP_XY_azorica$noprior_paralog > 0.8),
    sum(SD_POP_XY_azorica$noprior_xhemizygote > 0.8),
    sum(SD_POP_XY_azorica$noprior_xy > 0.8)
  )
)
counts_azo

counts_nob <- data.frame(
  category = c("autosomal", "haploid", "paralog", "xhemizygote", "xy"),
  n_genes = c(
    sum(SD_POP_XY_nobilis_all$noprior_autosomal > 0.8),
    sum(SD_POP_XY_nobilis_all$noprior_haploid > 0.8),
    sum(SD_POP_XY_nobilis_all$noprior_paralog > 0.8),
    sum(SD_POP_XY_nobilis_all$noprior_xhemizygote > 0.8),
    sum(SD_POP_XY_nobilis_all$noprior_xy > 0.8)
  )
)
counts_nob

counts_novo <- data.frame(
  category = c("autosomal", "haploid", "paralog", "xhemizygote", "xy"),
  n_genes = c(
    sum(SD_POP_XY_novocanariensis$noprior_autosomal > 0.8),
    sum(SD_POP_XY_novocanariensis$noprior_haploid > 0.8),
    sum(SD_POP_XY_novocanariensis$noprior_paralog > 0.8),
    sum(SD_POP_XY_novocanariensis$noprior_xhemizygote > 0.8),
    sum(SD_POP_XY_novocanariensis$noprior_xy > 0.8)
  )
)
counts_novo

counts_novo_CHR6 <- data.frame(
  category = c("autosomal", "haploid", "paralog", "xhemizygote", "xy"),
  n_genes = c(
    sum(SD_POP_XY_novocanariensis_CHR6$noprior_autosomal > 0.8),
    sum(SD_POP_XY_novocanariensis_CHR6$noprior_haploid > 0.8),
    sum(SD_POP_XY_novocanariensis_CHR6$noprior_paralog > 0.8),
    sum(SD_POP_XY_novocanariensis_CHR6$noprior_xhemizygote > 0.8),
    sum(SD_POP_XY_novocanariensis_CHR6$noprior_xy > 0.8)
  )
)
counts_novo_CHR6

count_genes <- function(df, species, chr = "all") {
  data.frame(
    species = species,
    chr = chr,
    category = c("autosomal", "haploid", "paralog", "xhemizygote", "xy"),
    n_genes = c(
      sum(df$noprior_autosomal > 0.8, na.rm = TRUE),
      sum(df$noprior_haploid > 0.8, na.rm = TRUE),
      sum(df$noprior_paralog > 0.8, na.rm = TRUE),
      sum(df$noprior_xhemizygote > 0.8, na.rm = TRUE),
      sum(df$noprior_xy > 0.8, na.rm = TRUE)
    )
  )
}

df_counts <- rbind(
  
  # ALL genome
  count_genes(SD_POP_XY_azorica, "azorica", "all"),
  count_genes(SD_POP_XY_nobilis_all, "nobilis", "all"),
  count_genes(SD_POP_XY_novocanariensis, "novocanariensis", "all"),
  
  # CHR6
  count_genes(SD_POP_XY_azorica_CHR6, "azorica", "CHR6"),
  count_genes(SD_POP_XY_nobilis_all_CHR6, "nobilis", "CHR6"),
  count_genes(SD_POP_XY_novocanariensis_CHR6, "novocanariensis", "CHR6")
)


count_genes_max <- function(df, species, chr = "all") {
  
  probs <- df[, c("noprior_autosomal",
                  "noprior_haploid",
                  "noprior_paralog",
                  "noprior_xhemizygote",
                  "noprior_xy")]
  
  # catégorie max par ligne
  max_cat <- apply(probs, 1, function(x) {
    names(x)[which.max(x)]
  })
  
  # nettoyer les noms
  max_cat <- gsub("noprior_", "", max_cat)
  
  # compter
  counts <- as.data.frame(table(max_cat))
  colnames(counts) <- c("category", "n_genes")
  
  # ajouter metadata
  counts$species <- species
  counts$chr <- chr
  
  return(counts)
}


df_counts_max <- rbind(
  
  # genome entier
  count_genes_max(SD_POP_XY_azorica, "azorica", "all"),
  count_genes_max(SD_POP_XY_nobilis_all, "nobilis", "all"),
  count_genes_max(SD_POP_XY_novocanariensis, "novocanariensis", "all"),
  
  # CHR6
  count_genes_max(SD_POP_XY_azorica_CHR6, "azorica", "CHR6"),
  count_genes_max(SD_POP_XY_nobilis_all_CHR6, "nobilis", "CHR6"),
  count_genes_max(SD_POP_XY_novocanariensis_CHR6, "novocanariensis", "CHR6")
)


#####################################
## Mapping all contig VS 12CHR only
####################################



SD_POP_XY_nobilis_all_CHR6_sex_genes <- SD_POP_XY_nobilis_all %>%
  filter(CHR_NUM == 6) %>%
  filter(noprior_xy > 0.8)


SD_POP_XY_nobilis_all_12CHR_CHR6_sexe_genes <- SD_POP_XY_nobilis_all_12CHR %>%
  filter(CHR_NUM == 6) %>%
  filter(noprior_xy > 0.8)



par(mfrow = c(2, 2))
#NOBILIS
manhattan(
  subset(SD_POP_XY_nobilis_all_12CHR, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "XY nobilis all mapping 12CHR (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)

manhattan(
  subset(SD_POP_XY_nobilis_all_12CHR, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "autosomal nobilis all maping 12CHR (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)


#NOBILIS_ALL
manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_xy",
  snp = "gene",
  logp = FALSE,
  main = "XY nobilis all mapping all contigs (CHR 6)",
  ylab = "noprior_xy"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)
manhattan(
  subset(SD_POP_XY_nobilis_all, CHR_NUM == 6),
  chr = "CHR_NUM",
  bp  = "start",
  p   = "noprior_autosomal",
  snp = "gene",
  logp = FALSE,
  main = "Autosomal nobilis all mapping all contigs (CHR 6)",
  ylab = "noprior_autosomal"
)
abline(h = 0.8, col = "red", lty = 2, lwd = 1)







