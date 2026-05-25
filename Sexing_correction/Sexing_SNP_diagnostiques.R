##################################################################
## SNP diagnostic (pour verifier le sexage des individus) #####
##################################################################

library(vcfR)

count_genotypes <- function(x) {
  c(
    hom_ref = sum(x == "0/0", na.rm = TRUE),
    het     = sum(x %in% c("0/1", "1/0"), na.rm = TRUE),
    hom_alt = sum(x == "1/1", na.rm = TRUE),
    N.A      = sum(is.na(x) | x == "./.")
  )
}

#L_nobilis_all_diagnostic_snp_L_nobilis_POR
vcf_L_nobilis_all_diagnostic_snp_L_nobilis_POR <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/Sexing_correction/SNP_diagnostic_L_nobilis_POR/L_nobilis_all_diagnostic_snp_L_nobilis_POR.vcf")
gt_L_nobilis_all_diagnostic_snp_L_nobilis_POR <- extract.gt(vcf_L_nobilis_all_diagnostic_snp_L_nobilis_POR, element = "GT")

#L_novocanariensis_diagnostic_snp_L_nobilis_POR
vcf_L_novocanariensis_diagnostic_snp_L_nobilis_POR <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/Sexing_correction/SNP_diagnostic_L_nobilis_POR/L_novocanariensis_diagnostic_snp_L_nobilis_POR.vcf")
gt_L_novocanariensis_diagnostic_snp_L_nobilis_POR <- extract.gt(vcf_L_novocanariensis_diagnostic_snp_L_nobilis_POR, element = "GT")

#L_nobilis_all_diagnostic_snp_L_azorica
vcf_L_nobilis_all_diagnostic_snp_L_azorica <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/Sexing_correction/SNP_diagnostic_L_azorica/L_nobilis_all_diagnostic_snp_L_azorica.vcf")
gt_L_nobilis_all_diagnostic_snp_L_azorica <- extract.gt(vcf_L_nobilis_all_diagnostic_snp_L_azorica, element = "GT")

#L_novocanariensis_diagnostic_snp_L_azorica
vcf_L_novocanariensis_diagnostic_snp_L_azorica <- read.vcfR("/Users/cibio/Desktop/Laurus/SDPOP/Sexing_correction/SNP_diagnostic_L_azorica/L_novocanariensis_diagnostic_snp_L_azorica.vcf")
gt_L_novocanariensis_diagnostic_snp_L_azorica <- extract.gt(vcf_L_novocanariensis_diagnostic_snp_L_azorica, element = "GT")


## Count from SNP diagnostiques L_nobilis_POR

#Nobilis ALL

counts_L_nobilis_all_diagnostic_snp_L_nobilis_POR <- apply(gt_L_nobilis_all_diagnostic_snp_L_nobilis_POR, 2, count_genotypes)
counts_L_nobilis_all_diagnostic_snp_L_nobilis_POR <- t(counts_L_nobilis_all_diagnostic_snp_L_nobilis_POR)

#Novocanariensis
counts_L_novocanariensis_diagnostic_snp_L_nobilis_POR <- apply(gt_L_novocanariensis_diagnostic_snp_L_nobilis_POR, 2, count_genotypes)
counts_L_novocanariensis_diagnostic_snp_L_nobilis_POR <- t(counts_L_novocanariensis_diagnostic_snp_L_nobilis_POR)

## Count from SNP diagnostiques L_azorica

#Nobilis_all
counts_L_nobilis_all_diagnostic_snp_L_azorica<- apply(gt_L_nobilis_all_diagnostic_snp_L_azorica, 2, count_genotypes)
counts_L_nobilis_all_diagnostic_snp_L_azorica <- t(counts_L_nobilis_all_diagnostic_snp_L_azorica)

#Novocanariensis
counts_L_novocanariensis_diagnostic_snp_L_azorica<- apply(gt_L_novocanariensis_diagnostic_snp_L_azorica, 2, count_genotypes)
counts_L_novocanariensis_diagnostic_snp_L_azorica <- t(counts_L_novocanariensis_diagnostic_snp_L_azorica)



## SAVE TABLES
setwd("/Users/cibio/Desktop/Laurus/SDPOP/Sexing_correction")

write.csv(counts_L_nobilis_all_diagnostic_snp_L_azorica, "counts_L_nobilis_all_diagnostic_snp_L_azorica.csv", row.names = TRUE)
write.csv(counts_L_novocanariensis_diagnostic_snp_L_azorica, "counts_L_novocanariensis_diagnostic_snp_L_azorica.csv", row.names = TRUE)
write.csv(counts_L_nobilis_all_diagnostic_snp_L_nobilis_POR, "counts_L_nobilis_all_diagnostic_snp_L_nobilis_POR.csv", row.names = TRUE)
write.csv(counts_L_novocanariensis_diagnostic_snp_L_nobilis_POR, "counts_L_novocanariensis_diagnostic_snp_L_nobilis_POR.csv", row.names = TRUE)

