
## HET ##

L_azorica_L_azorica_het <- read.table("/Users/cibio/Desktop/Laurus/nQuire/L_azorica_Lazorica_VCFALLSITES_C2_REPMONO_NBESTALLELES2.het.het", header = TRUE)
L_azorica_L_nobilis_het <- read.table("/Users/cibio/Desktop/Laurus/nQuire/L_azorica_Lnobilis_VCFALLSITES_C2_REPMONO_NBESTALLELES2_norm.recode.het.het", header = TRUE)

L_nobilis_L_nobilis_het <- read.table("/Users/cibio/Desktop/Laurus/nQuire/L_nobilis_Lnobilis_VCFALLSITES_C2_REPMONO_NBESTALLELES2_norm.recode.het.het", header = TRUE)
L_nobilis_L_azorica_het <- read.table("/Users/cibio/Desktop/Laurus/nQuire/L_nobilis_Lazorica_VCFALLSITES_C2_REPMONO_NBESTALLELES2.recode.het.het", header = TRUE)

L_noblis_het <- merge(L_nobilis_L_azorica_het,L_nobilis_L_nobilis_het,by.x = 1,by.y = 1,suffixes = c("_ref_az", "_ref_nob"))
L_azorica_het <- merge(L_azorica_L_azorica_het,L_azorica_L_nobilis_het,by.x=1,by.y=1,suffixes = c("_ref_az", "_ref_nob"))

## nQuire ##

nQuire_L_azorica <- read.table("/Users/cibio/Desktop/Laurus/nQuire/L_azorica_L_azorica/all_lrd_denoised.txt", header = TRUE)
nQuire_L_nobilis <- read.table("/Users/cibio/Desktop/Laurus/nQuire/L_nobilis_L_nobilis/all_lrd_denoised.txt", header = TRUE)


add_model_columns <- function(df) {

  # Find best model
  df$BestModel <- apply(df[, c("d_dip", "d_tri", "d_tet")], 1, function(x) {
    c("dip", "tri", "tet")[which.min(x)]
  })

  df$DeltaBest2 <- apply(df[, c("d_dip", "d_tri", "d_tet")], 1, function(x) {
    sorted <- sort(x)
    sorted[2] - sorted[1]
  }) / df$free
  

  return(df)
}

# Appliquer aux deux tableaux
nQuire_L_azorica <- add_model_columns(nQuire_L_azorica)
nQuire_L_nobilis <- add_model_columns(nQuire_L_nobilis)
``
