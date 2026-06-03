library(vcfR)
library(dplyr)
library(tidyr)
library(ggplot2)




####################
## L. azorica   ####
####################

setwd("/Users/cibio/Desktop/Laurus/SDPOP/Dosage_composation/L_azorica")

#Size of the reads library to normalize the exprssion per individual
library_sizes_L_azorica <- data.frame(
  sample = c(
    "LaF1","LaF2","LaF3","LaF4","LaF5","LaF6","LaF7","LaF8",
    "LaF9","LaF10","LaF11","LaF12","LaF13","LaF14","LaF15","LaF16",
  
    "LaM1","LaM3","LaM4","LaM5","LaM6","LaM7","LaM8","LaM9",
    "LaM10","LaM11","LaM12","LaM13","LaM14","LaM15","LaM16","LaM17","LaM2"
  ),
  size = c(
    # Females
    40843522, 42437979, 57707495, 36563630,
    40038370, 43701962, 49142313, 41721075,
    45007002, 33466785, 49301325, 37445885,
    45477146, 43567566, 54427627, 37396138,
    # Males
    41669962, 43303216, 54480000, 43717315,
    40061684, 44478014, 39693810, 38958918,
    52843891, 45076192, 43520579, 48375786,
    41463911, 48340243, 48904242, 40029867,
    45076192
  )
)

setwd("/Users/cibio/Desktop/Laurus/SDPOP/Dosage_composation/L_azorica")

# =========================================================
# Import du VCF et extraction des AD
# =========================================================

vcf_L_azorica <- read.vcfR("SNP_gametologues.vcf")

allelic_depth_L_azorica <- extract.gt(vcf_L_azorica, element = "AD")

variant_info_L_azorica <- getFIX(vcf_L_azorica)

df_allele_expression_L_azorica <- cbind(variant_info_L_azorica,allelic_depth_L_azorica)

# Suppression des colonnes inutiles
df_allele_expression_L_azorica <- df_allele_expression_L_azorica[,-ncol(df_allele_expression_L_azorica)]
df_allele_expression_L_azorica <- df_allele_expression_L_azorica[,!colnames(df_allele_expression_L_azorica) %in%c("QUAL", "FILTER", "ID")]

# =========================================================
# Dataset males L. azorica
# =========================================================

male_samples <- grep("^LaM",colnames(df_allele_expression_L_azorica),value = TRUE)
df_males_L_azorica <- df_allele_expression_L_azorica[,c("CHROM", "POS", "REF", "ALT", male_samples)]
df_males_L_azorica <- as.data.frame(df_males_L_azorica)

for (sample in male_samples) {
  lib_size <- library_sizes_L_azorica$size[match(sample, library_sizes_L_azorica$sample)]
  split_counts <- df_males_L_azorica %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),sep = ",",convert = TRUE)
  # NORMALISATION
  split_counts <- split_counts / lib_size
  df_males_L_azorica <- cbind(df_males_L_azorica, split_counts)
}

# sans normalisation
for (sample in male_samples) {
  split_counts <- df_males_L_azorica %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),sep = ",",convert = TRUE)
  df_males_L_azorica <- cbind(df_males_L_azorica, split_counts)
}



# Suppression des colonnes originales
df_males_L_azorica <- df_males_L_azorica[,!names(df_males_L_azorica) %in% male_samples]

male_ref_cols <- grep("_ref$",names(df_males_L_azorica),value = TRUE)
male_alt_cols <- grep("_alt$",names(df_males_L_azorica),value = TRUE)

# Résumé des comptages males
male_summary_L_azorica <- df_males_L_azorica %>%
  group_by(CHROM) %>%
  summarise(
    POS = first(as.numeric(POS)),
    X_male = sum(across(all_of(male_ref_cols)), na.rm = TRUE),
    Y_male = sum(across(all_of(male_alt_cols)), na.rm = TRUE)
  ) %>%
  mutate(Y_X_ratio = Y_male / X_male)

male_summary_L_azorica <- df_males_L_azorica %>%
  group_by(CHROM) %>%
  summarise(
    POS = first(as.numeric(POS)),
    X_male = sum(across(all_of(male_ref_cols)), na.rm = TRUE),
    Y_male = sum(across(all_of(male_alt_cols)), na.rm = TRUE)
  ) %>%
  mutate(
    X_tmp = X_male, X_male = if_else(grepl("16591", CHROM), Y_male, X_male),
    Y_male = if_else(grepl("16591", CHROM), X_tmp, Y_male)  #On change le X/Y de gene 16591 car maping Y
  ) %>%
  select(-X_tmp) %>%
  mutate(Y_X_ratio = Y_male / X_male)


# =========================================================
# Dataset females
# =========================================================

female_samples <- grep("^LaF",colnames(df_allele_expression_L_azorica),value = TRUE)
df_females_L_azorica <- df_allele_expression_L_azorica[,c("CHROM", "POS", "REF", "ALT", female_samples)]
df_females_L_azorica <- as.data.frame(df_females_L_azorica)

# Séparation REF / ALT + normalisation
for (sample in female_samples) {
  lib_size <- library_sizes_L_azorica$size[ match(sample, library_sizes_L_azorica$sample)]
  split_counts <- df_females_L_azorica %>%
    select(all_of(sample)) %>%
    separate(sample, into = c(paste0(sample, "_ref"), paste0(sample, "_alt")), sep = ",", convert = TRUE)
  # normalisation par taille de library
  split_counts <- split_counts / lib_size
  df_females_L_azorica <- cbind( df_females_L_azorica, split_counts)
}

# Séparation REF / ALT 
for (sample in female_samples) {
  split_counts <- df_females_L_azorica %>%
    select(all_of(sample)) %>%
    separate(sample, into = c(paste0(sample, "_ref"), paste0(sample, "_alt")), sep = ",", convert = TRUE)
  df_females_L_azorica <- cbind( df_females_L_azorica, split_counts)
}

# Suppression des colonnes originales
df_females_L_azorica <- df_females_L_azorica[,!names(df_females_L_azorica) %in% female_samples]

female_ref_cols <- grep("_ref$",names(df_females_L_azorica),value = TRUE)
female_alt_cols <- grep("_alt$",names(df_females_L_azorica),value = TRUE)

# Résumé des comptages femelles

female_summary_L_azorica <- df_females_L_azorica %>%
  group_by(CHROM) %>%
  summarise(
    female_alt = sum(across(all_of(female_alt_cols)), na.rm = TRUE),
    female_ref = sum(across(all_of(female_ref_cols)), na.rm = TRUE),
    X_female = if_else(female_alt > 10 * female_ref,female_alt, female_ref) #Si il ya bcp de compte alt alors maping sur Y donc X=alt
  )

# =========================================================
# Comparaison males vs femelles
# =========================================================

n_males <- length(male_samples)
n_females <- length(female_samples)

sex_ratio_L_azorica <- male_summary_L_azorica %>%
  select(CHROM, POS, X_male, Y_male, Y_X_ratio) %>%
  left_join(female_summary_L_azorica, by = "CHROM") %>%
  mutate(X_ratio = (X_male / n_males) / (X_female / n_females),
    Sex_biased_ratio =((X_male + Y_male) / n_males) / (X_female / n_females))



############ Option ########

ggplot(sex_ratio_L_azorica, aes(x = Y_X_ratio, y = X_ratio)) +
  geom_point() +
  coord_cartesian(
    xlim = c(0, 1.5),
    ylim = c(0, 1.5)
  ) +
  labs(
    x = "Y/X expression ratio",
    y = "X male / X female"
  ) +
  theme_classic()

#genes_to_remove <- c("g16039", "g15850", "g16591","g16795",	"g16738", "g16352",	"g16351")

sex_ratio_L_azorica <- sex_ratio_L_azorica %>%
  filter(!CHROM %in% genes_to_remove)

ggplot(sex_ratio_L_azorica, aes(x = Y_X_ratio, y = X_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  coord_cartesian(
    xlim = c(0, 1.5),
    ylim = c(0, 2)
  ) +
  labs(
    x = "Y/X expression ratio",
    y = "X male / X female"
  ) +
  theme_classic()


model <- lm(X_ratio ~ Y_X_ratio, data = sex_ratio_L_azorica)
summary(model)

r2 <- summary(model)$r.squared
r2
########################



##########################
### L. novocanariensis ###
##########################


#Size of the reads library to normalize the exprssion per individual
library_sizes_L_novocanariensis <- data.frame(
  sample = c(
    "LN-T1_M","LN-T2_F","LN-T3_F","LN-T4_M","LN-T5_F","LN-T6_M",
    "LN-T7_M","LN-T10_M","LN-T11_M","LN-T12_M","LN-T13_F","LN-T14_M",
    "LN-T15_M","LN-T16_M","LN-T17b_M","LN-T18b_M","LN-T19_M","LN-T20_M",
    "LN-T21_M","LN-T22_M","LN-T23_M","LN-T24_M","LN-T25","LN-T26_F",
    "LN-T27_F","LN-T28_F"
  ),
  size = c(
    59620507,48062093,100388264,70028228,103057344,67572207,
    78038790,65955540,62680258,67013269,56698441,69915724,
    52537066,45167542,76797668,75558706,66779529,62860343,
    68047303,56944283,36758301,65603044,50258522,68324796,
    81417156,61725748
  )
)

setwd("/Users/cibio/Desktop/Laurus/SDPOP/Dosage_composation/L_novocanariensis")

# =========================================================
# Import du VCF et extraction des AD
# =========================================================

vcf_L_novocanariensis <- read.vcfR("SNP_gametologues.vcf")

allelic_depth_L_novocanariensis <- extract.gt(vcf_L_novocanariensis, element = "AD")

variant_info_L_novocanariensis <- getFIX(vcf_L_novocanariensis)

df_allele_expression_L_novocanariensis <- cbind(variant_info_L_novocanariensis,allelic_depth_L_novocanariensis)

# Suppression des colonnes inutiles
df_allele_expression_L_novocanariensis <- df_allele_expression_L_novocanariensis[,-ncol(df_allele_expression_L_novocanariensis)]
df_allele_expression_L_novocanariensis <- df_allele_expression_L_novocanariensis[,!colnames(df_allele_expression_L_novocanariensis) %in%c("QUAL", "FILTER", "ID")]

# =========================================================
# Dataset males L. novocanariensis
# =========================================================

male_samples <- c("LN-T22_M","LN-T20_M","LN-T18b_M","LN-T17b_M","LN-T23_M","LN-T16_M","LN-T24_M","LN-T1_M","LN-T14_M","LN-T11_M","LN-T21_M","LN-T6_M","LN-T7_M","LN-T10_M","LN-T15_M","LN-T12_M","LN-T25")

df_males_L_novocanariensis <- df_allele_expression_L_novocanariensis[,c("CHROM", "POS", "REF", "ALT", male_samples)]
df_males_L_novocanariensis <- as.data.frame(df_males_L_novocanariensis)


for (sample in male_samples) {
  lib_size <- library_sizes_L_novocanariensis$size[match(sample, library_sizes_L_novocanariensis$sample)]
  split_counts <- df_males_L_novocanariensis %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),sep = ",",convert = TRUE)
  # NORMALISATION
  split_counts <- split_counts / lib_size
  df_males_L_novocanariensis <- cbind(df_males_L_novocanariensis, split_counts)
}


### Sans normalisation ####
for (sample in male_samples) {
  split_counts <- df_males_L_novocanariensis %>%
    select(all_of(sample)) %>%
    separate(sample, into = c(paste0(sample, "_ref"), paste0(sample, "_alt")), sep = ",", convert = TRUE)
  df_males_L_novocanariensis <- cbind(df_males_L_novocanariensis, split_counts)
}


# Suppression des colonnes originales
df_males_L_novocanariensis <- df_males_L_novocanariensis[,!names(df_males_L_novocanariensis) %in% male_samples]

male_ref_cols <- grep("_ref$",names(df_males_L_novocanariensis),value = TRUE)
male_alt_cols <- grep("_alt$",names(df_males_L_novocanariensis),value = TRUE)

# Résumé des comptages males
male_summary_L_novocanariensis <- df_males_L_novocanariensis %>%
  group_by(CHROM) %>%
  summarise(
    POS = first(as.numeric(POS)),
    X_male = sum(across(all_of(male_ref_cols)), na.rm = TRUE),
    Y_male = sum(across(all_of(male_alt_cols)), na.rm = TRUE)
  ) %>%
  mutate(
    Y_X_ratio = Y_male / X_male
  )

# =========================================================
# Dataset females L. novocanariensis
# =========================================================

female_samples <- c("LN-T28_F","LN-T27_F","LN-T26_F","LN-T2_F","LN-T3_F","LN-T5_F","LN-T13_F","LN-T4_M","LN-T19_M")

df_females_L_novocanariensis <- df_allele_expression_L_novocanariensis[,c("CHROM", "POS", "REF", "ALT", female_samples)]
df_females_L_novocanariensis <- as.data.frame(df_females_L_novocanariensis)

# Séparation REF / ALT + normalisation
for (sample in female_samples) {
  lib_size <- library_sizes_L_novocanariensis$size[ match(sample, library_sizes_L_novocanariensis$sample)]
  split_counts <- df_females_L_novocanariensis %>%
    select(all_of(sample)) %>%
    separate(sample, into = c(paste0(sample, "_ref"), paste0(sample, "_alt")), sep = ",", convert = TRUE)
  
  # normalisation par taille de library
  split_counts <- split_counts / lib_size
  df_females_L_novocanariensis <- cbind( df_females_L_novocanariensis, split_counts)
}


## Sans normalisation ##
for (sample in female_samples) {
  split_counts <- df_females_L_novocanariensis %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),sep = ",", convert = TRUE)
  df_females_L_novocanariensis <- cbind(df_females_L_novocanariensis, split_counts)
}

# Suppression des colonnes originales
df_females_L_novocanariensis <- df_females_L_novocanariensis[,!names(df_females_L_novocanariensis) %in% female_samples]

female_ref_cols <- grep("_ref$",names(df_females_L_novocanariensis),value = TRUE)
female_alt_cols <- grep("_alt$",names(df_females_L_novocanariensis),value = TRUE)

# Résumé des comptages femelles
female_summary_L_novocanariensis <- df_females_L_novocanariensis %>%
  group_by(CHROM) %>%
  summarise( X_female = sum(across(all_of(female_ref_cols)), na.rm = TRUE),
             female_ref = sum(across(all_of(female_ref_cols)), na.rm = TRUE),
             female_alt = sum(across(all_of(female_alt_cols)), na.rm = TRUE))

# =========================================================
# Comparaison males vs femelles L.novocanariensis
# =========================================================

n_males <- length(male_samples)
n_females <- length(female_samples)

sex_ratio_L_novocanariensis <- male_summary_L_novocanariensis %>%
  select(CHROM, POS, X_male, Y_male , Y_X_ratio) %>%
  left_join(female_summary_L_novocanariensis, by = "CHROM") %>%
  mutate(
    X_ratio = (X_male / n_males) / (X_female / n_females),
    Sex_biased_ratio = ((X_male + Y_male) / n_males) / (X_female / n_females))








##########################
### L. nobilis ###
##########################



#Size of the reads library to normalize the exprssion per individual
library_sizes_L_nobilis <- data.frame(
  sample = c(
    "GRE1-2_5b","GRE3-4_1b","GRE3-4_2b","GRE5_2b","GRE6-7_1",
    
    "Lm-F02Cb","Lm-F04Cbf","Lm-F05Cbf","Lm-F06Cbf","Lm-F07b",
    "Lm-F08Cbf","Lm-F09Cf","Lm-F10Cb","Lm-F11b","Lm-F12Cf",
    "Lm-F13Cbf","Lm-F14Cbf","Lm-F15Cbf","Lm-F16b",
    
    "Lm-M01Cf","Lm-M04b","Lm-M05b","Lm-M06Cb","Lm-M07Cf",
    "Lm-M08Cb","Lm-M09b","Lm-M10b","Lm-M11Cb","Lm-M12Cb",
    "Lm-M13Cf","Lm-M14Cb","Lm-M15Cb",
    
    "OR01F","OR01M","OR02M","OR03Fb","OR03M","OR04F","OR04Mb",
    "OR05F","OR05M","OR06Fb","OR07M","OR08F","OR09F","OR09M",
    "OR10M",
    
    "TUNI06-M","TUNI10-Mb","TUNI11-Fb","TUNI13-Fb","TUNI14-Fb"
  ),
  size = c(
    58577068,56462041,54599420,42813987,69012978,
    
    30385015,29020333,30270985,25734927,31145520,
    34493366,27682951,34068045,19259937,30496436,
    27470390,28488534,23729328,20485129,
    
    27089321,29350846,37929660,30907785,31448928,
    29744411,32074884,34734030,28940112,35449851,
    31291298,31940364,37869950,
    
    69321187,65093878,53968476,45971855,55540404,
    63669247,51049599,54386955,52600002,44902607,
    62013234,63299179,58298812,64120730,
    76435009,
    
    22259628,38271898,57438553,63952965,56605231
  )
)

setwd("/Users/cibio/Desktop/Laurus/SDPOP/Dosage_composation/L_nobilis")

# =========================================================
# Import du VCF et extraction des AD
# =========================================================

vcf_L_nobilis <- read.vcfR("SNP_gametologues_L_nobilis_all.vcf")

allelic_depth_L_nobilis <- extract.gt(vcf_L_nobilis, element = "AD")

variant_info_L_nobilis <- getFIX(vcf_L_nobilis)

df_allele_expression_L_nobilis <- cbind(variant_info_L_nobilis,allelic_depth_L_nobilis)

# Suppression des colonnes inutiles
df_allele_expression_L_nobilis <- df_allele_expression_L_nobilis[,-ncol(df_allele_expression_L_nobilis)]
df_allele_expression_L_nobilis <- df_allele_expression_L_nobilis[,!colnames(df_allele_expression_L_nobilis) %in%c("QUAL", "FILTER", "ID")]


# =========================================================
# Dataset males L. nobilis
# =========================================================

female_samples <- c("Lm-F16b","Lm-F13Cbf","Lm-F12Cf","Lm-F09Cf","Lm-F08Cbf","Lm-F07b","Lm-F05Cbf","OR05F","Lm-F02Cb","Lm-F04Cbf","Lm-F06Cbf","OR01F",
  "Lm-F10Cb","Lm-F03Cb","Lm-F14Cbf","Lm-F15Cbf","Lm-F11b","OR03Fb","OR04F","OR06Fb","OR09F","OR08F","OR03M","OR09M",
  "OR02M","OR04Mb","OR05M","Lm-F01Cb","TUNI06-M","TUNI10-Mb"
)

male_samples <- c("Lm-M14Cb","Lm-M11Cb","Lm-M09b","Lm-M15Cb","Lm-M04b","Lm-M05b","Lm-M13Cf","Lm-M12Cb","Lm-M07Cf","Lm-M10b","Lm-M08Cb","Lm-M06Cb",
  "Lm-M01Cf","OR07M","OR10M","OR01M","TUNI13-Fb","TUNI14-Fb","TUNI11-Fb"
)


df_males_L_nobilis <- df_allele_expression_L_nobilis[,c("CHROM", "POS", "REF", "ALT", male_samples)]
df_males_L_nobilis <- as.data.frame(df_males_L_nobilis)


for (sample in male_samples) {
  lib_size <- library_sizes_L_nobilis$size[match(sample, library_sizes_L_nobilis$sample)]
  split_counts <- df_males_L_nobilis %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),sep = ",",convert = TRUE)
  # NORMALISATION
  split_counts <- split_counts / lib_size
  df_males_L_nobilis <- cbind(df_males_L_nobilis, split_counts)
}


### Sans normalisation ####
for (sample in male_samples) {
  split_counts <- df_males_L_nobilis %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),  sep = ",", convert = TRUE)
  df_males_L_nobilis <- cbind(df_males_L_nobilis, split_counts)
}


# Suppression des colonnes originales
df_males_L_nobilis <- df_males_L_nobilis[,!names(df_males_L_nobilis) %in% male_samples]

male_ref_cols <- grep("_ref$",names(df_males_L_nobilis),value = TRUE)
male_alt_cols <- grep("_alt$",names(df_males_L_nobilis),value = TRUE)

# Résumé des comptages males
male_summary_L_nobilis <- df_males_L_nobilis %>%
  group_by(CHROM) %>%
  summarise(
    POS = first(as.numeric(POS)),
    X_male = sum(across(all_of(male_ref_cols)), na.rm = TRUE),
    Y_male = sum(across(all_of(male_alt_cols)), na.rm = TRUE)
  ) %>%
  mutate(
    Y_X_ratio = Y_male / X_male
  )


# =========================================================
# Dataset females L. nobilis
# =========================================================

female_samples <- c("Lm-F16b","Lm-F13Cbf","Lm-F12Cf","Lm-F09Cf","Lm-F08Cbf","Lm-F07b","Lm-F05Cbf","OR05F","Lm-F02Cb","Lm-F04Cbf","Lm-F06Cbf","OR01F",
                    "Lm-F10Cb","Lm-F03Cb","Lm-F14Cbf","Lm-F15Cbf","Lm-F11b","OR03Fb","OR04F","OR06Fb","OR09F","OR08F","OR03M","OR09M",
                    "OR02M","OR04Mb","OR05M","TUNI06-M","TUNI10-Mb")

df_females_L_nobilis <- df_allele_expression_L_nobilis[, c("CHROM", "POS", "REF", "ALT", female_samples)]
df_females_L_nobilis <- as.data.frame(df_females_L_nobilis)

# Séparation REF / ALT + normalisation
for (sample in female_samples) {
  lib_size <- library_sizes_L_nobilis$size[ match(sample, library_sizes_L_nobilis$sample)]
  split_counts <- df_females_L_nobilis %>%
    select(all_of(sample)) %>%
    separate(sample, into = c(paste0(sample, "_ref"), paste0(sample, "_alt")), sep = ",", convert = TRUE)
  
  # normalisation par taille de library
  split_counts <- split_counts / lib_size
  df_females_L_nobilis <- cbind( df_females_L_nobilis, split_counts)
}


## Sans normalisation ##
for (sample in female_samples) {
  split_counts <- df_females_L_nobilis %>%
    select(all_of(sample)) %>%
    separate(sample,into = c(paste0(sample, "_ref"), paste0(sample, "_alt")),sep = ",",convert = TRUE)
  df_females_L_nobilis <- cbind(df_females_L_nobilis, split_counts)
}


# Suppression des colonnes originales
df_females_L_nobilis <- df_females_L_nobilis[,!names(df_females_L_nobilis) %in% female_samples]

female_ref_cols <- grep("_ref$",names(df_females_L_nobilis),value = TRUE)
female_alt_cols <- grep("_alt$",names(df_females_L_nobilis),value = TRUE)

# Résumé des comptages femelles
female_summary_L_nobilis <- df_females_L_nobilis %>%
  group_by(CHROM) %>%
  summarise( X_female = sum(across(all_of(female_ref_cols)), na.rm = TRUE),
             female_ref = sum(across(all_of(female_ref_cols)), na.rm = TRUE),
             female_alt = sum(across(all_of(female_alt_cols)), na.rm = TRUE))


# =========================================================
# Comparaison males vs femelles
# =========================================================

n_males <- length(male_samples)
n_females <- length(female_samples)

sex_ratio_L_nobilis <- male_summary_L_nobilis %>%
  select(CHROM, POS, X_male, Y_male, Y_X_ratio) %>%
  left_join(female_summary_L_nobilis, by = "CHROM") %>%
  mutate(
    X_ratio = (X_male / n_males) / (X_female / n_females),
    Sex_biased_ratio =((X_male + Y_male) / n_males) /(X_female / n_females)
    )



###### PLOT FINAUX #######

####################

bg_alpha <- 0.4

bg_azorica <- data.frame(
  xmin = c(1.93e6, 23.78e6, 34.39e6, 45e6, 94.5e6),
  xmax = c(23.78e6, 34.39e6, 45e6, 94.5e6, 103.78e6),
  fill = c("#D6E6FA", "#8FB6E8", "#F8D9D6", "#F28B82", "#8FB6E8")
)

bg_novo <- data.frame(
  xmin = c(16.8e6, 34.39e6, 45e6),
  xmax = c(34.39e6, 45e6, 94.5e6),
  fill = c("#D6E6FA", "#F8D9D6", "#F28B82")
)

bg_nobilis <- data.frame(
  xmin = 45e6,
  xmax = 94.5e6,
  fill = "#F28B82"
)

Gene_shared <- c(
  "g16665", "g16678", "g16681", "g16684", "g16687",
  "g16692", "g16729", "g16731", "g16732", "g16736",
  "g16747", "g16750", "g16774", "g16806", "g16702",
  "g16717", "g16721", "g16742", "g16746", "g16765",
  "g16766", "g16767", "g16789", "g16791", "g16793",
  "g16804", "g16813", "g16814"
)

sex_ratio_L_azorica$shared_status <- ifelse(
  sex_ratio_L_azorica$CHROM %in% Gene_shared,
  "shared",
  "unique"
)

sex_ratio_L_novocanariensis$shared_status <- ifelse(
  sex_ratio_L_novocanariensis$CHROM %in% Gene_shared,
  "shared",
  "unique"
)

sex_ratio_L_nobilis$shared_status <- ifelse(
  sex_ratio_L_nobilis$CHROM %in% Gene_shared,
  "shared",
  "unique"
)





plot_ratio <- function(data, yvar, y_lab = NULL, title, bg_df) {
  
  p <- ggplot(data, aes(x = POS, y = .data[[yvar]])) +
    # BACKGROUND UNIQUE
    geom_rect(
      data = bg_df,
      aes(xmin = xmin, xmax = xmax,
          ymin = -Inf, ymax = Inf,
          fill = fill),
      inherit.aes = FALSE,
      alpha = bg_alpha
    ) +
    scale_fill_identity() +
    # GRILLE
    geom_hline(yintercept = seq(0, 2, by = 0.5),
               color = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = seq(0, 110e6, by = 10e6),
               color = "grey80", linewidth = 0.3) +
    # POINTS
    geom_point() 
    # AXES
    coord_cartesian(ylim = c(0, 2)) +
    scale_x_continuous(
      limits = c(0, 110e6),
      breaks = seq(0, 110e6, by = 10e6),
      labels = function(x) x / 1e6
    ) +
    labs(
      title = title,
      x = "Position (Mb)",
      y = y_lab
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA)
    )
  # =====================================================
  
  # LIGNES DES ATTENDUES DIPLOIDES/TETRAPLOIDES SELON LE RATIO
  
  # =====================================================

  if (yvar == "Y_X_ratio") {
    p <- p +
      geom_hline(yintercept = 1, color = "black", linewidth = 0.5 , linetype = "dashed") +
      geom_hline(yintercept = 0.33, color = "black", linewidth = 0.5)
  }
  
  if (yvar == "X_ratio") {
    p <- p +
      geom_hline(yintercept = 0.5, color = "black", linewidth = 0.5, linetype = "dashed") +
      geom_hline( yintercept = 0.75, color = "black", linewidth = 0.5)
  }
  return(p)
}



library(patchwork)

p1 <- plot_ratio(male_summary_L_azorica, "Y_X_ratio", expression(italic("L. azorica")), y_lab = "Ratio expression Y / X", bg_azorica)
p2 <- plot_ratio(male_summary_L_novocanariensis, "Y_X_ratio", expression(italic("L. novocanariensis")), y_lab = "Ratio expression Y / X", bg_novo)
p3 <- plot_ratio(male_summary_L_nobilis, "Y_X_ratio", expression(italic("L. nobilis")), y_lab = "Ratio expression Y / X",  bg_nobilis)

p4 <- plot_ratio(sex_ratio_L_azorica, "X_ratio", expression(italic("L. azorica")), y_lab = "Ratio expression Xm / Xf" ,bg_azorica)
p5 <- plot_ratio(sex_ratio_L_novocanariensis,"X_ratio", expression(italic("L. novocanariensis")), y_lab = "Ratio expression Xm / Xf", bg_novo) 
p6 <- plot_ratio(sex_ratio_L_nobilis,"X_ratio",expression(italic("L. nobilis")) , y_lab = "Ratio expression Xm / Xf", bg_nobilis) 

p7 <- plot_ratio(sex_ratio_L_azorica, "Sex_biased_ratio", expression(italic("L. azorica")), y_lab = "Ratio expression (Xm+Ym) / Xf" ,bg_azorica)
p8 <- plot_ratio(sex_ratio_L_novocanariensis,"Sex_biased_ratio", expression(italic("L. novocanariensis")), y_lab = "Ratio expression (Xm+Ym) / Xf", bg_novo) 
p9 <- plot_ratio(sex_ratio_L_nobilis,"Sex_biased_ratio",expression(italic("L. nobilis")) , y_lab = "Ratio expression (Xm+Ym) / Xf", bg_nobilis) 


final_plot <-
  (p1 | p4) /
  (p2 | p5) /
  (p3 | p6)


print(final_plot)

final_plot_plus <-
  (p1 | p4 | p7) /
  (p2 | p5 | p8) /
  (p3 | p6 | p9)

print(final_plot_plus)





#### STATS ######

library(dplyr)

assign_strata <- function(data, bg_df) {
  data$strata <- NA
  for(i in 1:nrow(bg_df)) {
    data$strata[
      data$POS >= bg_df$xmin[i] &
        data$POS <=  bg_df$xmax[i]
    ] <- paste0("Strate_", i)
  }
  return(data)
}


#Azorica
sex_ratio_L_azorica <- assign_strata(sex_ratio_L_azorica, bg_azorica)
stats_azorica <- sex_ratio_L_azorica %>%
  group_by(strata) %>%
  summarise(
    mean_X_ratio = mean(X_ratio, na.rm = TRUE),
    median_X_ratio = median(X_ratio, na.rm = TRUE),
    mean_Y_ratio = mean(Y_X_ratio, na.rm = TRUE),
    median_Y_ratio = median(Y_X_ratio, na.rm = TRUE)
    
  )

stats_azorica

#Novocanariensis
sex_ratio_L_novocanariensis <- assign_strata(sex_ratio_L_novocanariensis, bg_novo)

stats_novo <- sex_ratio_L_novocanariensis %>%
  group_by(strata) %>%
  summarise(
    mean_X_ratio = mean(X_ratio, na.rm = TRUE),
    median_X_ratio = median(X_ratio, na.rm = TRUE),
    mean_Y_ratio = mean(Y_X_ratio, na.rm = TRUE),
    median_Y_ratio = median(Y_X_ratio, na.rm = TRUE)
    
  )

stats_novo


#Nobilis
sex_ratio_L_nobilis <- assign_strata(sex_ratio_L_nobilis, bg_nobilis)
stats_nobilis <- sex_ratio_L_nobilis %>%
  group_by(strata) %>%
  summarise(
    mean_X_ratio = mean(X_ratio, na.rm = TRUE),
    median_X_ratio = median(X_ratio, na.rm = TRUE),
    mean_Y_ratio = mean(Y_X_ratio, na.rm = TRUE),
    median_Y_ratio = median(Y_X_ratio, na.rm = TRUE)
  )

stats_nobilis




### SANS STRATS

#Azorica
sex_ratio_L_azorica <- assign_strata(sex_ratio_L_azorica, bg_azorica)
stats_azorica <- sex_ratio_L_azorica %>%
  filter(!CHROM %in% c("g16580", "g16039")) %>%
  summarise(
    mean_X_ratio = mean(X_ratio, na.rm = TRUE),
    median_X_ratio = median(X_ratio, na.rm = TRUE),
    mean_Y_ratio = mean(Y_X_ratio, na.rm = TRUE),
    median_Y_ratio = median(Y_X_ratio, na.rm = TRUE)
    
  )

stats_azorica

#Novocanariensis
sex_ratio_L_novocanariensis <- assign_strata(sex_ratio_L_novocanariensis, bg_novo)

stats_novo <- sex_ratio_L_novocanariensis %>%
  summarise(
    mean_X_ratio = mean(X_ratio, na.rm = TRUE),
    median_X_ratio = median(X_ratio, na.rm = TRUE),
    mean_Y_ratio = mean(Y_X_ratio, na.rm = TRUE),
    median_Y_ratio = median(Y_X_ratio, na.rm = TRUE)
    
  )

stats_novo


#Nobilis
sex_ratio_L_nobilis <- assign_strata(sex_ratio_L_nobilis, bg_nobilis)
stats_nobilis <- sex_ratio_L_nobilis %>%
  summarise(
    mean_X_ratio = mean(X_ratio, na.rm = TRUE),
    median_X_ratio = median(X_ratio, na.rm = TRUE),
    mean_Y_ratio = mean(Y_X_ratio, na.rm = TRUE),
    median_Y_ratio = median(Y_X_ratio, na.rm = TRUE)
  )

stats_nobilis


### COMMUN DC ##

DC_candidate_azo <- sex_ratio_L_azorica %>% 
  filter(X_ratio > 0.85) %>%
  mutate(species = "azorica")

DC_candidate_nov <- sex_ratio_L_novocanariensis %>% 
  filter(X_ratio > 0.85) %>%
  mutate(species = "novocanariensis")

DC_candidate_nob <- sex_ratio_L_nobilis %>% 
  filter(X_ratio > 0.85) %>%
  mutate(species = "nobilis")

# merge all together

DC_candidates_all <- bind_rows(DC_candidate_azo, DC_candidate_nov, DC_candidate_nob)
