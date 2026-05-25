library(dplyr)
library(tidyr)
library(ggplot2)


#############
### Pixy ####
#############

Pixy_Dxy <- read.table("/Users/cibio/Desktop/Laurus/Pop_gen/Diversity/Pixy/pixy_dxy.txt", header = TRUE)
Pixy_Fst <- read.table("/Users/cibio/Desktop/Laurus/Pop_gen/Diversity/Pixy/pixy_fst.txt", header = TRUE)
Pixy_Pi <- read.table("/Users/cibio/Desktop/Laurus/Pop_gen/Diversity/Pixy/pixy_pi.txt", header = TRUE)

pop <- c(
  "L_novocanariensis",
  "L_azorica",
  "L_nobilis_POR",
  "L_nobilis_ITA",
  "L_nobilis_GRE",
  "L_nobilis_TUN"
)



Pixy_Pi_genome <- Pixy_Pi %>%
  group_by(pop) %>%
  summarise(
    genome_pi = sum(count_diffs, na.rm = TRUE) /
      sum(count_comparisons, na.rm = TRUE),
    total_diffs = sum(count_diffs, na.rm = TRUE),
    total_comparisons = sum(count_comparisons, na.rm = TRUE),
    n_windows = n()
  )


#### FST ####

fst_genome <- Pixy_Fst %>%
  group_by(pop1, pop2) %>%
  summarise(
    fst_genome = sum(avg_wc_fst * no_snps, na.rm = TRUE) /
      sum(no_snps, na.rm = TRUE),
    total_snps = sum(no_snps, na.rm = TRUE),
    n_windows = n()
  )


fst_genome <- fst_genome %>%
  mutate(
    popA = pmin(pop1, pop2),
    popB = pmax(pop1, pop2)
  ) %>%
  group_by(popA, popB) %>%
  summarise(fst = mean(fst_genome), .groups = "drop")

all_pops <- sort(unique(c(fst_genome$popA, fst_genome$popB)))

fst_matrix <- expand.grid(pop1 = all_pops, pop2 = all_pops)

fst_matrix <- merge(
  fst_matrix,
  fst_genome,
  by.x = c("pop1", "pop2"),
  by.y = c("popA", "popB"),
  all.x = TRUE
)

fst_matrix$fst[is.na(fst_matrix$fst)] <- 0


ggplot(fst_matrix, aes(pop1, pop2, fill = fst)) +
  geom_tile() +
  geom_text(aes(label = round(fst, 3)), size = 3) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(fill = "FST")

#### DXY ####


Pixy_Dxy_genome <- Pixy_Dxy %>%
  group_by(pop1, pop2) %>%
  summarise(
    genome_Dxy = sum(count_diffs, na.rm = TRUE) /
      sum(count_comparisons, na.rm = TRUE),
    total_diffs = sum(count_diffs, na.rm = TRUE),
    total_comparisons = sum(count_comparisons, na.rm = TRUE),
    n_windows = n()
  )


Dxy_genome <- Pixy_Dxy_genome %>%
  mutate(
    popA = pmin(pop1, pop2),
    popB = pmax(pop1, pop2)
  ) %>%
  group_by(popA, popB) %>%
  summarise(Dxy = mean(genome_Dxy), .groups = "drop")

all_pops <- sort(unique(c(Dxy_genome$popA, Dxy_genome$popB)))

Dxy_matrix <- expand.grid(pop1 = all_pops, pop2 = all_pops)

Dxy_matrix <- merge(
  Dxy_matrix,
  Dxy_genome,
  by.x = c("pop1", "pop2"),
  by.y = c("popA", "popB"),
  all.x = TRUE
)

Dxy_matrix$fst[is.na(Dxy_matrix$Dxy)] <- 0

ggplot(Dxy_matrix, aes(pop1, pop2, fill = Dxy)) +
  geom_tile() +
  geom_text(aes(label = round(Dxy, 5)), size = 3) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  labs(fill = "Dxy")




########################
###   VCFTOOLS PI   ####
########################

setwd("/Users/cibio/Desktop/Laurus/Pop_gen/Diversity/pi_vcftools")

#AZORICA
denominator <- 13626787

Pi_azorica_SYN <- read.table("synonymous.L_azorica.sites.pi", header = TRUE)
Pi_azorica_NON_SYN <- read.table("nonsynonymous.L_azorica.sites.pi", header = TRUE)
Pi_azorica_global <- read.table("global.L_azorica.sites.pi", header = TRUE)

Pi_azorica_SYN_tot <- sum(Pi_azorica_SYN$PI, na.rm = TRUE) / denominator*0.7
Pi_azorica_NONSYN_tot <- sum(Pi_azorica_NON_SYN$PI, na.rm = TRUE) / denominator*0.3
Pi_azorica_global_tot <- sum(Pi_azorica_global$PI, na.rm = TRUE) / denominator

#NOBILIS GRE
Pi_nobilis_GRE_SYN <- read.table("synonymous.L_nobilis_GRE.sites.pi", header = TRUE)
Pi_nobilis_GRE_NON_SYN <- read.table("nonsynonymous.L_nobilis_GRE.sites.pi", header = TRUE)
Pi_nobilis_GRE_global <- read.table("global.L_nobilis_GRE.sites.pi", header = TRUE)

Pi_nobilis_GRE_SYN_tot <- sum(Pi_nobilis_GRE_SYN$PI, na.rm = TRUE) / denominator*0.7
Pi_nobilis_GRE_NONSYN_tot <- sum(Pi_nobilis_GRE_NON_SYN$PI, na.rm = TRUE) / denominator*0.3
Pi_nobilis_GRE_global_tot <- sum(Pi_nobilis_GRE_global$PI, na.rm = TRUE) / denominator

#NOBILIS ITA
Pi_nobilis_ITA_SYN <- read.table("synonymous.L_nobilis_ITA.sites.pi", header = TRUE)
Pi_nobilis_ITA_NON_SYN <- read.table("nonsynonymous.L_nobilis_ITA.sites.pi", header = TRUE)
Pi_nobilis_ITA_global <- read.table("global.L_nobilis_ITA.sites.pi", header = TRUE)

Pi_nobilis_ITA_SYN_tot <- sum(Pi_nobilis_ITA_SYN$PI, na.rm = TRUE) / denominator*0.7
Pi_nobilis_ITA_NONSYN_tot <- sum(Pi_nobilis_ITA_NON_SYN$PI, na.rm = TRUE) / denominator*0.3
Pi_nobilis_ITA_global_tot <- sum(Pi_nobilis_ITA_global$PI, na.rm = TRUE) / denominator

#NOBILIS POR
Pi_nobilis_POR_SYN <- read.table("synonymous.L_nobilis_POR.sites.pi", header = TRUE)
Pi_nobilis_POR_NON_SYN <- read.table("nonsynonymous.L_nobilis_POR.sites.pi", header = TRUE)
Pi_nobilis_POR_global <- read.table("global.L_nobilis_POR.sites.pi", header = TRUE)

Pi_nobilis_POR_SYN_tot <- sum(Pi_nobilis_POR_SYN$PI, na.rm = TRUE) / denominator*0.7
Pi_nobilis_POR_NONSYN_tot <- sum(Pi_nobilis_POR_NON_SYN$PI, na.rm = TRUE) / denominator*0.3
Pi_nobilis_POR_global_tot <- sum(Pi_nobilis_POR_global$PI, na.rm = TRUE) / denominator

#NOBILIS TUN
Pi_nobilis_TUN_SYN <- read.table("synonymous.L_nobilis_TUN.sites.pi", header = TRUE)
Pi_nobilis_TUN_NON_SYN <- read.table("nonsynonymous.L_nobilis_TUN.sites.pi", header = TRUE)
Pi_nobilis_TUN_global <- read.table("global.L_nobilis_TUN.sites.pi", header = TRUE)

Pi_nobilis_TUN_SYN_tot <- sum(Pi_nobilis_TUN_SYN$PI, na.rm = TRUE) / denominator*0.7
Pi_nobilis_TUN_NONSYN_tot <- sum(Pi_nobilis_TUN_NON_SYN$PI, na.rm = TRUE) / denominator*0.3
Pi_nobilis_TUN_global_tot <- sum(Pi_nobilis_TUN_global$PI, na.rm = TRUE) / denominator

#NOBLIS CANARIENSIS
Pi_novocanariensis_SYN <- read.table("synonymous.L_novocanariensis.sites.pi", header = TRUE)
Pi_novocanariensis_NON_SYN <- read.table("nonsynonymous.L_novocanariensis.sites.pi", header = TRUE)
Pi_novocanariensis_global <- read.table("global.L_novocanariensis.sites.pi", header = TRUE)

Pi_novocanariensis_SYN_tot <- sum(Pi_novocanariensis_SYN$PI, na.rm = TRUE) / denominator*0.7
Pi_novocanariensis_NONSYN_tot <- sum(Pi_novocanariensis_NON_SYN$PI, na.rm = TRUE) / denominator*0.3
Pi_novocanariensis_global_tot <- sum(Pi_novocanariensis_global$PI, na.rm = TRUE) / denominator


Pi_summary <- data.frame(
  species = c(
    "L_azorica",
    "L_nobilis_GRE",
    "L_nobilis_ITA",
    "L_nobilis_POR",
    "L_nobilis_TUN",
    "L_novocanariensis"
  ),
  
  pi_synonymous = c(
    Pi_azorica_SYN_tot,
    Pi_nobilis_GRE_SYN_tot,
    Pi_nobilis_ITA_SYN_tot,
    Pi_nobilis_POR_SYN_tot,
    Pi_nobilis_TUN_SYN_tot,
    Pi_novocanariensis_SYN_tot
  ),
  
  pi_nonsynonymous = c(
    Pi_azorica_NONSYN_tot,
    Pi_nobilis_GRE_NONSYN_tot,
    Pi_nobilis_ITA_NONSYN_tot,
    Pi_nobilis_POR_NONSYN_tot,
    Pi_nobilis_TUN_NONSYN_tot,
    Pi_novocanariensis_NONSYN_tot
  ),
  
  pi_global = c(
    Pi_azorica_global_tot,
    Pi_nobilis_GRE_global_tot,
    Pi_nobilis_ITA_global_tot,
    Pi_nobilis_POR_global_tot,
    Pi_nobilis_TUN_global_tot,
    Pi_novocanariensis_global_tot
  )
)

# afficher le tableau
print(Pi_summary)


# transformer en matrice (sans la colonne species)
mat <- as.matrix(Pi_summary[, -1])

# barplot groupé
barplot(t(mat),
        beside = TRUE,
        names.arg = Pi_summary$species,
        las = 2,
        ylab = "Pi",
        xlab = "Population",
        legend.text = c("Synonymous", "Nonsynonymous", "Global"))

barplot(t(mat),
        beside = TRUE,
        names.arg = Pi_summary$species,
        las = 2,
        ylab = "Pi",
        col = c("lightblue", "salmon", "lightgreen"),
        legend.text = c("Synonymous", "Nonsynonymous", "Global"),
        args.legend = list(x = "topright", bty = "n"))



#############################
#### PiNPiS analyses. #######
#############################


####################
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)


setwd("/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/STATS_R")

###################
## 1. FUNCTION   ##
###################

read_clean_pi_data <- function(file,pi_prefix) {
    read.table(file, header = TRUE) %>%
    # nettoyage header suite au cat
    filter(Contig_name != "Contig_name") %>%
    # split chromosome/gene/CDS
    separate(Contig_name,
             into = c("CDS", "GENE", "CHR"),
             sep = "_",
             remove = FALSE) %>%
    # conversion des colonnes utiles en numeric
    #mutate(across(
      #c(paste0(pi_prefix, c("_piN", "_piS", "_nb_complete_site")),
        #"BiAllelic_SNP"),
     # ~ as.numeric(.)
      mutate(across(
        -c(Contig_name, CDS, GENE, CHR, Status),
        ~ as.numeric(.)
      ))
}

compute_pi_stats <- function(df, chr_remove = "OZ345845.1", piN_col, piS_col, site_col) {
  df %>%
    filter(CHR != chr_remove) %>%
    summarise(
      nb_genes = n_distinct(GENE),
      nb_CDS = n_distinct(CDS),
      nb_sites = sum(.data[[site_col]], na.rm = TRUE),
      nb_SNP = sum(.data[["BiAllelic_SNP"]] , na.rm = TRUE),
      piN = sum(.data[[piN_col]] * .data[[site_col]], na.rm = TRUE) / sum(.data[[site_col]], na.rm = TRUE),
      piS = sum(.data[[piS_col]] * .data[[site_col]], na.rm = TRUE) / sum(.data[[site_col]], na.rm = TRUE)
    ) %>%
    mutate(piN_piS = piN / piS)
}

bootstrap_pi <- function(df, piN_col, piS_col, site_col, B = 1000, seed = 1) {
  set.seed(seed)
  n <- nrow(df)
  res <- replicate(B, {
    idx <- sample(1:n, size = n, replace = TRUE)
    d <- df[idx, ]
    
    piN <- sum(d[[piN_col]] * d[[site_col]], na.rm = TRUE) /
      sum(d[[site_col]], na.rm = TRUE)
    piS <- sum(d[[piS_col]] * d[[site_col]], na.rm = TRUE) /
      sum(d[[site_col]], na.rm = TRUE)
    
    ratio <- piN / piS
    
    c(piN = piN, piS = piS, ratio = ratio)
  })
  
  res <- t(res)
  
  data.frame(
    piN_mean = mean(res[, "piN"]),
    piS_mean = mean(res[, "piS"]),
    ratio_mean = mean(res[, "ratio"]),
    
    piN_CI_low = quantile(res[, "piN"], 0.025),
    piN_CI_high = quantile(res[, "piN"], 0.975),
    
    piS_CI_low = quantile(res[, "piS"], 0.025),
    piS_CI_high = quantile(res[, "piS"], 0.975),
    
    ratio_CI_low = quantile(res[, "ratio"], 0.025),
    ratio_CI_high = quantile(res[, "ratio"], 0.975)
  )
}

compute_FIT_stats <- function(df, chr_remove = "OZ345845.1", fit_col, wc_fit_col, site_col) {
  df %>%
    mutate(across(all_of(c(fit_col, wc_fit_col, site_col)), as.numeric)) %>%
    # On suprimmes les lignes -999 
    filter(!is.na(.data[[fit_col]]) & !is.na(.data[[wc_fit_col]]) &
        .data[[fit_col]] != -999 & .data[[wc_fit_col]] != -999) %>%
    filter(CHR != chr_remove) %>%
    summarise(
      nb_genes = n_distinct(GENE),
      nb_CDS = n_distinct(CDS),
      nb_sites = sum(.data[[site_col]], na.rm = TRUE),
      # moyenne pondérée FIT (simple)
      FIT = sum(.data[[fit_col]] * .data[[site_col]], na.rm = TRUE) /
        sum(.data[[site_col]], na.rm = TRUE),
      # moyenne pondérée Weir & Cockerham FIT
      FIT_WC = sum(.data[[wc_fit_col]] * .data[[site_col]], na.rm = TRUE) /
        sum(.data[[site_col]], na.rm = TRUE)
    )
}

compute_genomic_heterozygosity <- function(df,site_col, het_prefix) {
  # colonnes d'hétérozygotie
  het_cols <- grep(paste0("^", het_prefix), names(df), value = TRUE)
  df_clean <- df %>%
    mutate(across(all_of(het_cols),
                  ~ ifelse(. %in% c(-999, -1), NA, .)))
  # moyenne pondérée par individu
  res <- sapply(het_cols, function(ind) {
    sum(df_clean[[ind]] * df_clean[[site_col]], na.rm = TRUE) /
      sum(df_clean[[site_col]], na.rm = TRUE)
  })
  # nettoyage des noms
  individual_clean <- het_cols %>%
    gsub("^H\\.", "", .) %>%   # enlève "H."
    gsub("\\.$", "", .)       # enlève "." final
  
  data.frame(
    individual = individual_clean,
    H_genome = res
  )
}

##################
## 2. PATH.  #####
##################

files <- list(
  L_nobilis_POR = "/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/L_nobilis_POR/ALL_CHR_forward_reverse_no_indiv_with_N_renamed.out",
  L_nobilis_ITA = "/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/L_nobilis_ITA/ALL_CHR_forward_reverse_no_indiv_with_N_renamed.out",
  L_nobilis_TUN = "/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/L_nobilis_TUN/ALL_CHR_forward_reverse_no_indiv_with_N_renamed.out",
  L_nobilis_GRE = "/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/L_nobilis_GRE/ALL_CHR_forward_reverse_no_indiv_with_N_renamed.out",
  L_azorica = "/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/L_azorica/ALL_CHR_forward_reverse_no_indiv_with_N_renamed.out",
  L_novocanariensis = "/Users/cibio/Desktop/Laurus/Pop_gen/PiNPiS/RESULTS/L_novocanariensis/ALL_CHR_forward_reverse_no_indiv_with_N_renamed.out"
)

##################################
### DOWNLOAD DATA IN R ENV. ######
##################################

data_list <- setNames(lapply(names(files), function(pop) {
  df <- read_clean_pi_data(files[[pop]])
  df$population <- pop
  df
}),
names(files)
)

list2env(data_list, envir = .GlobalEnv)


#####################
### STAT PIN/PIS ####
#####################

results_PiNPiS <- lapply(names(files), function(pop) {
  df <- read_clean_pi_data(
    files[[pop]],
    pi_prefix = pop
  )
  stats <- compute_pi_stats(
    df,
    chr_remove = "OZ345845.1",
    piN_col = paste0(pop, "_piN"),
    piS_col = paste0(pop, "_piS"),
    site_col = paste0(pop, "_nb_complete_site")
  )
  
  stats$pop <- pop
  stats
})

results_PiNPiS_df <- bind_rows(results_PiNPiS)

#write_csv(results_PiNPiS_df, "results_PiNPiS_df.csv")

##########################
### BOOTSTRAP PIN/PIS ####
##########################

results_bootstrap_PiNPiS <- lapply(names(files), function(pop) {
  df <- read_clean_pi_data(
    files[[pop]],
    pi_prefix = pop
  )
  res <- bootstrap_pi(
    df,
    piN_col = paste0(pop, "_piN"),
    piS_col = paste0(pop, "_piS"),
    site_col = paste0(pop, "_nb_complete_site"),
    B = 1000
  )
  res$pop <- pop
  res
})

results_bootstrap_PiNPiS_df <- do.call(rbind, results_bootstrap_PiNPiS)
results_bootstrap_PiNPiS_df$pop <- gsub("^L_", "L.", results_bootstrap_PiNPiS_df$pop)
#write_csv(bootrap_PiNPiS_df, "bootrap_PiNPiS_df.csv")

###PLOT###

results_bootstrap_PiNPiS_df$pop <- factor( results_bootstrap_PiNPiS_df$pop, 
    levels = c("L.azorica","L.novocanariensis","L.nobilis_POR","L.nobilis_ITA","L.nobilis_TUN","L.nobilis_GRE"))

ggplot(results_bootstrap_PiNPiS_df, aes(x = pop, y = piS_mean, color = pop)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    ymin = piS_CI_low,
    ymax = piS_CI_high
  ), width = 0.15) +
  scale_color_manual(values = c(
    "L.azorica" = "#AC0C2F",
    "L.novocanariensis" = "#01558D",
    "L.nobilis_POR" = "#E6AF2E",
    "L.nobilis_ITA" = "#017939",
    "L.nobilis_TUN" = "#31AFD4",
    "L.nobilis_GRE" = "#9467BD"
  )) +
  scale_x_discrete(labels = c(
    "L.azorica" = expression(italic("L. azorica")),
    "L.novocanariensis" = expression(italic("L. novocanariensis")),
    "L.nobilis_POR" = expression(italic("L. nobilis") ~ "(POR)"),
    "L.nobilis_ITA" = expression(italic("L. nobilis") ~ "(ITA)"),
    "L.nobilis_TUN" = expression(italic("L. nobilis") ~ "(TUN)"),
    "L.nobilis_GRE" = expression(italic("L. nobilis") ~ "(GRE)")
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "Population",
    y = expression(pi[S])
  )


ggplot(results_bootstrap_PiNPiS_df, aes(x = pop, y = piN_mean, color = pop)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    ymin = piN_CI_low,
    ymax = piN_CI_high
  ), width = 0.15) +
  scale_color_manual(values = c(
    "L.azorica" = "#AC0C2F",
    "L.novocanariensis" = "#01558D",
    "L.nobilis_POR" = "#E6AF2E",
    "L.nobilis_ITA" = "#017939",
    "L.nobilis_TUN" = "#31AFD4",
    "L.nobilis_GRE" = "#9467BD"
  )) +
  scale_x_discrete(labels = c(
    "L.azorica" = expression(italic("L. azorica")),
    "L.novocanariensis" = expression(italic("L. novocanariensis")),
    "L.nobilis_POR" = expression(italic("L. nobilis") ~ "(POR)"),
    "L.nobilis_ITA" = expression(italic("L. nobilis") ~ "(ITA)"),
    "L.nobilis_TUN" = expression(italic("L. nobilis") ~ "(TUN)"),
    "L.nobilis_GRE" = expression(italic("L. nobilis") ~ "(GRE)")
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "Population",
    y = expression(pi[N])
  )
  



#####################
### STAT FIT ####
#####################

results_FIT <- lapply(names(files), function(pop) {
  df <- read_clean_pi_data(
    files[[pop]],
    pi_prefix = pop
  )
  stats <- compute_FIT_stats(
    df,
    chr_remove = "OZ345845.1",
    fit_col = paste0(pop, "_Fit"),
    wc_fit_col =  paste0(pop, "_WeirCockerham84_Fit"),
    site_col = paste0(pop, "_nb_complete_site")
  )
  stats$pop <- pop
  stats
})


results_FIT_df <- bind_rows(results_FIT)

#write_csv(results_FIT_df,"results_FIT_df.csv")

bootstrap_FIT <- function(df, fit_col, wc_fit_col, site_col, chr_remove = "OZ345845.1", B = 1000) {
  df_clean <- df %>%
    mutate(across(all_of(c(fit_col, wc_fit_col, site_col)), as.numeric)) %>%
    filter(
      !is.na(.data[[fit_col]]) &
        !is.na(.data[[wc_fit_col]]) &
        .data[[fit_col]] != -999 &
        .data[[wc_fit_col]] != -999
    ) %>%
    filter(CHR != chr_remove)
  n <- nrow(df_clean)
  boot_res <- replicate(B, {
    
    idx <- sample(1:n, replace = TRUE)
    d <- df_clean[idx, ]
    
    sum(d[[fit_col]] * d[[site_col]], na.rm = TRUE) /
      sum(d[[site_col]], na.rm = TRUE)
  })
  
  data.frame(
    FIT_mean = mean(boot_res, na.rm = TRUE),
    FIT_CI_low = quantile(boot_res, 0.025, na.rm = TRUE),
    FIT_CI_high = quantile(boot_res, 0.975, na.rm = TRUE)
  )
}

results_bootstrap_FIT <- lapply(names(files), function(pop) {
  
  df <- read_clean_pi_data(
    files[[pop]],
    pi_prefix = pop
  )
  
  res <- bootstrap_FIT(
    df,
    fit_col = paste0(pop, "_Fit"),
    wc_fit_col = paste0(pop, "_WeirCockerham84_Fit"),
    site_col = paste0(pop, "_nb_complete_site"),
    chr_remove = "OZ345845.1",
    B = 10
  )
  
  res$pop <- pop
  res
})

#write_csv(results_bootstrap_FIT_df,"results_bootstrap_FIT_df.csv")
results_bootstrap_FIT_df <- do.call(rbind, results_bootstrap_FIT)
results_bootstrap_FIT_df$pop <- gsub("^L_", "L.", results_bootstrap_FIT_df$pop)


### PLOT ###
results_bootstrap_FIT_df$pop <- factor( results_bootstrap_PiNPiS_df$pop, 
                                           levels = c("L.azorica","L.novocanariensis","L.nobilis_POR","L.nobilis_ITA","L.nobilis_TUN","L.nobilis_GRE"))

ggplot(results_bootstrap_FIT_df, aes(x = pop, y = FIT_mean, color = pop)) +
  geom_point(size = 2) +
  geom_errorbar(aes(
    ymin = FIT_CI_low,
    ymax = FIT_CI_high
  ), width = 0.15) +
  scale_color_manual(values = c(
    "L.azorica" = "#AC0C2F",
    "L.novocanariensis" = "#01558D",
    "L.nobilis_POR" = "#E6AF2E",
    "L.nobilis_ITA" = "#017939",
    "L.nobilis_TUN" = "#31AFD4",
    "L.nobilis_GRE" = "#9467BD"
  )) +
  scale_x_discrete(labels = c(
    "L.azorica" = expression(italic("L. azorica")),
    "L.novocanariensis" = expression(italic("L. novocanariensis")),
    "L.nobilis_POR" = expression(italic("L. nobilis") ~ "(POR)"),
    "L.nobilis_ITA" = expression(italic("L. nobilis") ~ "(ITA)"),
    "L.nobilis_TUN" = expression(italic("L. nobilis") ~ "(TUN)"),
    "L.nobilis_GRE" = expression(italic("L. nobilis") ~ "(GRE)")
  )) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "Population",
    y = expression(F[IT])
  )






###########################
####  HETEROZYGOTIE   #####
###########################

compute_genomic_heterozygosity <- function(df,site_col, het_prefix) {
  # colonnes d'hétérozygotie
  het_cols <- grep(paste0("^", het_prefix), names(df), value = TRUE)
  df_clean <- df %>%
    mutate(across(all_of(het_cols),
                  ~ ifelse(. %in% c(-999, -1), NA, .)))
  # moyenne pondérée par individu
  res <- sapply(het_cols, function(ind) {
    sum(df_clean[[ind]] * df_clean[[site_col]], na.rm = TRUE) /
      sum(df_clean[[site_col]], na.rm = TRUE)
  })
  # nettoyage des noms
  individual_clean <- het_cols %>%
    gsub("^H\\.", "", .) %>%   # enlève "H."
    gsub("\\.$", "", .)       # enlève "." final
  
  data.frame(
    individual = individual_clean,
    H_genome = res
  )
}

Ho_azorica <- compute_genomic_heterozygosity(
  df = L_azorica,
  site_col = "L_azorica_nb_complete_site",
  het_prefix = "H\\."
)

Ho_novocanariensis <- compute_genomic_heterozygosity(
  df = L_novocanariensis,
  site_col = "L_novocanariensis_nb_complete_site",
  het_prefix = "H\\."
)

Ho_L_nobilis_POR <- compute_genomic_heterozygosity(
  df = L_nobilis_POR,
  site_col = "L_nobilis_POR_nb_complete_site",
  het_prefix = "H\\."
)

Ho_L_nobilis_ITA <- compute_genomic_heterozygosity(
  df = L_nobilis_ITA,
  site_col = "L_nobilis_ITA_nb_complete_site",
  het_prefix = "H\\."
)

Ho_L_nobilis_TUN <- compute_genomic_heterozygosity(
  df = L_nobilis_TUN,
  site_col = "L_nobilis_TUN_nb_complete_site",
  het_prefix = "H\\."
)

Ho_L_nobilis_GRE <- compute_genomic_heterozygosity(
  df = L_nobilis_GRE,
  site_col = "L_nobilis_GRE_nb_complete_site",
  het_prefix = "H\\."
)


Ho_azorica$pop <- "L.azorica"
Ho_novocanariensis$pop <- "L.novocanariensis"
Ho_L_nobilis_POR$pop <- "L.nobilis_POR"
Ho_L_nobilis_ITA$pop <- "L.nobilis_ITA"
Ho_L_nobilis_TUN$pop <- "L.nobilis_TUN"
Ho_L_nobilis_GRE$pop <- "L.nobilis_GRE"

Ho_all <- bind_rows(
  Ho_azorica,
  Ho_novocanariensis,
  Ho_L_nobilis_POR,
  Ho_L_nobilis_ITA,
  Ho_L_nobilis_TUN,
  Ho_L_nobilis_GRE
)

#write_csv(Ho_all, "Observed_heterozygosity_Autosmes.csv")

## PLOT Ho ##

Ho_all$pop <- factor(
  Ho_all$pop,
  levels = c(
    "L.azorica",
    "L.novocanariensis",
    "L.nobilis_POR",
    "L.nobilis_ITA",
    "L.nobilis_TUN",
    "L.nobilis_GRE"
  )
)


ggplot(Ho_all, aes(x = pop, y = H_genome, fill = pop)) +
  geom_violin(trim = FALSE, alpha = 0.6) +
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.1, size = 1, alpha = 0.6) +
  
  scale_fill_manual(values = c(
    "L.azorica" = "#AC0C2F",
    "L.novocanariensis" = "#01558D",
    "L.nobilis_POR" = "#E6AF2E",
    "L.nobilis_ITA" = "#017939",
    "L.nobilis_TUN" = "#31AFD4",
    "L.nobilis_GRE" = "#9467BD"
  )) +
  
  scale_x_discrete(labels = c(
    "L.azorica" = expression(italic("L. azorica")),
    "L.novocanariensis" = expression(italic("L. novocanariensis")),
    "L.nobilis_POR" = expression(italic("L. nobilis") ~ "(POR)"),
    "L.nobilis_ITA" = expression(italic("L. nobilis") ~ "(ITA)"),
    "L.nobilis_TUN" = expression(italic("L. nobilis") ~ "(TUN)"),
    "L.nobilis_GRE" = expression(italic("L. nobilis") ~ "(GRE)")
  )) +
  
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "",
    y = expression("Genome-wide " ~ H[o])
  )
  





##########################
###.  PAR GENE.     ######
##########################


df_gene <- L_azorica_all_CHR %>%
  mutate(across(c(L_azorica_piN,
                  L_azorica_piS,
                  L_azorica_nb_complete_site),
                as.numeric)) %>%
  mutate(across(c(L_azorica_piN, L_azorica_piS), ~na_if(., -999))) %>%
  
  group_by(GENE, CHR) %>%
  summarise(
    Nsitegene = sum(L_azorica_nb_complete_site, na.rm = TRUE),
    
    piN = sum(L_azorica_piN * L_azorica_nb_complete_site, na.rm = TRUE) /
      sum(L_azorica_nb_complete_site, na.rm = TRUE),
    
    piS = sum(L_azorica_piS * L_azorica_nb_complete_site, na.rm = TRUE) /
      sum(L_azorica_nb_complete_site, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  filter(Nsitegene > 500) %>%
  
  mutate(piN_piS = piN / piS)

