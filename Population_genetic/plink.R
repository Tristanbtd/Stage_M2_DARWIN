
#################  SCRIPT FOR PLINK VIZUALISATION  #######################################

# Charger les packages nécessaires
library(ggplot2)
library(ggpubr)  # Pour ggarrange
library(readr)   # Pour read_delim
library(data.table)
library(qqman)
library(dplyr)
library(tidyr)

########################################################
###################### ACP #############################
########################################################

setwd("/Users/cibio/Desktop/Laurus/Pop_gen/PLINK/PCA")
# Définir le répertoire de travail

# Charger valeurs propres
eigenvalues <- scan("laurus_pca.eigenval")

# Variance expliquée (%)
variance_explained <- eigenvalues / sum(eigenvalues) * 100

variance_explained_df <- data.frame(
  PC = factor(paste0("PC", 1:length(variance_explained)),
              levels = paste0("PC", 1:length(variance_explained))),
  Variance = variance_explained
)

# Plot
ggplot(variance_explained_df, aes(x = PC, y = Variance)) +
  geom_bar(
    stat = "identity",
    fill = "grey70",        # gris doux
    color = "black",        # bordure noire
    linewidth = 0.5
  ) +
  labs(
    x = "Composantes principales",
    y = "Variance expliquée (%)")+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(face = "bold"),
    axis.line = element_line(color = "black"),
    #panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1
  ))



# Lire données
eigenvec <- read.table("laurus_pca.eigenvec", header = FALSE) %>%
  select(-V1)

colnames(eigenvec) <- c("Indiv", paste0("PC", 1:10))

# Espèce
eigenvec <- eigenvec %>%
  mutate(Espece = case_when(
    grepl("^La", Indiv) ~ "L. azorica",
    grepl("^LN", Indiv) ~ "L. novocanariensis",
    grepl("Lm", Indiv) ~ "L. nobilis (POR)",
    grepl("OR", Indiv) ~ "L. nobilis (ITA)",
    grepl("TUN", Indiv) ~ "L. nobilis (TUN)",
    grepl("GRE", Indiv) ~ "L. nobilis (GRE)"
  ))

eigenvec$Espece <- factor(
  eigenvec$Espece,
  levels = c(
    "L. azorica",
    "L. novocanariensis",
    "L. nobilis (POR)",
    "L. nobilis (ITA)",
    "L. nobilis (TUN)",
    "L. nobilis (GRE)"
  )
)
#########################
#### PLOT PC1 VS PC2 ####
#########################

ggplot(eigenvec, aes(PC1, PC2, color = Espece)) +
  geom_point(size = 3) +
  # ellipse (PROPRE et stable)
  #stat_ellipse(aes(group = Espece), linewidth = 0.8, linetype = 2) +
  # labels (IMPORTANT : ne pas trop restreindre)
  #geom_text_repel(
    #aes(label = Indiv),
    #size = 3,
    #max.overlaps = 20   # pour voir tous les labels)+
  scale_color_manual(
    values = c(
      "L. azorica" = "#01558D",
      "L. novocanariensis" = "#AC0C2F",
      "L. nobilis (POR)" = "#E6AF2E",
      "L. nobilis (ITA)" = "#017939",
      "L. nobilis (TUN)" = "#31AFD4",
      "L. nobilis (GRE)" = "#9467BD"),
    labels = c(
      expression(italic("L. azorica")),
      expression(italic("L. novocanariensis")),
      expression(italic("L. nobilis") ~ "(POR)"),
      expression(italic("L. nobilis") ~ "(ITA)"),
      expression(italic("L. nobilis") ~ "(TUN)"),
      expression(italic("L. nobilis") ~ "(GRE)") )
  ) +
  labs(
    x = "PC1 (21.12%)",
    y = "PC2 (13.15%)",
    title = "PCA of Laurus populations",
    color = "Species"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect( color = "black", fill = NA, linewidth = 1)
  )

#########################
#### PLOT PC2 VS PC3 ####
#########################
ggplot(eigenvec, aes(PC2, PC3, color = Espece)) +
  geom_point(size = 3) +
  # ellipse (PROPRE et stable)
  #stat_ellipse(aes(group = Espece), linewidth = 0.8, linetype = 2) +
  # labels (IMPORTANT : ne pas trop restreindre)
  #geom_text_repel(
    #aes(label = Indiv),
    #size = 3,
    #max.overlaps = 20   # pour voir tous les labels
  #) +
  scale_color_manual(
    values = c(
      "L. azorica" = "#01558D",
      "L. novocanariensis" = "#AC0C2F",
      "L. nobilis (POR)" = "#E6AF2E",
      "L. nobilis (ITA)" = "#017939",
      "L. nobilis (TUN)" = "#31AFD4",
      "L. nobilis (GRE)" = "#9467BD"),
    labels = c(
      expression(italic("L. azorica")),
      expression(italic("L. novocanariensis")),
      expression(italic("L. nobilis") ~ "(POR)"),
      expression(italic("L. nobilis") ~ "(ITA)"),
      expression(italic("L. nobilis") ~ "(TUN)"),
      expression(italic("L. nobilis") ~ "(GRE)") )
  ) +
  labs(
    x = "PC2 (13.15%)",
    y = "PC3 (12.56)",
    title = "PCA of Laurus populations",
    color = "Species"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.line = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.border = element_rect( color = "black", fill = NA, linewidth = 1)
  )

#########################
#### PLOT PC1 VS PC3 ####
#########################

ggplot(eigenvec, aes(PC1, PC3, color = Espece)) +
  geom_point(size = 3) +
  # ellipse (PROPRE et stable)
  #stat_ellipse(aes(group = Espece), linewidth = 0.8, linetype = 2) +
  # labels (IMPORTANT : ne pas trop restreindre)
  #geom_text_repel(
    #aes(label = Indiv),
    #size = 3,
    #max.overlaps = 20   # pour voir tous les labels
  #) +
  scale_color_manual(values = c(
    "L. azorica" = "#01558D",
    "L. novocanariensis" = "#AC0C2F",
    "L. nobilis (POR)" = "#E6AF2E",
    "L. nobilis (ITA)" = "#017939",
    "L. nobilis (TUN)" = "#31AFD4",
    "L. nobilis (GRE)" = "#9467BD"
  )) +
  labs(
    x = "PC1 (21.12%)",
    y = "PC3 (12.56)",
    title = "PCA of Laurus populations",
    color = "Species"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5)
  )

##### ACP 3D ######
library(plotly)
library(dplyr)

plot_ly(
  data = eigenvec,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  color = ~Espece,
  colors = c(
    "L. azorica" = "#AC0C2F",
    "L. novocanariensis" = "#01558D",
    "L. nobilis (POR)" = "#E6AF2E",
    "L. nobilis (ITA)" = "#017939",
    "L. nobilis (TUN)" = "#31AFD4",
    "L. nobilis (GRE)" = "#9467BD"
  ),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 4)
) %>%
  layout(
    title = "3D PCA of Laurus populations",
    scene = list(
      xaxis = list(title = "PC1"),
      yaxis = list(title = "PC2"),
      zaxis = list(title = "PC3")
    )
  )



####################################################
################# ADMIXTURE ########################
####################################################


# Set the working directory
setwd("/Users/cibio/Desktop/Laurus/Pop_gen/PLINK/ADMIXTURE")

# Read the sample list file
data_name <- read_delim("/Users/cibio/Desktop/Laurus/Pop_gen/PLINK/laurus_snps_final_LD_prunned.fam", col_names = FALSE, delim = "\t")
# Convert the sample list into a tibble
samplelist <- tibble(sample = data_name[[2]])


data_cross_val <- read_delim("cv_error.txt",col_names = FALSE)
# Remove suffixes from the names


# Initialize an empty tibble to store all data
all_data <- tibble(sample = character(),
                   k = numeric(),
                   Q = character(),
                   value = numeric())

# Loop through values of k
for (k in 4) {
  # Read the corresponding data file
  data <- read_delim(paste0("laurus_snps_final_LD_prunned.", k, ".Q"),
                     col_names = paste0("Q", seq(1, k)),
                     delim = " ")
  
  # Add the sample and k columns
  data <- data %>%
    mutate(sample = samplelist$sample,
           k = k)
  
  # Convert data from wide to long format
  data_long <- data %>%
    pivot_longer(cols = starts_with("Q"), names_to = "Q", values_to = "value")
  
  # Combine the current data with the overall data
  all_data <- bind_rows(all_data, data_long)
}

sample_order <- samplelist %>%
  mutate(order = case_when(
    grepl("^La", sample) ~ 1,
    grepl("^LN", sample) ~ 2,
    grepl("^Lm", sample) ~ 3,
    grepl("^OR", sample) ~ 4,
    grepl("^TUN|^TUNI", sample) ~ 5,
    grepl("^GRE", sample) ~ 6,
    TRUE ~ 7
  )) %>%
  arrange(order) %>%
  select(-order)              

# Convert sample to factor with specified order
all_data <- all_data %>%
  mutate(sample = factor(sample, levels = sample_order$sample))

all_data <- all_data %>%
  arrange(sample, k)

palette <- c(
  "#01558D",  # k=3 et k=4 Novocanariensis
  "#AC0C2F",  # k=3 et k=4 Azorica
  "#E6AF2E",  # k=3 Nobilis | k=4 Nobilis POR
  "#017939",  # NA | k=4 Nobilis IT/TUN/GRE
  "#796400",  # brun (121,100,0)
  "#D62728",  # rouge
  "#9467BD",  # violet
  "#8C564B",  # marron
  "#17BECF",  # cyan
  "#BCBD22"   # jaune-vert
)

#pour k=3
palette <- c(
  "#01558D",  # k=3 et k=4 Novocanariensis
  "#AC0C2F",  # k=3 et k=4 Azorica
  "#E6AF2E" # k=3 Nobilis | k=4 Nobilis POR
)

#pour k=4
palette <- c(
  "#AC0C2F",  # k=3 et k=4 Novocanariensis  # k=3 et k=4 Azorica  # k=3 Nobilis | k=4 Nobilis POR
  "#017939",
  "#01558D" ,
  "#E6AF2E"
)

all_data %>%
  ggplot(.,aes(x=sample,y=value,fill=factor(Q))) + 
  geom_bar(stat="identity",position="stack") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1 ,size=8)) +
  scale_fill_manual(values = palette, name = "K", labels = seq(1:5)) +
  facet_wrap(~k,ncol=1)

library(ggplot2)

ggplot(data_cross_val, aes(x = X1, y = X2)) +
  geom_point(color = "black", size = 3) +  # Points en bleu avec une taille de 3
  geom_line(color = "black", linetype = "dashed", size = 1) +  # Ligne rouge en pointillé
  labs( x = "K",
        y = "CV Error") +
  theme_minimal(base_size = 15) +  # Augmenter la taille de la police pour une meilleure lisibilité
  theme(
    axis.title = element_text(face = "bold"),  # Mettre en gras les titres des axes
    axis.text = element_text(color = "black"),  # Changer la couleur du texte des axes
    panel.border = element_rect(color = "black", fill = NA, size = 1),  # Ajouter un cadre autour du graphique
    panel.grid.major = element_line(color = "gray80"),  # Ajouter des lignes de grille principales
    panel.grid.minor = element_line(color = "gray80")   # Ajouter des lignes de grille mineures
  )



