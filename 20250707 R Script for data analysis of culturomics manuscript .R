

# Title: Targeted cultivation and isolation of probiotic candidate taxa from the cecal microbiota of broiler chickens using culturomics
# Date:2025-07-07
# Author: Manhong Wang
# Mississippi State University


# R Script for data analysis of "Targeted cultivation and isolation of probiotic candidate taxa from the cecal microbiota of broiler chickens using culturomics"




###Fig. 2A Venn diagrams showed the numbers of unique and shared bacterial species identified from CD and CI approaches


install.packages("eulerr")
install.packages("grid")
library(eulerr)
library(grid)

data <- c(
  CD = 190,
  CI = 342,
  "CD&CI" = 160
)

fit <- euler(data)

grid.newpage()
p <- plot(fit,
          fills = list(fill = c("green", "blue"), alpha = 0.5),
          edges = list(col = "black"),
          quantities = FALSE, 
          labels = FALSE) 

grid.draw(p)

grid.text("CD\n(Total: 350)", x = 0.10, y = 0.45, gp = gpar(fontsize = 18, fontface = "bold")) 
grid.text("CI\n(Total: 502)", x = 0.90, y = 0.45, gp = gpar(fontsize = 18, fontface = "bold")) 

grid.text("190", x = 0.25, y = 0.50,
          gp = gpar(fontsize = 18, fontface = "bold", family = "Times New Roman"))
grid.text("160", x = 0.45, y = 0.50,
          gp = gpar(fontsize = 18, fontface = "bold", family = "Times New Roman"))
grid.text("342", x = 0.68, y = 0.50,
          gp = gpar(fontsize = 18, fontface = "bold", family = "Times New Roman"))








###Fig. 2B Stacked bar plot showed the relative abundance of bacterial species identified by the CD and CI approaches


library(ggplot2)
library(dplyr)

df_species <- read.csv("Species abundance in CD&CI.csv")

all_combinations <- expand.grid(
  Group = unique(df_species$Group),
  Species = unique(df_species$Species)
)

df_species <- full_join(df_species, all_combinations, by = c("Group", "Species")) %>%
  mutate(Abundance = ifelse(is.na(Abundance), 0, Abundance))

df_species$Species <- factor(df_species$Species, levels = unique(df_species$Species))

df_species$Abundance <- ifelse(df_species$Abundance == 0, 0.0001, df_species$Abundance)

species_colors <- setNames(colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(length(unique(df_species$Species))),
                           unique(df_species$Species))
unique(df_species$Species)

highlight_species <- c(
  "Escherichia coli" = "#377EB8",
  "Proteus mirabilis" = "#4DAF4A",
  "Limosilactobacillus reuteri" = "#E41A1C",
  "Enterococcus faecalis" = "#FF7F00",
  "Ligilactobacillus salivarius" = "#984EA3",
  "Faecalibacterium prausnitzii" = "#A65628",
  "Eisenbergiella massiliensis" = "#FFFF33",
  "Negativibacillus massiliensis" = "#D95F02",
  "Lacrimispora saccharolytica" = "#7570B3",
  "Flintibacter sp. KGMB00164" = "#54278F",
  "Others<0.5% in CI" = "#E7298A",
  
  "Dysosmobacter welbionis" = "#1B9E77",
  "Fournierella massiliensis" = "#E6AB02",
  "Hungatella hathewayi" = "#A6761D",
  "Butyricicoccus pullicaecorum" = "#E7298A",
  "Lachnoclostridium sp. YL32" = "#66A61E",
  "Enterocloster bolteae" = "#E78AC3",
  "Mediterraneibacter glycyrrhizinilyticus" = "#A6CEE3",
  "Blautia sp. YL58" = "#FC8D62",
  "Ruthenibacterium lactatiformans" = "#8DA0CB",
  "Intestinimonas butyriciproducens" = "#B3DE69",
  "Paludicola psychrotolerans" = "#FFD92F",
  "Pseudoflavonifractor capillosus" = "#E5C494",
  "Lachnoclostridium phocaense" = "#A6D854",
  "Anaerostipes butyraticus" = "#FDB462",
  "Gemmiger formicilis" = "#BEAED4",
  "Pseudoclostridium thermosuccinogenes" = "#CAB2D6",
  "Faecalicatena fissicatena" = "#BC80BD",
  "Anaeromassilibacillus senegalensis" = "#FB9A99",
  "Eisenbergiella tayi" = "#FDBF6F",
  "Faecalinomas umbilicata" = "#FFED6F",
  "Papillibacter cinnamivorans" = "#F0027F",
  "Clostridium sp. BNL1100" = "#6A3D9A",
  "Agathobaculum butyriciproducens" = "#B15928",
  "Neglecta timonensis" = "#FF7F50",
  "Roseburia hominis" = "#40E0D0",
  "Anaerotruncus rubinifantis" = "#BC8F8F",
  "Fusicatenibacter saccharivorans" = "#FF1493"
)

other_species <- setdiff(unique(df_species$Species), names(highlight_species))
other_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(other_species))

species_colors <- c(highlight_species, setNames(other_colors, other_species))

print(highlight_species["Limosilactobacillus reuteri"])

p <- ggplot(df_species, aes(x = Group, y = Abundance, fill = Species)) +
  geom_bar(stat = "identity", position = "fill") +  # Stacked Bar Plot
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "25", "50", "75", "100")
  ) +
  scale_fill_manual(values = species_colors, drop = FALSE) +
  theme_minimal(base_family = "Times New Roman") +  
  labs(x = "", y = "Relative Abundance (%)", fill = "Species") +
  theme(
    axis.text.x = element_text(
      size = 14, face = "bold", family = "Times New Roman",
      margin = margin(t = 0)  
    ),
    axis.text.y = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title.y = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.title = element_text(size = 12, face = "bold", family = "Times New Roman"),
    legend.text = element_text(size = 11, face = "bold.italic", family = "Times New Roman"),
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right",
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.2, "cm"),
    plot.title = element_text(face = "bold", family = "Times New Roman", hjust = 0.5),
    plot.margin = margin(2, 10, 2, 10)  
  )+
  guides(fill = guide_legend(ncol = 1))


print(p)




###Fig.3A Alpha diversity indices of bacterial communities across different culture media 


install.packages("readr")
insinstall.packages("vegan")  
install.packages("tidyr")

library(vegan)             
library(readr)
library(ggplot2)
library(dplyr)
library(FSA)        
library(rstatix)     
library(ggsignif)   
library(ggtext)
library(tidyr)


otu <- read.csv("CD abundance table.csv", row.names = 1, check.names = FALSE)   
meta <- read.csv("Metadata CD.csv", row.names = 1, check.names = FALSE)         


otu_t <- as.data.frame(t(otu))


chao1    <- estimateR(otu_t)["S.chao1", ]               
shannon  <- diversity(otu_t, index = "shannon")        
simpson  <- diversity(otu_t, index = "simpson")        

alpha_df <- data.frame(
  SampleID = rownames(otu_t),
  Chao1 = chao1,
  Shannon = shannon,
  Simpson = simpson
) %>%
  pivot_longer(cols = -SampleID, names_to = "Measure", values_to = "Value")

meta$SampleID <- rownames(meta)
alpha_merged <- left_join(alpha_df, meta, by = "SampleID")

alpha_res <- list()
alpha_res$data_alpha <- alpha_merged

selected_measures <- c("Chao1", "Shannon", "Simpson")

plot_data <- alpha_res$data_alpha %>%
  filter(Measure %in% selected_measures)

annotation_all <- list()
kw_pvalues <- c()

for (index in selected_measures) {
  subset_data <- plot_data %>% filter(Measure == index)
  
  
  kw_res <- kruskal.test(Value ~ Medium, data = subset_data)
  kw_pvalues[index] <- format.pval(kw_res$p.value, digits = 3, eps = 0.001)
  
 
  dunn_res <- dunnTest(Value ~ Medium, data = subset_data, method = "bonferroni")$res
 
  sig_res <- dunn_res %>%
    filter(P.adj < 0.05) %>%
    mutate(
      group1 = sub(" - .*", "", Comparison),
      group2 = sub(".*- ", "", Comparison),
      p_signif = case_when(
        P.adj < 0.001 ~ "***",
        P.adj < 0.01 ~ "**",
        P.adj < 0.05 ~ "*"
      ),
      Measure = index
    )
  
  
  if (nrow(sig_res) > 0) {
    y_max <- max(subset_data$Value, na.rm = TRUE)
    sig_res$y.position <- seq(y_max * 1.05, by = y_max * 0.1, length.out = nrow(sig_res))
    annotation_all[[index]] <- sig_res
  }
}


annotation_df <- do.call(rbind, annotation_all)


kw_labels <- mapply(function(name, value) {
  if (grepl("^<", value)) {
    paste0(name, "<br>KW <i>P</i> ", value)
  } else {
    paste0(name, "<br>KW <i>P</i> = ", value)
  }
}, name = names(kw_pvalues), value = kw_pvalues)

measure_labels <- setNames(kw_labels, names(kw_pvalues))

p <- ggplot(plot_data, aes(x = Medium, y = Value, fill = Medium)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(
    ~ Measure,
    scales = "free_y",
    labeller = labeller(Measure = measure_labels)
  ) +
  theme_minimal(base_family = "Times New Roman") +  
  labs(
    title = "Alpha Diversity by Medium with Dunn Test",
    x = "Medium",
    y = "Diversity Index"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90, hjust = 1, face = "bold", size = 12,
      family = "Times New Roman"
    ),
    axis.text.y = element_text(
      face = "bold", size = 12,
      family = "Times New Roman"
    ),
    axis.title.x = element_text(
      face = "bold", size = 14,
      margin = margin(t = 10, b = 20),
      family = "Times New Roman"
    ),
    axis.title.y = element_text(
      face = "bold", size = 14,
      family = "Times New Roman"
    ),
    strip.text = element_markdown(
      face = "bold", size = 12,
      family = "Times New Roman"
    ),
    strip.background = element_rect(fill = "#D9D9D9", color = NA),
    plot.title = element_text(
      hjust = 0.5, size = 16, face = "bold",
      margin = margin(b = 14),
      family = "Times New Roman"
    ),
    legend.position = "none"
  )

if (!is.null(annotation_df)) {
  p <- p + geom_signif(
    data = annotation_df,
    aes(xmin = group1, xmax = group2, annotations = p_signif, y_position = y.position),
    manual = TRUE,
    tip_length = 0.01,
    textsize = 4.2,
    inherit.aes = FALSE
  )
}

print(p)

#output Chao1", "Shannon", "Simpson" data as CSV file

library(dplyr)
library(FSA)  # dunnTest
selected_measures <- c("Chao1", "Shannon", "Simpson")

all_dunn_results <- list()
for (index in selected_measures) {
  subset_data <- plot_data %>% filter(Measure == index)
  
  # Dunn Test with Bonferroni correction
  dunn_res <- dunnTest(Value ~ Medium, data = subset_data, method = "bonferroni")$res
  dunn_res$Measure <- index
  
  all_dunn_results[[index]] <- dunn_res
}
combined_dunn_df <- bind_rows(all_dunn_results)
write.csv(combined_dunn_df, "dunn_pairwise_comparisons.csv", row.names = FALSE)
getwd()



### Fig.3B Principal coordinate analysis (PCoA) based on Bray–Curtis dissimilarity showing microbial community differences by culture medium 


install.packages(c("vegan", "ggtext", "tidyverse"))

library(vegan)
library(ggtext)
library(tidyverse)


otu <- read.csv("CD abundance table.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Metadata CD.csv", row.names = 1, check.names = FALSE)

otu_t <- t(otu)
otu_t <- otu_t[rownames(meta), ]


bray_dist <- vegdist(otu_t, method = "bray")


adonis_result <- adonis2(bray_dist ~ Medium, data = meta)
adonis_r2 <- adonis_result$R2[1]
adonis_p <- adonis_result$`Pr(>F)`[1]

adonis_r2_fmt <- format(round(adonis_r2, 4), nsmall = 4)
adonis_p_fmt <- ifelse(adonis_p < 0.001, "< 0.001", format.pval(adonis_p, digits = 3))
adonis_title <- paste0(
  "PERMANOVA: R² = ", adonis_r2_fmt, ", <i>P</i> = ", adonis_p_fmt
)


pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- cbind(pcoa_df, meta[rownames(pcoa_df), ])


var_explained <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)


ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Medium)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = adonis_title,
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  ) +
  scale_color_manual(values = c(
    "BBE" = "#D55E00", "BHI" = "#E69F00", "CBA" = "#56B4E9",
    "CBAB" = "#009E73", "CNA" = "#F0E442", "CNAB" = "#0072B2",
    "FAA" = "#CC79A7", "KVLB" = "#999999", "M9" = "#999933",
    "M9I" = "#6699CC", "MAC" = "#FF9966", "MRS" = "#33CC99",
    "MSA" = "#9966CC", "TSY" = "#FF6666"
  )) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = ggtext::element_markdown(
      hjust = 0.5, face = "bold", size = 16,
      margin = margin(b = 10), family = "Times New Roman"
    ),
    legend.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    legend.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    axis.text = element_text(size = 14, face = "bold", family = "Times New Roman")
  )




### Fig.3C Alpha diversity indices of bacterial communities across different air conditions


install.packages(c("vegan", "FSA", "rstatix", "ggsignif", "ggtext"))

library(vegan)
library(FSA)
library(rstatix)
library(ggsignif)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)


otu_raw <- read.csv("CD abundance table.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Metadata CD.csv", row.names = 1, check.names = FALSE)


otu_t <- as.data.frame(t(otu_raw))
otu_t <- otu_t[rownames(meta), ]  

chao1 <- estimateR(otu_t)["S.chao1", ]
shannon <- diversity(otu_t, index = "shannon")
simpson <- diversity(otu_t, index = "simpson")


data_alpha <- data.frame(
  SampleID = rownames(otu_t),
  Chao1 = chao1,
  Shannon = shannon,
  Simpson = simpson
) %>%
  pivot_longer(cols = -SampleID, names_to = "Measure", values_to = "Value")


meta$SampleID <- rownames(meta)
plot_data <- left_join(data_alpha, meta, by = "SampleID")


selected_measures <- c("Chao1", "Shannon", "Simpson")
plot_data <- plot_data %>% filter(Measure %in% selected_measures)


annotation_all <- list()
kw_pvalues <- c()

for (index in selected_measures) {
  subset_data <- plot_data %>% filter(Measure == index)
  
 
  kw_res <- kruskal.test(Value ~ `Air condition`, data = subset_data)
  kw_pvalues[index] <- format.pval(kw_res$p.value, digits = 3, eps = 0.001)
  
  dunn_res <- dunnTest(Value ~ `Air condition`, data = subset_data, method = "bonferroni")$res
  
  sig_res <- dunn_res %>%
    filter(P.adj < 0.05) %>%
    mutate(
      group1 = sub(" - .*", "", Comparison),
      group2 = sub(".*- ", "", Comparison),
      p_signif = case_when(
        P.adj < 0.001 ~ "***",
        P.adj < 0.01 ~ "**",
        P.adj < 0.05 ~ "*"
      ),
      Measure = index
    )
  
  if (nrow(sig_res) > 0) {
    y_max <- max(subset_data$Value, na.rm = TRUE)
    sig_res$y.position <- seq(y_max * 1.05, by = y_max * 0.1, length.out = nrow(sig_res))
    annotation_all[[index]] <- sig_res
  }
}

annotation_df <- do.call(rbind, annotation_all)


kw_labels <- mapply(function(name, value) {
  if (grepl("^<", value)) {
    paste0(name, "<br>KW <i>P</i> ", value)
  } else {
    paste0(name, "<br>KW <i>P</i> = ", value)
  }
}, name = names(kw_pvalues), value = kw_pvalues)

measure_labels <- setNames(kw_labels, names(kw_pvalues))


p <- ggplot(plot_data, aes(x = `Air condition`, y = Value, fill = `Air condition`)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Measure, scales = "free_y", labeller = labeller(Measure = measure_labels)) +
  theme_minimal(base_family = "Times New Roman") +
  labs(
    title = "Alpha Diversity by Air condition with Dunn Test",
    x = "Air condition",
    y = "Diversity Index"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, face = "bold", size = 12, family = "Times New Roman"),
    axis.text.y = element_text(face = "bold", size = 12, family = "Times New Roman"),
    axis.title.x = element_text(face = "bold", size = 14, margin = margin(t = 10, b = 20), family = "Times New Roman"),
    axis.title.y = element_text(face = "bold", size = 14, family = "Times New Roman"),
    strip.text = ggtext::element_markdown(face = "bold", size = 12, family = "Times New Roman"),
    strip.background = element_rect(fill = "#D9D9D9", color = NA),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 14), family = "Times New Roman"),
    legend.position = "none"
  )


if (!is.null(annotation_df) && nrow(annotation_df) > 0) {
  p <- p + geom_signif(
    data = annotation_df,
    aes(xmin = group1, xmax = group2, annotations = p_signif, y_position = y.position),
    manual = TRUE,
    tip_length = 0.01,
    textsize = 4.2,
    inherit.aes = FALSE
  )
}


print(p)



#output Chao1", "Shannon", "Simpson" data as CSV file

library(dplyr)
library(FSA) 

selected_measures <- c("Chao1", "Shannon", "Simpson")

all_dunn_results <- list()

for (index in selected_measures) {
  subset_data <- plot_data %>% filter(Measure == index)
  
  
  dunn_res <- dunnTest(Value ~ Aircondition, data = subset_data, method = "bonferroni")$res
  
  dunn_res$Measure <- index
  
  all_dunn_results[[index]] <- dunn_res
}

combined_dunn_df <- bind_rows(all_dunn_results)

write.csv(combined_dunn_df, "dunn_pairwise_comparisons_aircondition.csv", row.names = FALSE)

getwd()







### Fig.3D Principal coordinate analysis (PCoA) based on Bray–Curtis dissimilarity showing microbial community differences by air condition



install.packages(c("vegan", "ggplot2", "ggtext"))
library(vegan)
library(ggplot2)
library(ggtext)

otu <- read.csv("CD abundance table.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Metadata CD.csv", row.names = 1, check.names = FALSE)

otu_t <- t(otu)
otu_t <- otu_t[rownames(meta), ]

bray_dist <- vegdist(otu_t, method = "bray")

pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- cbind(pcoa_df, meta[rownames(pcoa_df), ])


adonis_result <- adonis2(bray_dist ~ `Air condition`, data = meta)
adonis_r2 <- adonis_result$R2[1]
adonis_p <- adonis_result$`Pr(>F)`[1]


adonis_r2_formatted <- format(round(adonis_r2, 4), nsmall = 4)
adonis_p_formatted <- ifelse(adonis_p < 0.001, "< 0.001", format.pval(adonis_p, digits = 3))
adonis_caption <- paste0(
  "PCoA of Bray-Curtis Distance by Air condition<br>",
  "PERMANOVA: R² = ", adonis_r2_formatted, ", <i>P</i> = ", adonis_p_formatted
)


var_explained <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)


ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = `Air condition`)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(
    title = adonis_caption,
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  ) +
  scale_color_manual(values = c("AE" = "#E64B35FF", "AN" = "#4DBBD5FF")) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = ggtext::element_markdown(
      hjust = 0.5,
      face = "bold",
      size = 16,
      lineheight = 1.4,
      margin = margin(b = 10),
      family = "Times New Roman"
    ),
    axis.title.x = element_text(face = "bold", size = 14, family = "Times New Roman"),
    axis.title.y = element_text(face = "bold", size = 14, family = "Times New Roman"),
    axis.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.title = element_text(face = "bold", size = 16, family = "Times New Roman"),
    legend.text = element_text(size = 14, face = "bold", family = "Times New Roman")
  )



###Fig.3E Stacked bar plot showed the relative abundance of bacterial genera across 14 media under both AE and AN conditions



install.packages(c("ggplot2", "RColorBrewer", "randomcoloR", "dplyr", "ggtext"))
library(ggplot2)
library(RColorBrewer)
library(randomcoloR)
library(dplyr)
library(ggtext)
library(stringr)
library(grid)

df <- read.csv("Genus abundance under 28 conditions.csv", na.strings = c("", "NA"))


df$Abundance <- as.numeric(gsub("%", "", df$Abundance))

df$Genus <- df$Genus %>%
  str_trim() %>%
  str_replace_all("[^[:print:]]", "") %>%
  str_to_title()
df$Genus <- factor(df$Genus, levels = unique(df$Genus))  

df$Medium_Condition <- interaction(df$Medium, df$Condition, sep = " ")

custom_colors <- c(
  "Escherichia" = "#4B44F2", "Proteus" = "#008000", "Enterococcus" = "#984EA3",
  "Lactobacillus" = "#FF0000", "Limosilactobacillus" = "#FF9999", "Ligilactobacillus" = "#ff5500",
  "Shigella" = "#8B4513", "Laceyella" = "#66C2A5", "Enterobacter" = "#FC8D62",
  "Xenorhabdus" = "#8DA0CB", "Other genus>0.05%" = "#E78AC3", "Klebsiella" = "#B3B3B3",
  "Bacillus" = "#FFFF33", "Faecalibacterium" = "#E5C494", "Pediococcus" = "#B3B3B3",
  "Paenibacillus" = "#E41A1C", "Bacteroides" = "#993c00", "Acinetobacter" = "#4DAF4A",
  "Blautia" = "#FFD700", "Lacrimispora" = "#FF7F00", "Eubacterium" = "#FFFF33",
  "Pseudoflavonifractor" = "#A65628", "Anaerotignum" = "#F781BF", "Flintibacter" = "#999999",
  "Mediterraneibacter" = "#66C2A5", "Massilimicrobiota" = "#FC8D62", "Flavonifractor" = "#8DA0CB",
  "Drancourtella" = "#E78AC3", "Fournierella" = "#A6D854", "Staphylococcus" = "#FFD92F",
  "Sellimonas" = "#E5C494", "Lachnoclostridium" = "#A6D854", "Anaerotruncus" = "#E41A1C",
  "Enterocloster" = "#377EB8", "Faecalicatena" = "#4DAF4A", "Erysipelatoclostridium" = "#984EA3",
  "Massilistercora" = "#FF7F00", "Butyricicoccus" = "#FFFF33", "Alistipes" = "#A65628",
  "Agathobaculum" = "#F781BF"
)

missing_colors <- setdiff(levels(df$Genus), names(custom_colors))
if (length(missing_colors) > 0) {
  extra_colors <- distinctColorPalette(length(missing_colors))
  names(extra_colors) <- missing_colors
  custom_colors <- c(custom_colors, extra_colors)
}


levels(df$Genus) <- paste0("*", levels(df$Genus), "*")

names(custom_colors) <- paste0("*", names(custom_colors), "*")

p <- ggplot(df, aes(x = Medium, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title = "Bacterial community composition at genus level",
    x = "Medium",
    y = "Relative Abundance (%)"
  ) +
  scale_fill_manual(values = custom_colors) +
  facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    axis.text.x = element_text(
      angle = 90, vjust = 0.5, hjust = 1,
      size = 14, face = "bold", family = "Times New Roman"
    ),
    axis.text.y = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title.x = element_text(size = 16, face = "bold", family = "Times New Roman"),
    axis.title.y = element_text(size = 16, face = "bold", family = "Times New Roman"),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 20), family = "Times New Roman"),
    legend.position = "bottom", 
    legend.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.text = element_markdown(size = 12, face = "bold", family = "Times New Roman"),  
    legend.key.size = unit(0.5, "lines"),
    legend.box = "horizontal", 
    plot.margin = margin(t = 1, r = 2, b = 1, l = 1, unit = "cm"),
    panel.spacing = unit(2, "lines"),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text = element_text(size = 14, face = "bold", family = "Times New Roman")
  ) +
  guides(fill = guide_legend(ncol = 5))  

print(p)





### Fig.4A PCoA based on Bray–Curtis dissimilarity showing community structure across five selected culture media



library(vegan)
library(tidyverse)
library(ggtext)  

otu <- read.csv("CD abundance table.csv", row.names = 1, check.names = FALSE)
meta <- read.csv("Metadata CD.csv", row.names = 1, check.names = FALSE)


selected_mediums <- c("BBE", "CNA", "CNAB", "MRS", "MSA")
meta_filtered <- meta[meta$Medium %in% selected_mediums, ]

otu_t <- t(otu)
otu_t <- otu_t[rownames(meta_filtered), ]
bray_dist <- vegdist(otu_t, method = "bray")

adonis_result <- adonis2(bray_dist ~ Medium, data = meta_filtered)

adonis_r2 <- adonis_result$R2[1]
adonis_p <- adonis_result$`Pr(>F)`[1]
adonis_r2_formatted <- format(round(adonis_r2, 4), nsmall = 4)
adonis_p_formatted <- ifelse(adonis_p < 0.001, "&lt; 0.001", format.pval(adonis_p, digits = 3))


adonis_caption <- paste0(
  "PCoA of Bray-Curtis Distance (Selected Media)<br>",
  "PERMANOVA: R² = ", adonis_r2_formatted, ", <i>P</i> = ", adonis_p_formatted
)


pcoa_res <- cmdscale(bray_dist, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa_res$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$SampleID <- rownames(pcoa_df)
pcoa_df <- cbind(pcoa_df, meta_filtered[rownames(pcoa_df), ])


var_explained <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 1)


ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = Medium)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(
    aes(group = Medium),
    linetype = "dashed",
    size = 1,
    type = "t",
    level = 0.95
  ) +
  labs(
    title = adonis_caption,
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  ) +
  scale_color_manual(values = c(
    "BBE" = "#D55E00", "CNA" = "#F0E442", "CNAB" = "#0072B2",
    "MRS" = "#33CC99", "MSA" = "#9966CC"
  )) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    plot.title = ggtext::element_markdown(
      hjust = 0.5, face = "bold", size = 16, lineheight = 1.3, margin = margin(b = 10)
    ),
    legend.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    legend.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title = element_text(face = "bold", size = 14, family = "Times New Roman"),
    axis.text = element_text(size = 14, face = "bold", family = "Times New Roman")
  )





### Fig.4B UpSet plot showing the number of shared and unique species among the five media


library(tidyverse)
library(readr)
library(UpSetR)

metadata <- read_csv("Metadata CD.csv")
abundance <- read_csv("CD abundance table.csv")

abundance <- as.data.frame(abundance)
rownames(abundance) <- abundance$`#name`
abundance <- abundance[, -1]

abundance_rel <- sweep(abundance, 2, colSums(abundance), FUN = "/")

abundance_long <- abundance_rel %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Barcode", values_to = "RelAbundance")

metadata <- metadata %>% rename(Barcode = `#name`)
merged_df <- left_join(abundance_long, metadata, by = "Barcode")

selected_mediums <- c("MRS", "CNAB", "CNA", "BBE", "MSA")
merged_df <- merged_df %>%
  filter(Medium %in% selected_mediums)


colnames(merged_df) <- gsub("Aircondition", "Air condition", colnames(merged_df))

abundance_summary <- merged_df %>%
  group_by(Species, Medium, `Air condition`) %>%
  summarise(mean_abundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")

abundance_filtered <- abundance_summary %>%
  group_by(Species) %>%
  filter(any(mean_abundance >= 0.0001)) %>%
  ungroup()

top_species <- abundance_filtered %>%
  group_by(Species) %>%
  summarise(avg_abund = mean(mean_abundance)) %>%
  slice_max(avg_abund, n = 50) %>%
  pull(Species)

abundance_top50 <- abundance_summary %>%
  filter(Species %in% top_species)

otu_bin <- abundance_top50 %>%
  filter(mean_abundance >= 0.0001) %>%
  distinct(Species, Medium) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = Medium, values_from = value, values_fill = 0) %>%
  as.data.frame()

rownames(otu_bin) <- otu_bin$Species
otu_bin$Species <- NULL


write.csv(otu_bin, "Top50_PresenceAbsence_5Media_MergedAEAN.csv")

UpSetR::upset(
  otu_bin,
  sets = c("MRS", "CNAB", "CNA", "BBE", "MSA"),
  sets.bar.color = "#3CA370",
  keep.order = TRUE,
  point.size = 3,
  line.size = 0.7,
  order.by = "freq",
  mainbar.y.label = "Number of Species",
  sets.x.label = "Number of Species/medium",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(2.0, 2.0, 2.0, 2.0, 2.0, 2.0)  
)

colSums(otu_bin)






### Fig.4C UpSet plot showing the number of shared and unique species among air condition

library(tidyverse)
library(readr)
library(UpSetR)

metadata <- read_csv("Metadata CD.csv")
abundance <- read_csv("CD abundance table.csv")

abundance <- as.data.frame(abundance)
rownames(abundance) <- abundance$`#name`
abundance <- abundance[, -1]

abundance_rel <- sweep(abundance, 2, colSums(abundance), FUN = "/")


abundance_long <- abundance_rel %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Barcode", values_to = "RelAbundance")

metadata <- metadata %>% rename(Barcode = `#name`)
merged_df <- left_join(abundance_long, metadata, by = "Barcode")


selected_mediums <- c("MRS", "CNAB", "CNA", "BBE", "MSA")
merged_df <- merged_df %>%
  filter(Medium %in% selected_mediums)

colnames(merged_df) <- gsub("Aircondition", "Air condition", colnames(merged_df))


abundance_summary <- merged_df %>%
  group_by(Species, Medium, `Air condition`) %>%
  summarise(mean_abundance = mean(RelAbundance), .groups = "drop")


abundance_filtered <- abundance_summary %>%
  group_by(Species) %>%
  filter(any(mean_abundance >= 0.0001)) %>%
  ungroup()

top_species <- abundance_filtered %>%
  group_by(Species) %>%
  summarise(avg_abund = mean(mean_abundance)) %>%
  slice_max(avg_abund, n = 50) %>%
  pull(Species)

aean_summary <- merged_df %>%
  filter(Species %in% top_species) %>%
  group_by(Species, `Air condition`) %>%
  summarise(mean_abundance = mean(RelAbundance), .groups = "drop")

otu_bin <- aean_summary %>%
  filter(mean_abundance >= 0.0001) %>%
  select(Species, `Air condition`) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = `Air condition`, values_from = value, values_fill = 0) %>%
  distinct(Species, .keep_all = TRUE) %>%
  as.data.frame()

rownames(otu_bin) <- otu_bin$Species
otu_bin$Species <- NULL

write.csv(otu_bin, "Top50_PresenceAbsence_AEAN_Revised.csv")

png("upset_species_condition.png", width = 3500, height = 2600, res = 300)

upset(
  otu_bin,
  sets = c("AE", "AN"),
  sets.bar.color = "#3CA370",
  keep.order = TRUE,
  point.size = 3,
  line.size = 0.7,
  order.by = "freq",
  mainbar.y.label = "Number of Species",
  sets.x.label = "Number of Species/Air condition",
  mb.ratio = c(0.6, 0.4),
  text.scale = c(2.0, 2.0, 2.0, 2.0, 2.0, 3.0)  
)

colSums(otu_bin)






### Fig.4D Heatmap revealed distinct enrichment patterns among 10 different culture conditions (5 medium types x 2 air conditions).

library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(circlize)
metadata <- read_csv("Metadata CD.csv")
abundance <- read_csv("CD abundance table.csv")

abundance <- as.data.frame(abundance)
rownames(abundance) <- abundance$`#name`
abundance <- abundance[, -1]

abundance_rel <- sweep(abundance, 2, colSums(abundance), FUN = "/")

abundance_long <- abundance_rel %>%
  rownames_to_column("Species") %>%
  pivot_longer(-Species, names_to = "Barcode", values_to = "RelAbundance")

metadata <- metadata %>% rename(Barcode = `#name`)
merged_df <- left_join(abundance_long, metadata, by = "Barcode")

selected_mediums <- c("MRS", "CNAB", "CNA", "BBE", "MSA")
merged_df <- merged_df %>%
  filter(Medium %in% selected_mediums)

colnames(merged_df) <- gsub("Aircondition", "Air condition", colnames(merged_df))


abundance_summary <- merged_df %>%
  group_by(Species, Medium, `Air condition`) %>%
  summarise(mean_abundance = mean(RelAbundance, na.rm = TRUE), .groups = "drop")


abundance_filtered <- abundance_summary %>%
  group_by(Species) %>%
  filter(any(mean_abundance >= 0.0001)) %>%
  ungroup()


top_species <- abundance_filtered %>%
  group_by(Species) %>%
  summarise(avg_abund = mean(mean_abundance)) %>%
  slice_max(avg_abund, n = 50) %>%
  pull(Species)


abundance_top <- abundance_summary %>%
  filter(Species %in% top_species) %>%
  unite("Medium_Condition", Medium, `Air condition`, sep = "_") %>%
  pivot_wider(names_from = Medium_Condition, values_from = mean_abundance, values_fill = 0) %>%
  column_to_rownames("Species")


write.csv(abundance_top, "Top50_RelAbundance_5Media.csv")


log_mat <- log10(as.matrix(abundance_top) + 1e-6)

col_fun <- colorRamp2(
  c(-6, -3, -1, 0),
  c("ivory", "#abd9e9", "#fdae61", "#d7191c")
)


column_groups <- ifelse(grepl("_AE$", colnames(log_mat)), "AE", "AN")
column_order <- colnames(log_mat)[order(column_groups, colnames(log_mat))]
log_mat <- log_mat[, column_order]
column_groups <- ifelse(grepl("_AE$", colnames(log_mat)), "AE", "AN")


library(ComplexHeatmap)
library(circlize)
library(grid)

ht <- Heatmap(
  log_mat,
  name = "log10 (abundance + 1e-6)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  column_split = factor(column_groups, levels = c("AN", "AE")),
  show_column_names = TRUE,
  column_names_rot = 90,
  row_dend_width = unit(3, "cm"),
  
  row_names_gp = gpar(fontsize = 9, fontface = "bold.italic", fontfamily = "Times New Roman"),
  column_names_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Times New Roman"),
  
  column_title_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Times New Roman"),
 
  heatmap_legend_param = list(
    title_gp = gpar(fontface = "bold", fontfamily = "Times New Roman"),
    labels_gp = gpar(fontsize = 12, fontface = "bold", fontfamily = "Times New Roman")
  ),
  

  layer_fun = function(j, i, x, y, w, h, fill) {
    grid.rect(x = x, y = y, width = w, height = h,
              gp = gpar(col = "white", lwd = 0.5, fill = NA))
  }
)

draw(
  ht,
  heatmap_legend_side = "right",
  padding = unit(c(10, 10, 30, 10), "mm")  
)

png("Top50_Heatmap_SpeciesItalic.png", width = 3000, height = 2600, res = 300)
draw(ht, heatmap_legend_side = "right", padding = unit(c(8, 8, 15, 8), "mm"))
dev.off()








### Fig.5 Sankey diagram for 20 beneficial species across 5 medium and 2 air conditions

library(tidyverse)
library(ggalluvial)
library(scales)

df <- read.csv("20 beneficial species for sankey.csv", check.names = FALSE)

df_long <- df %>%
  pivot_longer(
    cols = -Species,
    names_to = "Condition",
    values_to = "Abundance"
  )

df_long <- df_long %>%
  mutate(
    Abundance = as.numeric(gsub("%", "", Abundance)),
    Aircondition = sub(".*_", "", Condition),
    Medium = sub("_.*", "", Condition)
  )

df_long <- df_long %>% filter(Abundance > 0)

species_colors <- c(
  "Limosilactobacillus reuteri" = "#E41A1C",
  "Enterococcus faecalis" = "#4DAF4A",
  "Ligilactobacillus salivarius" = "#984EA3",
  "Enterococcus cecorum" = "#FFFF33",
  "Enterococcus sp. FDAARGOS_375" = "#FF9900",
  "Lactobacillus johnsonii" = "#f4cccc",
  "Bacteroides fragilis" = "#993c00",
  "Bacillus pumilus" = "#00FFF6",
  "Lactobacillus vaginalis" = "#cc3399",
  "Lactobacillus crispatus" = "#9900cc",
  "Bacillus subtilis" = "#666666",
  "Enterococcus hirae" = "#979300",
  "Enterococcus sp. CR-Ec1" = "#FF7700",
  "Eubacterium callanderi" = "#CC0000",
  "Enterococcus faecium" = "#FFAA00",
  "Paenibacillus silagei" = "#003300",
  "Lacrimispora saccharolytica" = "#6600cc",
  "Bacillus altitudinis" = "#991c00",
  "Fournierella massiliensis" = "#66cc00",
  "Lacrimispora amygdalina" = "#0000cc"
)


label_fill <- c(
  "AE" = "#c6dfb5", "AN" = "#b6d7e8",
  "BBE" = "#f4cccc", "CNA" = "#cfe2f3", "CNAB" = "#d9d2e9",
  "MRS" = "#ead1dc", "MSA" = "#fff2cc"
)

all_colors <- c(species_colors, label_fill)

species_order <- df_long %>%
  group_by(Species) %>%
  summarise(Total = sum(Abundance), .groups = "drop") %>%
  arrange(desc(Total)) %>%
  pull(Species)

df_long$Species <- factor(df_long$Species, levels = species_order)

ggplot(df_long,
       aes(axis1 = Aircondition, axis2 = Species, axis3 = Medium, y = Abundance)) +
  geom_alluvium(aes(fill = Species), width = 1/12, alpha = 0.85) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.3, color = "black") +
  scale_fill_manual(
    values = all_colors,
    breaks = species_order
  ) +
  scale_x_discrete(
    limits = c("Air condition", "Species", "Medium"),
    expand = c(0.05, 0.05)
  ) +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(
      size = 14,
      face = "bold",
      family = "Times New Roman",
      margin = margin(t = -5)
    ),
    axis.title.x = element_text(
      size = 14,
      face = "bold",
      family = "Times New Roman"
    ),
    axis.title.y = element_text(
      size = 14,
      face = "bold",
      family = "Times New Roman",
      margin = margin(r = 15)  
    ),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(
      size = 16,
      face = "bold.italic",
      family = "Times New Roman"
    ),
    legend.key.size = unit(1.2, "lines"), 
    plot.title = element_text(
      size = 16,
      face = "bold",
      hjust = 0.5,
      family = "Times New Roman"
    )
  ) +
  guides(fill = guide_legend(ncol = 1, keywidth = 1, keyheight = 1.2)) +
  labs(
    title = "Sankey Diagram: 20 Beneficial Species",
    x = NULL,
    y = "Relative Abundance (%)"
  )




### Fig.6A Bar plot showing the proportion of single-colony cultures classified as either single-species (>95% relative abundance) or multiple-species based on full-length 16S rRNA sequencing


library(ggplot2)
library(ggplot2)

df <- data.frame(
  Category = c("Single Species", "Multiple Species"),
  Count = c(150, 101)
)
df$Percent <- round(df$Count / sum(df$Count) * 100, 1)
df$Label <- paste0(df$Count, " (", df$Percent, "%)")

ggplot(df, aes(x = Category, y = Percent, fill = Category)) +
  geom_bar(stat = "identity", width = 0.3) +
  geom_text(
    aes(label = Label),
    vjust = -0.5,
    fontface = "bold",
    family = "Times New Roman",  
    size = 5
  ) +
  scale_fill_manual(values = c("Single Species" = "#66c2a5", "Multiple Species" = "#fc8d62")) +
  theme_minimal(base_family = "Times New Roman") + 
  labs(
    x = " ",
    y = "Number of Colonies",
    title = "Proportion of Single vs. Mixed-Species Colonies"
  ) +
  theme(
    axis.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    plot.title = element_text(hjust = 0, face = "bold", family = "Times New Roman"),
    plot.title.position = "plot",  
    legend.position = "none",
    plot.margin = margin(t = 25, r = 10, b = 10, l = 10)
  )





###Fig.6B Heatmap demonstrated the distribution of these 150 single-species isolates across 10 selected culture conditions.


library(tidyverse)
library(readr)
library(ComplexHeatmap)
library(circlize)
library(grid)

df <- read_csv("150 single species abundance table.csv")
df$Condition <- paste(df$Medium_name, df$Oxygen, sep = "_")

heatmap_df <- df %>%
  group_by(species, Condition) %>%
  summarise(Count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Condition, values_from = Count, values_fill = 0)

heatmap_matrix <- as.matrix(heatmap_df[, -1])
rownames(heatmap_matrix) <- heatmap_df$species


col_fun <- colorRamp2(c(0, 10, 30), c("white", "orange", "red"))

probiotic_species <- c(
  "Bacillus pumilus",
  "Pediococcus pentosaceus",
  "Ligilactobacillus salivarius",
  "Bacillus safensis",
  "Enterococcus faecalis",
  "Enterococcus hirae",
  "Bacteroides fragilis",
  "Lactobacillus crispatus",
  "Eubacterium callanderi"
)

row_label_colors <- ifelse(rownames(heatmap_matrix) %in% probiotic_species, "#377eb8", "black")


ht <- Heatmap(
  heatmap_matrix,
  name = "Number of Colonies",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  
  row_names_gp = gpar(
    fontsize = 12,
    fontface = "bold.italic",            
    fontfamily = "Times New Roman",
    col = row_label_colors
  ),
  column_names_gp = gpar(
    fontsize = 12,
    fontface = "bold",
    fontfamily = "Times New Roman"
  ),
  heatmap_legend_param = list(
    title_gp = gpar(
      fontsize = 12,
      fontface = "bold",
      fontfamily = "Times New Roman"
    ),
    labels_gp = gpar(
      fontsize = 12,
      fontface = "bold",
      fontfamily = "Times New Roman"
    )
  ),
  rect_gp = gpar(col = "grey80", lwd = 0.5),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(
      sprintf("%d", heatmap_matrix[i, j]),
      x, y,
      gp = gpar(
        fontsize = 12,
        fontface = "bold",
        fontfamily = "Times New Roman"
      )
    )
  }
)

grid.newpage()
draw(ht, heatmap_legend_side = "right", padding = unit(c(5, 3, 5, 5), "mm"))







###Fig S1A Venn diagram showed the numbers of unique and shared bacterial taxa identified from CD and CI approaches at the phylum level

library(eulerr)
library(grid)
data <- c("CD&CI" = 4, "CI" = 5, "CD" = 0)  
fit <- euler(data)
grid.newpage()
p <- plot(fit,
          fills = list(fill = c("#BDC9E1", "#A1D99B"), alpha = 0.6),
          edges = list(col = "black"),
          quantities = FALSE,
          labels = FALSE)
grid.draw(p)

grid.text("CD\n(Total: 4)", 
          x = 0.30, y = 0.50,
          gp = gpar(fontsize = 20, fontface = "bold"))

grid.text("CI\n(Total: 9)", 
          x = 0.75, y = 0.50,
          gp = gpar(fontsize = 20, fontface = "bold"))

grid.text("4", x = 0.45, y = 0.53,
          gp = gpar(fontsize = 20, fontface = "bold", family = "Times New Roman"))

grid.text("5", x = 0.65, y = 0.53,
          gp = gpar(fontsize = 20, fontface = "bold", family = "Times New Roman"))



###Fig S1B Venn diagram showed the numbers of unique and shared bacterial taxa identified from CD and CI approaches at the genus level

install.packages("eulerr")
install.packages("grid")
library(eulerr)
library(grid)

data <- c(
  CD = 43,
  CI = 152,
  "CD&CI" = 83
)

fit <- euler(data)

grid.newpage()
p <- plot(fit,
          fills = list(fill = c("#66C2A5", "#FC8D62"), alpha = 0.5),
          edges = list(col = "black"),
          quantities = FALSE, 
          labels = FALSE)

grid.draw(p)

grid.text("CD\n(Total: 126)", x = 0.10, y = 0.45, gp = gpar(fontsize = 18, fontface = "bold")) 
grid.text("CI\n(Total: 235)", x = 0.90, y = 0.45, gp = gpar(fontsize = 18, fontface = "bold")) 

grid.text("43", x = 0.22, y = 0.50,
          gp = gpar(fontsize = 18, fontface = "bold", family = "Times New Roman"))
grid.text("83", x = 0.40, y = 0.50,
          gp = gpar(fontsize = 18, fontface = "bold", family = "Times New Roman"))
grid.text("152", x = 0.68, y = 0.50,
          gp = gpar(fontsize = 18, fontface = "bold", family = "Times New Roman"))





###Fig S1C Stacked bar plot showed the relative abundance of bacterial taxa identified by the CD and CI approaches at the phylum level


library(ggplot2)
library(dplyr)
library(readr)
library(scales)  
library(grid)    

df <- read_delim("Phylum abundance in CD&CI.csv")

all_combinations <- expand.grid(
  Group = unique(df$Group),
  Phylum = c("Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", 
             "Tenericutes", "Synergistetes", "Fusobacteria", "Ignavibacteriae", "Aquificae")
)

df <- full_join(df, all_combinations, by = c("Group", "Phylum")) %>%
  mutate(
    Abundance = ifelse(is.na(Abundance), 0, Abundance),
    Abundance = ifelse(Abundance == 0, 0.0001, Abundance)
  )

df$Phylum <- factor(df$Phylum, levels = c(
  "Proteobacteria", "Firmicutes", "Bacteroidetes", "Actinobacteria", 
  "Tenericutes", "Synergistetes", "Fusobacteria", "Ignavibacteriae", "Aquificae"
))

phylum_colors <- c(
  "Proteobacteria" = "#377EB8",
  "Firmicutes" = "#4DAF4A",
  "Bacteroidetes" = "#E41A1C",
  "Actinobacteria" = "#999999",
  "Tenericutes" = "#984EA3",
  "Synergistetes" = "#FF7F00",
  "Fusobacteria" = "#A65628",
  "Ignavibacteriae" = "#F781BF",
  "Aquificae" = "#FFFF33"
)


p <- ggplot(df, aes(x = Group, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = c("0", "25", "50", "75", "100")
  ) +
  scale_fill_manual(values = phylum_colors, drop = FALSE) +
  labs(x = NULL, y = "Relative Abundance (%)", fill = "Phylum") +
  theme_minimal(base_family = "Times New Roman") +
  theme(
    axis.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title.y = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.text = element_text(size = 14, face = "bold.italic", family = "Times New Roman"),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))



print(p)






###Fig S1D Stacked bar plot showed the relative abundance of bacterial taxa identified by the CD and CI approaches at the genus level



install.packages("ggplot2")
install.packages("RColorBrewer")

library(ggplot2)
library(RColorBrewer)
library(readr)
library(dplyr)


df_genus <- read.csv("Genus abundance in CD&CI.csv")


genus_abundance <- aggregate(Abundance ~ Genus, data = df_genus, sum)


top_10_genus <- genus_abundance[order(-genus_abundance$Abundance), "Genus"][1:10]


highlight_genus <- c(
  "Escherichia" = "#377EB8",
  "Proteus" = "#4DAF4A",
  "Limosilactobacillus" = "#E41A1C",
  "Enterococcus" = "#FF7F00",
  "Ligilactobacillus" = "#984EA3",
  "Faecalibacterium" = "#A65628",
  "Eisenbergiella" = "#FFFF33",
  "Negativibacillus" = "#D95F02",
  "Lacrimispora" = "#7570B3",
  "Flintibacter" = "#54278F",
  "Others<0.5% in CI" = "#E7298A",
  "Blautia" = "#008B8B",
  "Bacillus" = "#8A2BE2"
)


other_genus <- setdiff(unique(df_genus$Genus), names(highlight_genus))
other_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(length(other_genus))


genus_colors <- c(highlight_genus, setNames(other_colors, other_genus))

p <- ggplot(df_genus, aes(x = Group, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = c("0", "25", "50", "75", "100")
  ) +
  scale_fill_manual(values = genus_colors, drop = FALSE) +
  theme_minimal(base_family = "Times New Roman") +  
  labs(x = "", y = "Relative Abundance (%)", fill = "Genus") +
  theme(
    axis.text = element_text(size = 14, face = "bold", family = "Times New Roman"),
    axis.title.y = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.title = element_text(size = 14, face = "bold", family = "Times New Roman"),
    legend.text = element_text(size = 12, face = "bold.italic", family = "Times New Roman"),  
    legend.key.size = unit(0.3, "cm"),
    legend.position = "right",
    legend.spacing.x = unit(0.2, "cm"),
    legend.spacing.y = unit(0.2, "cm")
  ) +
  guides(fill = guide_legend(ncol = 1))


print(p)
