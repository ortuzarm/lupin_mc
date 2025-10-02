#Microbiome data processing and visualization script (diversity and community)
#maiteortuzar@usal.es


#Libraries
library(tidyverse)
library(vegan)
library(ape)
library(rstatix)
library(RColorBrewer)
library(multcompView)
library(colorspace)
library(readr)

#Data
counts <- read_tsv("NC.tsv", show_col_types = FALSE)
taxonomy <- read_tsv("taxonomy.tsv", show_col_types = FALSE) %>%
  rename(FeatureID = `Feature_ID`, Taxon = Taxon, Confidence = Confidence)
meta <- read_tsv("metadatos.tsv", show_col_types = FALSE)

#Merge counts with taxonomy
if(!"ASV" %in% names(counts)) stop("Counts file must contain 'ASV' column.")

counts_tax <- counts %>%
  left_join(taxonomy %>% select(FeatureID, Taxon), by = c("ASV" = "FeatureID")) %>%
  mutate(Taxon = if_else(is.na(Taxon) | Taxon=="" , "Unassigned", Taxon))

sample_cols <- setdiff(names(counts_tax), c("ASV","Taxon"))
counts_tax[sample_cols] <- lapply(counts_tax[sample_cols], as.numeric)

counts_by_taxon <- counts_tax %>%
  group_by(Taxon) %>%
  summarise(across(all_of(sample_cols), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup()

write_tsv(counts_by_taxon, "counts_by_taxonomy_NC.tsv")

#Relative Abundance
counts_by_taxon_rel <- counts_by_taxon %>%
  rowwise() %>%
  mutate(across(all_of(sample_cols), as.numeric)) %>%
  ungroup() %>%
  mutate(across(all_of(sample_cols), ~ .x / sum(.x), .names = "{.col}"))

write_tsv(counts_by_taxon_rel, "NC_rel_abundance_by_taxonomy.tsv")

otu_tax <- counts_by_taxon %>%
  column_to_rownames("Taxon") %>% t() %>% as.data.frame()

common_samples <- intersect(rownames(otu_tax), meta$ID)
meta_sub <- meta %>% filter(ID %in% common_samples) %>% arrange(match(ID, common_samples))
otu_tax <- otu_tax[meta_sub$ID, , drop = FALSE]

#Bray-Curtis distance and PCoA
bc_dist <- vegdist(otu_tax, method = "bray")
pcoa_res <- pcoa(bc_dist)
pcoa_axes <- as.data.frame(pcoa_res$vectors[,1:6])
var_exp <- 100 * (pcoa_res$values$Eigenvalues / sum(pcoa_res$values$Eigenvalues))
pct1 <- round(var_exp[1], 1)
pct2 <- round(var_exp[2], 1)

PCOA_NC <- meta_sub %>%
  bind_cols(pcoa_axes %>% select(1,2) %>% rename(PC1 = 1, PC2 = 2)) %>%
  mutate(Compartment = factor(Compartment, levels = c("Soil", "Rhizosphere","Roots")),
         types = paste(Treatment, Compartment))

#PERMANOVA
adonis_res <- adonis2(bc_dist ~ Treatment * Compartment, data = meta_sub, permutations = 999, by = "margin")
permanova_label <- paste0("PERMANOVA:\nTreatment R2=", round(adonis_res["Treatment","R2"],3)," p=", signif(adonis_res["Treatment","Pr(>F)"],3),
                          "\nCompartment R2=", round(adonis_res["Compartment","R2"],3)," p=", signif(adonis_res["Compartment","Pr(>F)"],3))

#PCoA plot
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual', ]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
n <- 30
if(length(unique(PCOA_NC$Treatment)) > n) stop("More than 30 treatments; increase n.")

treat_levels <- unique(PCOA_NC$Treatment)
colors_assigned <- setNames(col_vector[1:length(treat_levels)], treat_levels)

PLOT_NC <- ggplot(PCOA_NC, aes(x = PC1, y = PC2, color = Treatment, shape = Compartment)) +
  geom_point(size = 4) +
  labs(x = paste0("PC1 (", pct1, "%)"), y = paste0("PC2 (", pct2, "%)"), title = "PCoA - Nc") +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold.italic")) +
  scale_color_manual(values = colors_assigned) +
  scale_shape_manual(values = c(16, 3, 17)) +
  annotate("text", x = min(PCOA_NC$PC1, na.rm=TRUE), y = max(PCOA_NC$PC2, na.rm=TRUE),
           label = permanova_label, hjust = 0, vjust = 1, size = 3.5)

print(PLOT_NC)
ggsave("PCoA_NC.png", PLOT_NC, width = 8, height = 6, dpi = 300)
write_tsv(PCOA_NC, "PCoA_NC_data_for_plot.tsv")

#Alpha diversity
samples_all <- rownames(otu_tax)
otu_all <- otu_tax[samples_all, , drop = FALSE]
meta_all_sub <- meta %>% filter(ID %in% samples_all) %>% arrange(match(ID, samples_all))

alpha_df <- tibble(
  ID = rownames(otu_all),
  Observed = rowSums(otu_all > 0),
  Shannon = diversity(otu_all, index = "shannon"),
  Simpson = diversity(otu_all, index = "simpson")
) %>% left_join(meta_all_sub, by = "ID")

write_tsv(alpha_df, "alpha_diversity_all_samples_per_sample.tsv")

#Boxplot Shannon by compartment
p_shannon_comp <- ggplot(alpha_df, aes(x = Compartment, y = Shannon, fill = Compartment)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, alpha = 0.6) +
  labs(title = "Shannon (All Samples) by Compartment") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

print(p_shannon_comp)
write_tsv(alpha_df, "alpha_diversity_all_samples_data_for_plot.tsv")

#Boxplots by treatment within each compartment
for(metric in c("Observed","Shannon","Simpson")){
  for(comp in unique(alpha_df$Compartment)){
    dfc <- alpha_df %>% filter(Compartment == comp)
    if(nrow(dfc) == 0) next
    p <- ggplot(dfc, aes(x = Treatment, y = .data[[metric]], fill = Treatment)) +
      geom_boxplot() + geom_jitter(width = 0.2, alpha = 0.6) +
      labs(title = paste0(metric, " (All Samples) - ", comp), x = "Treatment", y = metric) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
    write_tsv(dfc, paste0("NC_alpha_", metric, "_all_samples_byTreatment_data_", comp, ".tsv"))
  }
}

#Relative abundance bar plots
data <- read_tsv("nc_barras.tsv", show_col_types = FALSE)
data_long <- data %>%
  pivot_longer(cols = Methylophilaceae:Others, names_to = "Genus", values_to = "Abundance")
top_genera <- data_long %>%
  group_by(Genus) %>%
  summarise(Total_Abundance = sum(Abundance, na.rm = TRUE)) %>%
  arrange(desc(Total_Abundance)) %>%
  filter(Genus != "Others") %>%
  slice_head(n = 29) %>%
  bind_rows(tibble(Genus = "Others", Total_Abundance = 0))
data_filtered <- data_long %>%
  filter(Genus %in% top_genera$Genus) %>%
  mutate(Genus = factor(Genus, levels = top_genera$Genus))
compartment_order <- c("Soil", "Rhizosphere", "Roots")
data_filtered$Compartment <- factor(data_filtered$Compartment, levels = compartment_order)
data_normalized <- data_filtered %>%
  group_by(Syncom, Compartment) %>%
  mutate(Abundance_rel = Abundance / sum(Abundance) * 100) %>%
  ungroup()

n_colors <- length(levels(data_normalized$Genus))
paleta_base <- colorRampPalette(brewer.pal(8, "Pastel2"))(n_colors)
paleta_intensa <- darken(paleta_base, amount = 0.3)
genus_colors <- c(setNames(paleta_intensa, levels(data_normalized$Genus)), "Others" = "#A0A0A0")

ggplot(data_normalized, aes(x = Syncom, y = Abundance_rel, fill = Genus)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Compartment, scales = "free_x") +
  scale_fill_manual(values = genus_colors) +
  labs(title = "Relative Abundance by Syncom and Compartment", x = "Syncom", y = "Relative Abundance (%)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(size = 12, face = "bold"))

write_tsv(data_normalized, "barras_normalized_data_for_plot.tsv")
