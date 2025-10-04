#Isolates recovery ratio
#maiteortuzar@usal.es

#Libraries
library(tidyverse)
library(readr)
library(randomcoloR)

#Read tables
metagen <- read_tsv("metagenomica_tr_2.tsv")
meta <- read_tsv("metadatos_tr_2.tsv")
aislados <- read_tsv("aislados_tr.tsv")

#Reshape metagenomics table into long format
metagen_long <- metagen %>%
  pivot_longer(cols = -(1:6), names_to = "ID", values_to = "rel_abund") %>%
  left_join(meta, by = "ID") %>%
  select(Genus, ID, rel_abund, Location, Compartment)

#Summarize isolates by genus
aislados_summary <- aislados %>%
  group_by(Genus = genus, Location, Compartment) %>%
  summarise(observed_isolates = n(), .groups = "drop")

#Calculate total number of isolates per sample
total_isolates_per_sample <- aislados %>%
  group_by(Location, Compartment) %>%
  summarise(N = n(), .groups = "drop")

#Combine relative abundance with isolate counts
df <- metagen_long %>%
  left_join(total_isolates_per_sample, by = c("Location","Compartment")) %>%
  left_join(aislados_summary, by = c("Genus","Location","Compartment")) %>%
  mutate(observed_isolates = replace_na(observed_isolates, 0),
         expected_isolates = rel_abund/100 * N,
         recovery_ratio = ifelse(expected_isolates > 0, observed_isolates / expected_isolates, NA))

#Scatter plot
df_scatter <- df %>%
  filter(!is.na(recovery_ratio) & recovery_ratio > 0)
genus_list <- unique(df_scatter$Genus)
n_genus <- length(genus_list)
pal <- distinctColorPalette(n_genus)
names(pal) <- genus_list
ggplot(df_scatter, aes(x = rel_abund, y = observed_isolates, color = Genus)) +
  geom_point(alpha = 0.7, size = 3) +
  scale_x_log10(labels = label_number()) +
  scale_y_log10(labels = label_number()) +
  scale_color_manual(values = pal) +
  labs(x = "Relative abundance (%)",
       y = "Number of isolates",
       color = "Genus",
       title = "Recovery rates") +
  facet_wrap(~ Location) +
  theme_bw() +
  guides(color = guide_legend(ncol = 2)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    strip.background = element_rect(fill = "grey90", color = "black")
  )

#Export results
write_tsv(df, "df_scatter_export.tsv")
df_scatter <- df %>%
  filter(!is.na(recovery_ratio) & recovery_ratio > 0)
df_scatter_export <- df_scatter %>%
  dplyr::select(rel_abund, observed_isolates, Genus, Location)
write.csv(df_scatter_export, "Scatterplot_data.csv", row.names = FALSE)
