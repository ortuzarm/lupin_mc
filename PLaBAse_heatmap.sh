#PLaBAse SynCom strain genome annotations (heatmap)
#Directory containing the TSV files
directory <- "C:/Users/maiteortuzar/Desktop/"

#Level to merge on (can be changed to other levels)
merge_level <- "level3"

#Libraries
library(ComplexHeatmap)
library(dplyr)
library(RColorBrewer)
library(reshape2)
library(circlize)

#Color vector
qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# Load all .txt files
files <- list.files(directory, pattern = "\\.txt$", full.names = TRUE)

dfs <- list() # Create empty list
for (file in files) {
  sample_name <- tools::file_path_sans_ext(basename(file))  # Extract sample name from file
  df <- read.table(file, sep = "\t", header = TRUE)
  colnames(df)[colnames(df) == "freq"] <- sample_name  # Rename freq column to sample name
  dfs[[length(dfs) + 1]] <- df
}

#Combine all dataframes
combined_df <- Reduce(function(x, y) merge(x, y, by = c("level1", "level2", "level3", "level4", "level5", "level6"), all = TRUE), dfs)
combined_df[is.na(combined_df)] <- 0 # Fill NA cells with 0

write.table(combined_df, file = file.path(directory, "combined_data.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
combined_data <- combined_df

merge_level_pos <- which(colnames(combined_data) == merge_level)
if (length(merge_level_pos) == 0) {
  stop(paste("Merge level", merge_level, "not found in the data."))
}

id_vars <- colnames(combined_data)[1:merge_level_pos]

numeric_vars <- colnames(combined_data)[sapply(combined_data, is.numeric)]
summed_data <- combined_data %>%
  group_by(across(all_of(id_vars))) %>%
  summarise(across(all_of(numeric_vars), sum, .names = "{.col}"), .groups = 'drop')

#Order data
summed_data <- summed_data %>%
  arrange(level1, level2, level3)

#Factorize categories
summed_data$level2 <- factor(summed_data$level2, levels = unique(summed_data$level2))
summed_data$level3 <- factor(summed_data$level3, levels = unique(summed_data$level3))

#Create heatmap matrix
heatmap_data <- summed_data %>% select(all_of(numeric_vars))
heatmap_matrix <- as.matrix(heatmap_data)
rownames(heatmap_matrix) <- summed_data[[merge_level]]
heatmap_matrix[heatmap_matrix > 100] <- 100 # Cap values at 100

# Assign colors dynamically for annotations
level3_colors <- setNames(col_vector[1:length(unique(summed_data$level3))], levels(summed_data$level3))
level2_colors <- setNames(col_vector[(length(level3_colors)+1):(length(level3_colors) + length(unique(summed_data$level2)))], levels(summed_data$level2))
col_fun <- colorRamp2(c(min(heatmap_matrix), 1,100), c("#FFDAB9", "#B0E2FF", "#00688B"))

#Function to format labels: remove underscores and lowercase
format_labels <- function(labels) {
  labels <- gsub("_", " ", labels)
  labels <- tolower(labels)
  return(labels)
}

#Condensed text annotations for level3
level3_labels <- format_labels(levels(summed_data$level3))
level3_text <- rep("", nrow(summed_data))
level3_mid_indices <- sapply(level3_labels, function(label) {
  indices <- which(format_labels(summed_data$level3) == label)
  return(indices[ceiling(length(indices) / 2)])
})
level3_mid_indices <- unlist(level3_mid_indices)
level3_text[level3_mid_indices] <- level3_labels

row_anno_level3 <- rowAnnotation(
  level3 = factor(summed_data$level3, levels = levels(summed_data$level3)),
  level3_text = anno_text(level3_text, gp = gpar(fontsize = 6)),
  col = list(level3 = level3_colors),
  show_legend = FALSE
)

#Condensed text annotations for level2
level2_labels <- format_labels(levels(summed_data$level2))
level2_text <- rep("", nrow(summed_data))
level2_mid_indices <- sapply(level2_labels, function(label) {
  indices <- which(format_labels(summed_data$level2) == label)
  return(indices[ceiling(length(indices) / 2)])
})
level2_mid_indices <- unlist(level2_mid_indices)
level2_text[level2_mid_indices] <- level2_labels

row_anno_level2 <- rowAnnotation(
  level2 = factor(summed_data$level2, levels = levels(summed_data$level2)),
  level2_text = anno_text(level2_text, gp = gpar(fontsize = 8)),
  col = list(level2 = level2_colors),
  show_legend = FALSE
)

#Create heatmaps split by level3 and level2
ht_level3 <- Heatmap(
  heatmap_matrix,
  name = "genes found",
  row_title = "sample",
  show_row_dend = FALSE,
  cluster_rows = FALSE,
  col = col_fun,
  show_row_names = FALSE,
  heatmap_legend_param = list(at = c(0, 50, 100), labels = c("0", "50", ">100")),
  border = TRUE,
  split = summed_data$level3,
  column_names_gp = grid::gpar(fontsize = 6)
)

ht_level2 <- Heatmap(
  heatmap_matrix,
  name = "genes found",
  row_title = "sample",
  show_row_dend = FALSE,
  cluster_rows = FALSE,
  col = col_fun,
  show_row_names = FALSE,
  heatmap_legend_param = list(at = c(0, 50, 100), labels = c("0", "50", ">100")),
  border = TRUE,
  split = summed_data$level2,
  column_names_gp = grid::gpar(fontsize = 6)
)

#Draw heatmaps with annotations
ht1 <- ht_level3 + row_anno_level3
ht2 <- ht_level2 + row_anno_level2

draw(ht1, heatmap_legend_side = "right")
draw(ht2, heatmap_legend_side = "right")
