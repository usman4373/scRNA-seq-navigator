# ==============================================================================
# Load Required Libraries
# ==============================================================================
library(paletteer)
library(ggplot2)
library(Seurat)

# ==============================================================================
# BOX PLOT - CONDITION-BASED EXPRESSION
# ==============================================================================

# Load the data if not already
data <- readRDS("Annotated_data.rds")

# Create a data frame with gene expression and conditions instead of clusters
boxplot_df <- data.frame(gene = GetAssayData(data, assay = "RNA", layer = "data")["HSP90AA1",], 
                         condition = data$condition)

custom_colors <- c("LM" = "#FF3200FF", "Normal" = "#0076BBFF") #, "Metastasis" = "#2A9D8F")

box <- ggplot(boxplot_df, aes(x = condition, y = gene)) + 
  geom_jitter(aes(color = condition), width = 0.2, size = 1.5, alpha = 0.6, show.legend = FALSE) +  
  geom_boxplot(aes(fill = condition), width = 0.7, outlier.shape = NA, alpha = 0.6, color = "black") +  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "HSP90AA1", x = "Condition", y = "Expression") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.x.top    = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.ticks.x.top   = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.ticks.x.bottom   = element_line(),
    axis.ticks.y.left = element_line(),
    panel.grid.major = element_line(color = "grey80", size = 0.4)
  )
box
# Save the plot
ggsave("HSP90AA1-condition.png", plot = box, bg = "white", dpi = 600)

# ==============================================================================
# BOX PLOT - CELL TYPE SPECIFIC EXPRESSION
# ==============================================================================

# Specify the cell type of interest
target_celltype <- "Memory CD4+ T cells"

# Logical filter: keep rows that start with the cell type (e.g., "Cancer cells")
cell_indices <- startsWith(data$celltype.condition, paste0(target_celltype, "_"))

# Create filtered data frame for box plot
boxplot_df <- data.frame(
  gene = GetAssayData(data, assay = "RNA", layer = "data")["HSP90AA1", cell_indices], 
  condition = data$celltype.condition[cell_indices]
)

# Recalculate group levels for coloring
groups <- unique(boxplot_df$condition)

custom_colors <- setNames(
  ifelse(endsWith(groups, "_LM"), "#FF3200FF", "#0076BBFF"),
  groups
)

# Plot
box <- ggplot(boxplot_df, aes(x = condition, y = gene)) + 
  geom_jitter(aes(color = condition), width = 0.2, size = 1.5, alpha = 0.6, show.legend = FALSE) +  
  geom_boxplot(aes(fill = condition), width = 0.7, outlier.shape = NA, alpha = 0.6, color = "black") +  
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = paste("HSP90AA1"), x = "Condition", y = "Expression") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.line.x.top    = element_blank(),
    axis.line.y.right  = element_blank(),
    axis.ticks.x.top   = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.ticks.x.bottom   = element_line(),
    axis.ticks.y.left = element_line(),
    panel.grid.major = element_line(color = "grey80", size = 0.4)
  )
box
# Save the plot
ggsave("HSP90AA1-Memory-CD4-T-cells.png", plot = box, bg = "white", dpi = 600)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================