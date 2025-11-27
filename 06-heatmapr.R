# ==============================================================================
# Load Required Libraries
# ==============================================================================
library(Seurat)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(paletteer)

# ==============================================================================
# Load Data and Extract Gene Expression
# ==============================================================================
data <- readRDS("Annotated_data.rds")

# If you store expression in layers, ensure you joined them
data <- JoinLayers(data)

genes_of_interest <- c(
  "CD247","RAC1","CALM1","CD3D","HLA-DQA1",
  "HLA-DQA2","HLA-DPA1","HLA-DRA","CD74","CD3G"
)

# Extract expression (genes x cells)
expr_mat <- GetAssayData(data, assay = "RNA", layer = "data")[genes_of_interest, ]

# Transpose to cells x genes and make a data.frame
expr_df <- as.data.frame(t(expr_mat))

# Attach metadata (make sure these column names exist in your object)
metadata <- data@meta.data
expr_df$celltype <- metadata$sctype_classification
expr_df$condition <- metadata$condition

# ==============================================================================
# STEP 1 — Z-SCORE each gene across ALL CELLS (cell-level z-scores)
# ==============================================================================
# This centers & scales each gene using all cells: mean = 0, sd = 1
expr_df_zcell <- expr_df %>%
  mutate(across(all_of(genes_of_interest), ~ scale(.)[, 1]))

# ==============================================================================
# STEP 2 — AVERAGE z-scores WITHIN each celltype × condition group
# (So each group will become a single column in the heatmap)
# ==============================================================================
z_group_df <- expr_df_zcell %>%
  group_by(celltype, condition) %>%
  summarise(across(all_of(genes_of_interest), mean, .names = "{.col}"), .groups = "drop")

# Create a column name that combines celltype and condition (for rownames later)
z_group_df <- z_group_df %>%
  mutate(group = paste(celltype, condition, sep = "_"))

# Save the group-level z-scores table
write.csv(z_group_df, "Genes-z-scores.csv", row.names = FALSE)

# ==============================================================================
# PREPARE MATRIX FOR HEATMAP
# ==============================================================================
# Build matrix genes x groups (ComplexHeatmap expects rows = genes, cols = samples/groups)
expr_mat_group <- as.matrix(z_group_df[, genes_of_interest])  # rows = groups, cols = genes
rownames(expr_mat_group) <- z_group_df$group

# transpose so rows = genes and columns = groups (as in your original)
expr_mat_group <- t(expr_mat_group)  # now rows = genes, cols = celltype_condition

# ==============================================================================
# MANUALLY ORDER COLUMNS BY CONDITION THEN CELLTYPE
# ==============================================================================
# Get unique conditions and celltypes
conditions <- unique(z_group_df$condition)
celltypes <- unique(z_group_df$celltype)

# Create the desired column order: all CRC groups first, then all Control groups
desired_order <- c()
for(cond in conditions) {  # This will loop through conditions in their natural order
  for(ct in celltypes) {
    group_name <- paste(ct, cond, sep = "_")
    if(group_name %in% colnames(expr_mat_group)) {
      desired_order <- c(desired_order, group_name)
    }
  }
}

# Reorder the matrix columns
expr_mat_group <- expr_mat_group[, desired_order, drop = FALSE]

# ==============================================================================
# Prepare column annotation (for groups)
# ==============================================================================
# Extract condition and celltype per group (match order from reordered matrix)
annotation_col <- data.frame(
  condition = factor(sub(".*_", "", colnames(expr_mat_group))),
  celltype = factor(sub("_.*", "", colnames(expr_mat_group)))
)
rownames(annotation_col) <- colnames(expr_mat_group)

# ==============================================================================
# Annotation colors (customize as needed)
# ==============================================================================
# Define your full color vector
all_colors <- c(
  "#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b",
  "#abdda4","#66c2a5","#3288bd","#5e4fa2","#6090D8",
  "#4C825D","#8CAE9E","#D5AB85","#99610A","#C38F16",
  "#FFCD00","#C0532B","#B50200","#67322E","#B47E83",
  "#807B7F","#8C86A0","#674D53","#508CA7","#0E2A4D"
)

# Extract unique cell types in the order they appear in annotation_col
celltypes <- unique(annotation_col$celltype)

# Assign colors to cell types (truncate color vector if needed)
celltype_colors <- setNames(all_colors[1:length(celltypes)], celltypes)

# Build annotation color list
ann_colors <- list(
  condition = c(LM = "#2c98a0", Normal = "#67dba5"),  # adjust for more conditions if needed
  celltype = celltype_colors
)


# ==============================================================================
# Heatmap colors & plotting
# ==============================================================================
# Use a diverging color scale centered at 0 (z-scores)
hcolors <- rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging", n = 200))

ht <- Heatmap(
  expr_mat_group,
  name = "Z-score",
  col = hcolors,
  cluster_rows = TRUE,
  cluster_columns = TRUE,  # Turn off column clustering to use manual order, if not then turn on "TRUE"
  show_row_names = TRUE,
  show_column_names = FALSE,
  top_annotation = HeatmapAnnotation(df = annotation_col, col = ann_colors,
                                     # Add borders to annotation cells
                                     gp = gpar(col = "black", lwd = 0.8)),
  # Add borders to heatmap cells
  rect_gp = gpar(col = "black", lwd = 0.8),
  heatmap_legend_param = list(
    legend_width = unit(0.8, "cm"),
    direction = "vertical",
    title_position = "topleft"
  ),
  use_raster = FALSE
)
ht

# Save as PNG (high res)
png("Heatmap.png", width = 3000, height = 2000, res = 300)
draw(ht)
dev.off()
