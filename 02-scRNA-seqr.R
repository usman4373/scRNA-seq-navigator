# ==============================================================================
# Load Required Libraries
# ==============================================================================
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(openxlsx)
library(HGNChelper)
library(plotly)
#library(SeuratDisk) # Use when RAM is insufficient, it uses disk space to compensate

# ==============================================================================
# LOAD AND PREPARE DATA
# ==============================================================================

# Load merged data
rawdata <- merged_samples
dim(rawdata)  # 38224 features and 36289 cells
rm(merged_samples)

# calculate mitochondrial count
rawdata[["percent.mt"]] <- PercentageFeatureSet(rawdata, pattern = "^(?i)MT-")  # (?i) captures mt-Nd1, Mt-Nd1, and MT-Nd1

Idents(rawdata) <- "dataset"
#Idents(rawdata) <- data$orig.ident

# Visualize QC metrics as a violin plot
Plot1 <- VlnPlot(rawdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "condition", ncol = 3)
Plot1
ggsave(file = "raw_data.png", Plot1, dpi = 600, bg = "white")


# ==============================================================================
# DATA FILTRATION
# ==============================================================================
filtered_data <- subset(rawdata, subset = nFeature_RNA > 200
                        & nFeature_RNA < 6000 & nCount_RNA > 2000 & percent.mt < 10)
dim(filtered_data) # 38224 features and 16195 cells in filtered data

# visualize filtered data
filter_plt <- VlnPlot(filtered_data, features = c("nFeature_RNA", "nCount_RNA",
                                                  "percent.mt"),
                      group.by = "condition", ncol = 3)
filter_plt
ggsave(file = "filtered_data.png", filter_plt, dpi = 600, bg = "white")

# ==============================================================================
# DATA NORMALIZATION AND FEATURE SELECTION
# ==============================================================================
normalized_data <- NormalizeData(filtered_data, normalization.method = "LogNormalize", scale.factor = 10000)
normalized_data <- FindVariableFeatures(normalized_data, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(normalized_data), 10)
top10 # "S100A8" "SELENOP" "APOE" "APOC1" "C1QB" "G0S2" "CCL4L2" "IFI27" "IDO1" "SPP1" top 10 HVFs

# plot variable features with and without labels
plt1 <- VariableFeaturePlot(normalized_data)
plt2 <- LabelPoints(plot = plt1, points = top10, repel = TRUE)
variable_plt <- plt2
variable_plt
ggsave(file = "HVFs.png", variable_plt, dpi = 600, bg = "white")

# ==============================================================================
# DATA SCALING AND DIMENSIONALITY REDUCTION
# ==============================================================================

# ==============================================================================
# Scale Data and Perform PCA
# ==============================================================================
# Data Scaling
#all.genes <- rownames(normalized_data) # Run if want to use all genes of dataset
scaled_data <- ScaleData(normalized_data) # Add "features = all.genes" if ran above line

# Linear Dimensionality Reduction
scaled_data <- RunPCA(scaled_data, features = VariableFeatures(object = scaled_data))

# Examine and visualize PCA results a few different ways
print(scaled_data[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scaled_data, dims = 1:2, reduction = "pca")
pca_plot <- DimPlot(scaled_data, reduction = "pca")
pca_plot
ggsave(file = "PCA_plot.png", pca_plot, dpi = 600, bg = "white")

# ==============================================================================
# Determine Data Dimensions
# ==============================================================================
elbow <- ElbowPlot(scaled_data, ndim = 50)
elbow
ggsave(file = "Elbow_plot.png", elbow, dpi = 600, bg = "white")

# ==============================================================================
# CLUSTERING AND NON-LINEAR DIMENSIONAL REDUCTION
# ==============================================================================

# ==============================================================================
# Cluster Data and Run UMAP/t-SNE
# ==============================================================================
# Data Clustering
data <- FindNeighbors(scaled_data, dims = 1:20) #  graph.name = c("nn_graph", "snn_graph")
data <- FindClusters(data, resolution = 0.6, algorithm = 4) # Add "graph.name" argument if used in above line of code

levels(data) # Gives number of clusters generated

# Non-Linear Dimensional Reduction
data <- RunUMAP(data, dims = 1:20)
data <- RunTSNE(data, dims = 1:20)

# ==============================================================================
# Visualize Clustering Results
# ==============================================================================
# umap and tsne plots
umap <- DimPlot(data, reduction = "umap", label = TRUE)
umap
ggsave(file = "umap.png", umap, dpi = 600, bg = "white")

tSNE <- DimPlot(data, reduction = "tsne", label = TRUE)
tSNE
ggsave(file = "tsne.png", tSNE, dpi = 600, bg = "white")

# umap and tsne plots with conditions
umap_condition <- DimPlot(data, reduction = "umap", label = TRUE,
                          split.by = "condition")
umap_condition
ggsave(file = "umap_condition.png", umap_condition, dpi = 600, bg = "white")

tSNE_condition <- DimPlot(data, reduction = "tsne", label = TRUE,
                          split.by = "condition")
tSNE_condition
ggsave(file = "tsne_condition.png", tSNE_condition, dpi = 600, bg = "white")

# ==============================================================================
# CELL TYPE ANNOTATION USING SCTYPE
# ==============================================================================

# ==============================================================================
# Load scType Functions and Prepare Gene Sets
# ==============================================================================

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# DB file
db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue <- "Immune system" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, tissue)

# ==============================================================================
# Run scType Annotation
# ==============================================================================
# check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
seurat_package_v5 <- isFALSE('counts' %in% names(attributes(data[["RNA"]])));
print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))

# extract scaled scRNA-seq matrix
scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(data[["RNA"]]$scale.data) else as.matrix(data[["RNA"]]@scale.data)

# run ScType
es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. For raw (unscaled) count matrix set scaled = FALSE
# When using Seurat, we use "RNA" slot with 'scale.data' by default. Please change "RNA" to "SCT" for sctransform-normalized data,
# or to "integrated" for joint dataset analysis. To apply sctype with unscaled data, use e.g. pbmc[["RNA"]]$counts or pbmc[["RNA"]]@counts, with scaled set to FALSE.

# ==============================================================================
# Process scType Results and Assign Cell Types
# ==============================================================================
# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] <- "Unknown"
print(sctype_scores[,1:3])
write.csv(sctype_scores, file = "Celltypes.csv")

# Convert sctype_scores to a data frame if it isn't already
sctype_scores <- as.data.frame(sctype_scores)

# Extract the cluster and cell type columns
cluster_ids <- sctype_scores$cluster
cell_types <- sctype_scores$type

# Create a named vector for the new cluster identities
new_cluster_ids <- cell_types
names(new_cluster_ids) <- cluster_ids

# Rename the cluster identities in the object
data <- RenameIdents(data, new_cluster_ids)
levels(data)

# ==============================================================================
# Verify the Renaming of Cluster Ids to Annotated Cell Types
# ==============================================================================
renamed_clusters <- data.frame(
ClusterID = names(new_cluster_ids),
CellType = as.character(new_cluster_ids))

# Print the renamed clusters to verify
print(renamed_clusters) # Open Celltypes.csv file and match the cluster ids with the cell types

# ==============================================================================
# Add Cell Type Classification to Metadata
# ==============================================================================
# Overlay the identified cell types on UMAP plot
data@meta.data$sctype_classification = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  data@meta.data$sctype_classification[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

# Verify visually by comparing unannotated and annotated tsne plots and matching with Celltypes.csv file
tsne_annt <- DimPlot(data, reduction = "tsne", label = TRUE,
                     repel = TRUE, group.by = 'sctype_classification')

verf_annt <- tSNE + tsne_annt
verf_annt
ggsave(file = "annotation_verification.png", verf_annt, dpi = 600, bg = "white")


# ==============================================================================
# Count Cells by Cell Type and Condition
# ==============================================================================
# Number of cells in each condition
# Create the table
num_of_cells <- table(data@meta.data$sctype_classification, data@meta.data$condition)
# Convert to data frame
num_of_cells_df <- as.data.frame.matrix(num_of_cells)
num_of_cells_df
# Save as CSV
write.csv(num_of_cells_df, file = "number_of_cells_in_conditions.csv", row.names = TRUE)

# ==============================================================================
# VISUALIZE ANNOTATED DATA
# ==============================================================================

# umap and tsne plots
umap_annt <- DimPlot(data, reduction = "umap", label = FALSE,
                     repel = TRUE, group.by = 'sctype_classification')        
umap_annt
ggsave(file = "umap_annotated.png", umap_annt, dpi = 600, bg = "white")

tsne_annt <- DimPlot(data, reduction = "tsne", label = FALSE,
                     repel = TRUE, group.by = 'sctype_classification')        
tsne_annt
ggsave(file = "tsne_annotated.png", tsne_annt, dpi = 600, bg = "white")

umap_annt_condition <- DimPlot(data, reduction = "umap", label = FALSE,
                               repel = TRUE, group.by = 'sctype_classification',
                               split.by = "condition")        
umap_annt_condition
ggsave(file = "umap_annotated_condition.png", umap_annt_condition, dpi = 600, bg = "white")

tsne_annt_condition <- DimPlot(data, reduction = "tsne", label = FALSE,
                               repel = TRUE, group.by = 'sctype_classification',
                               split.by = "condition")        
tsne_annt_condition
ggsave(file = "tsne_annotated_condition.png", tsne_annt_condition, dpi = 600, bg = "white")

saveRDS(data, file = "Annotated_data.rds")

# ==============================================================================
# DATA SUBSETTING (Optional)
# ==============================================================================

# ==============================================================================
# Subset Specific Cell Types If Needed
# ==============================================================================
# Check available cell types in the metadata
num_of_cells_df
# cells to keep
cells <- c("Memory CD4+ T cells",
           "Naive CD4+ T cells",
           "Cancer cells",
           "CD8+ NKT-like cells",
           "Natural killer  cells")
# Subset based on metadata column
data_subset <- subset(data, subset = sctype_classification %in% cells)
# Verify the subset
levels(data_subset)
table(data_subset@meta.data$sctype_classification)
table(data_subset@meta.data$condition)
saveRDS(data_subset, file = "data_subset.rds")

# ==============================================================================
# DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ==============================================================================

# Join Layers before DGE if not already
data <- JoinLayers(data) # Adjust object name accordingly

# ==============================================================================
# Differential Expression Between Cell Types in Different Conditions
# ==============================================================================
data$celltype.condition <- paste(data$sctype_classification, data$condition, sep = "_")
Idents(data) <- "celltype.condition"
levels(data)

dff_markers <- FindMarkers(data, ident.1 = "Memory CD4+ T cells_LM",
                           ident.2 = "Memory CD4+ T cells_Normal", min.pct = 0.1,
                           test.use = "wilcox", verbose = TRUE, only.pos = FALSE)
# Apply FDR
dff_markers$fdr = p.adjust(dff_markers$p_val, method='fdr')
write.csv(dff_markers, file = "Memory CD4+ T cells_DEGs.csv") # Required for volcano plot

dff_up_markers <- subset(dff_markers, avg_log2FC > 0.25)
dff_dn_markers <- subset(dff_markers, avg_log2FC < -0.25)
# Save degs results
write.csv(dff_up_markers, file = "Memory CD4+ T cells_up_markers.csv")
write.csv(dff_dn_markers, file = "Memory CD4+ T cells_dn_markers.csv")

# ==============================================================================
# Differential Expression Between Different Cell Types Regardless of Condition
# ==============================================================================
Idents(data) <- "sctype_classification"
levels(data)
degs <- FindMarkers(data, ident.1 = "Cancer cells",
                    ident.2 = NULL, test.use = "wilcox",
                    only.pos = FALSE)
# Apply FDR
degs$fdr = p.adjust(degs$p_val, method='fdr')
write.csv(degs, file = "Cancer cells_DEGs.csv") # Required for volcano plot

degs_up <- subset(degs, avg_log2FC > 0.25)
degs_down <- subset(degs, avg_log2FC < -0.25)
# Save degs results
write.csv(degs_up, file = "Cancer_cells_up_markers.csv")
write.csv(degs_down, file = "Cancer_cells_dn_markers.csv")

# ==============================================================================
# Differential Expression Between Different Cell Types in Different Conditions
# ==============================================================================
data$celltype.condition <- paste(data$sctype_classification, data$condition, sep = "_")
Idents(data) <- "celltype.condition"
levels(data)
# DEGs of Cancer cells in diseased condition vs All other cell types in Normal condition
degs <- FindMarkers(data, ident.1 = "Cancer cells_LM",
                    ident.2 = c("Macrophages_Normal",
                                "Natural killer  cells_Normal",
                                "Basophils_Normal",
                                "Non-classical monocytes_Normal",
                                "CD8+ NKT-like cells_Normal",
                                "Memory CD4+ T cells_Normal",
                                "Myeloid Dendritic cells_Normal",
                                "Pre-B cells_Normal",
                                "Endothelial_Normal",
                                "Classical Monocytes_Normal",
                                "Plasmacytoid Dendritic cells_Normal",
                                "γδ-T cells_Normal",
                                "Memory B cells_Normal",
                                "Naive CD4+ T cells_Normal"
                    ), test.use = "wilcox",
                    only.pos = FALSE)
# Apply FDR
degs$fdr = p.adjust(degs$p_val, method='fdr')
write.csv(degs, file = "Cancer cells_DEGs.csv") # Required for volcano plot

degs_up <- subset(degs, avg_log2FC > 0.25)
degs_down <- subset(degs, avg_log2FC < -0.25)
# Save degs results
write.csv(degs_up, file = "Cancer_cells_up_markers.csv")
write.csv(degs_down, file = "Cancer_cells_dn_markers.csv")

# ==============================================================================
# END OF SCRIPT
# ==============================================================================