# ==============================================================================
# Load Required Libraries
# ==============================================================================
library(Seurat)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
library(openxlsx)
library(HGNChelper)
library(plotly)
#library(SeuratDisk) # Use when RAM is insufficient, it uses disk space to compensate

# ==============================================================================
# PROCESS ALL SAMPLES FROM ALL DATASETS AND PREPARE DATA
# ==============================================================================

# Dataset 01 Samples
s1 <- Read10X("samples/GSE231559/GSM7290763")
s1 <- CreateSeuratObject(s1, project = "GSE231559-CRC1")

s2 <- Read10X("samples/GSE231559/GSM7290769")
s2 <- CreateSeuratObject(s2, project = "GSE231559-CRC2")

s3 <- Read10X("samples/GSE231559/GSM7290770")
s3 <- CreateSeuratObject(s3, project = "GSE231559-N1")

s4 <- Read10X("samples/GSE231559/GSM7290772")
s4 <- CreateSeuratObject(s4, project = "GSE231559-N2")

# Dataset 02 Samples
s5 <- Read10X("samples/GSE261012/GSM6432445")
s5 <- CreateSeuratObject(s5, project = "GSE261012-CRC1")

s6 <- Read10X("samples/GSE261012/GSM6432446")
s6 <- CreateSeuratObject(s6, project = "GSE261012-CRC2")

s7 <- Read10X("samples/GSE261012/GSM6432447")
s7 <- CreateSeuratObject(s7, project = "GSE261012-N1")

s8 <- Read10X("samples/GSE261012/GSM6432448")
s8 <- CreateSeuratObject(s8, project = "GSE261012-N2")

# Add metadata labels
s1$dataset <- "GSE231559";  s1$condition <- "CRC"
s2$dataset <- "GSE231559";  s2$condition <- "CRC"
s3$dataset <- "GSE231559";  s3$condition <- "Control"
s4$dataset <- "GSE231559";  s4$condition <- "Control"

s5$dataset <- "GSE261012";  s5$condition <- "CRC"
s6$dataset <- "GSE261012";  s6$condition <- "CRC"
s7$dataset <- "GSE261012";  s7$condition <- "Control"
s8$dataset <- "GSE261012";  s8$condition <- "Control"

# combine all samples
combined <- merge(s1, y = c(s2,s3,s4,s5,s6,s7,s8),
                  add.cell.ids = c("GSE231559-CRC1","GSE231559-CRC2","GSE231559-N1","GSE231559-N2",
                                   "GSE261012-CRC1","GSE261012-CRC2","GSE261012-N1","GSE261012-N2"))
# Verify the merge
table(combined$condition)
table(combined$orig.ident)
dim(combined)

# Remove intermediate objects
rm(s1,s2,s3,s4,s5,s6,s7,s8)

# ==============================================================================
# Perform Quality Control
# ==============================================================================

# calculate mitochondrial count
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^(?i)MT-")  # (?i) captures mt-Nd1, Mt-Nd1, and MT-Nd1
# Visualize QC metrics as a violin plot
Plot1 <- VlnPlot(combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 group.by = "condition", ncol = 3)
Plot1
ggsave(file = "raw_data.png", Plot1, dpi = 600, bg = "white")

# Filter data
filtered_data <- subset(combined, subset = nFeature_RNA > 200
                        & nFeature_RNA < 6000 & nCount_RNA > 2000 & percent.mt < 10)
dim(filtered_data) 

# visualize filtered data
filter_plt <- VlnPlot(filtered_data, features = c("nFeature_RNA", "nCount_RNA",
                                                  "percent.mt"),
                      group.by = "condition", ncol = 3)
filter_plt
ggsave(file = "filtered_data.png", filter_plt, dpi = 600, bg = "white")

# ==============================================================================
# Perform Normalization and calculate HVFs
# ==============================================================================
normalized_data <- NormalizeData(filtered_data, normalization.method = "LogNormalize", scale.factor = 10000)
normalized_data <- FindVariableFeatures(normalized_data, selection.method = "vst", nfeatures = 2000)

# The 10 most highly variable genes
top10 <- head(VariableFeatures(normalized_data), 10)
top10

# plot variable features with and without labels
plt1 <- VariableFeaturePlot(normalized_data)
plt2 <- LabelPoints(plot = plt1, points = top10, repel = TRUE)
variable_plt <- plt2
variable_plt
ggsave(file = "HVFs.png", variable_plt, dpi = 600, bg = "white")

# ==============================================================================
# Perform Scaling
# ==============================================================================
scaled_data <- ScaleData(normalized_data)
# Linear Dimensionality Reduction
scaled_data <- RunPCA(scaled_data, features = VariableFeatures(object = scaled_data))
# visualize PCA results
pca_plot <- DimPlot(scaled_data, reduction = "pca")
pca_plot
ggsave(file = "PCA_plot.png", pca_plot, dpi = 600, bg = "white")

# Determine Data Dimensions
elbow <- ElbowPlot(scaled_data, ndim = 50)
elbow
ggsave(file = "Elbow_plot.png", elbow, dpi = 600, bg = "white")

# rename data
data <- scaled_data
# remove intermediate object
rm(combined, filtered_data, normalized_data, scaled_data)

# ==============================================================================
# Cluster batch-effect uncorrected data
# ==============================================================================
data <- FindNeighbors(data, dims = 1:20, reduction = "pca")
data <- FindClusters(data, resolution = 0.6, algorithm = 4,
                         cluster.name = "raw_clusters")

# Run umap and tsne on batch-effect uncorrected data
data <- RunUMAP(data, dims = 1:20, reduction = "pca",
                    reduction.name = "umap.raw")

data <- RunTSNE(data, dims = 1:20, reduction = "pca",
                    reduction.name = "tsne.raw")

# Visualize batch-effect uncorrected data
# umap
umap_data <- DimPlot(data, reduction = "umap.raw")
umap_data
ggsave(file = "umap_raw_clusters.png", umap_data, dpi = 600, bg = "white")

# tsne
tsne_data <- DimPlot(data, reduction = "tsne.raw")
tsne_data
ggsave(file = "tsne_raw_clusters.png", tsne_data, dpi = 600, bg = "white")

# Visualize with conditions
# umap
umap_condition_data <- DimPlot(data, reduction = "umap.raw", group.by = "condition")
umap_condition_data
ggsave(file = "umap_raw_clusters-condition.png", umap_condition_data, dpi = 600, bg = "white")

# tsne
tsne_condition_data <- DimPlot(data, reduction = "tsne.raw", group.by = "condition")
tsne_condition_data
ggsave(file = "tsne_raw_clusters.png", tsne_condition_data, dpi = 600, bg = "white")

# ==============================================================================
# Perform Batch-Effect Correction between datasets using Harmony
# ==============================================================================
crt_data <- RunHarmony(data, group.by.vars = "dataset")

# save corrected data
saveRDS(crt_data, file = "Batch_effect_corrected_data.rds")

# ==============================================================================
# Cluster Batch-Effect Corrected Data
# ==============================================================================
# Determine Batch-Effect Corrected Data Dimensions
hm_elbow <- ElbowPlot(crt_data, ndim = 50)
hm_elbow
ggsave(file = "Batch-effect_corrected-Elbow_plot.png", hm_elbow, dpi = 600, bg = "white")

# Now find clusters for batch-effect corrected data
crt_data <- FindNeighbors(crt_data, reduction = "harmony", dims = 1:20)
crt_data <- FindClusters(crt_data, resolution = 0.6, algorithm = 4,
                         cluster.name = "harmony_clusters")

# Run umap and tsne for corrected data
crt_data <- RunUMAP(crt_data, reduction = "harmony",
                    dims = 1:20, reduction.name = "umap.harmony")
crt_data <- RunTSNE(crt_data, reduction = "harmony", dims = 1:20,
                    reduction.name = "tsne.harmony", reduction.key = "tsneharmony_")

# Visualize batch-effect corrected data
# umap
umap_crt <- DimPlot(crt_data, reduction = "umap.harmony")
umap_crt
ggsave(file = "umap_corrected-clusters.png", umap_crt, dpi = 600, bg = "white")

# tsne
tsne_crt <- DimPlot(crt_data, reduction = "tsne.harmony")
tsne_crt
ggsave(file = "tsne_corrected-clusters.png", tsne_crt, dpi = 600, bg = "white")

# Visualize with conditions
# umap
umap_crt_condition <- DimPlot(crt_data, reduction = "umap.harmony", group.by = "condition")
umap_crt_condition
ggsave(file = "umap_corrected-clusters-condition.png", umap_crt_condition, dpi = 600, bg = "white")

# tsne
tsne_crt_condition <- DimPlot(crt_data, reduction = "tsne.harmony", group.by = "condition")
tsne_crt_condition
ggsave(file = "tsne_corrected-clusters-condition.png", tsne_crt_condition, dpi = 600, bg = "white")

# Visualize data before and after batch-effect correction together
# Add titles to individual UMAP plots
umap_before <- umap_data + ggtitle("Clusters Before Batch Correction")
umap_after  <- umap_crt  + ggtitle("Clusters After Batch Correction")

# umap
combine_umap <- umap_before + umap_after
combine_umap
ggsave(file = "umap-before-and-after-batch-effect-correction.png", combine_umap, dpi = 600, bg = "white")

# Add titles to individual t-SNE plots
tsne_before <- tsne_data + ggtitle("Clusters Before Batch Correction")
tsne_after  <- tsne_crt  + ggtitle("Clusters After Batch Correction")

# tsne
combine_tsne <- tsne_before + tsne_after
combine_tsne
ggsave(file = "tsne-before-and-after-batch-effect-correction.png", combine_tsne, dpi = 600, bg = "white")

# ==============================================================================
# CELL TYPE ANNOTATION USING SCTYPE
# ==============================================================================

# Rename objects
raw_data <- data # Rename the "data" to "raw_data" to keep the uncorrected data
data <- crt_data # Rename the "crt_data" to "data" to proceed with further analysis
                 # Keep crt_data in the environment for backup

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
cL_resutls <- do.call("rbind", lapply(unique(data@meta.data$harmony_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$harmony_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$harmony_clusters==cl)), 10)
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
  data@meta.data$sctype_classification[data@meta.data$harmony_clusters == j] = as.character(cl_type$type[1])
}

# Verify visually by comparing unannotated and annotated tsne plots and matching with Celltypes.csv file
tsne_annt <- DimPlot(data, reduction = "tsne.harmony", label = TRUE,
                     repel = TRUE, group.by = 'sctype_classification')

verf_annt <- tsne_crt + tsne_annt
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
umap_annt <- DimPlot(data, reduction = "umap.harmony", label = FALSE,
                     repel = TRUE, group.by = 'sctype_classification')        
umap_annt
ggsave(file = "umap_annotated.png", umap_annt, dpi = 600, bg = "white")

tsne_annt <- DimPlot(data, reduction = "tsne.harmony", label = FALSE,
                     repel = TRUE, group.by = 'sctype_classification')        
tsne_annt
ggsave(file = "tsne_annotated.png", tsne_annt, dpi = 600, bg = "white")

umap_annt_condition <- DimPlot(data, reduction = "umap.harmony", label = FALSE,
                               repel = TRUE, group.by = 'sctype_classification',
                               split.by = "condition")        
umap_annt_condition
ggsave(file = "umap_annotated_condition.png", umap_annt_condition, dpi = 600, bg = "white")

tsne_annt_condition <- DimPlot(data, reduction = "tsne.harmony", label = FALSE,
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