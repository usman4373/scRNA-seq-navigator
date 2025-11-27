# ==============================================================================
# REQUIRED LIBRARIES FOR 3D VISUALIZATION
# ==============================================================================
library(Seurat)
library(plotly)
library(viridis)
library(htmlwidgets)
library(webshot2)

# ==============================================================================
# 3D UMAP VISUALIZATION
# ==============================================================================

# Load the dataset if not already
data <- readRDS("Annotated_data.rds")

# Generate 3D UMAP Coordinates
# Re-run UMAPs that you have accurate calculations for all UMAP(s)
data_umap <- RunUMAP(data, dims = 1:20, n.components = 3L)

# ==============================================================================
# Prepare UMAP Data for Visualization
# ==============================================================================

# Visualize what headings are called so that you can extract them to form a dataframe
umap_embeddings <- Embeddings(object = data_umap, reduction = "umap")
colnames(umap_embeddings)

# Prepare a dataframe for cell plotting and change the "RNA_snn_res.0.6" value according to your set resolution during clustering 
plot.data_umap <- FetchData(object = data_umap, vars = c("umap_1",
                                                         "umap_2",
                                                         "umap_3",
                                                         "RNA_snn_res.0.6"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data_umap$label <- paste(rownames(plot.data_umap))

# ==============================================================================
# Create Interactive 3D UMAP Plot by Cluster
# ==============================================================================

# Plot your data
umap_3d <- plot_ly(data = plot.data_umap, 
                   x = ~umap_1, y = ~umap_2, z = ~umap_3, 
                   color = ~RNA_snn_res.0.6,
                   colors = viridis::turbo(25),
                   type = "scatter3d", 
                   mode = "markers", 
                   marker = list(size = 2, width=2), # controls size of points
                   text=~label, # Extra column for cell ID in hover
                   hoverinfo="text") %>%
  layout(
    scene = list(
      xaxis = list(title = "UMAP 1", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      yaxis = list(title = "UMAP 2", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      zaxis = list(title = "UMAP 3", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      bgcolor = "black" # Black scene background
    ),
    paper_bgcolor = "black",  # Black plot background
    plot_bgcolor = "black",   # Black figure background
    font = list(color = "white"))
umap_3d

# ==============================================================================
# Export 3D UMAP Plot in Multiple Formats
# ==============================================================================

# Save as interactive HTML file
saveWidget(umap_3d, "3D_umap.html", selfcontained = TRUE)

# Convert HTML to PNG (high resolution)
webshot("3D_umap.html", file = "3D_umap.png", vwidth = 1200, vheight = 800)

# ==============================================================================
# 3D GENE EXPRESSION VISUALIZATION ON UMAP
# ==============================================================================

# Prepare Gene Expression Data for 3D Visualization
# Set default assay and create gene expression dataframe
DefaultAssay(object = data_umap) <- "RNA"

# Create a dataframe for gene expression visualization
plot.genedata <- FetchData(object = data_umap,
                           vars = c("umap_1",
                                    "umap_2",
                                    "umap_3",
                                    "ACTB"), # Add gene name
                           layer = 'scale.data')   # scale.data or data

# ==============================================================================
# Adjust Expression Scale for Better Visualization
# ==============================================================================

# Adjust the scale so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'expression' column of your dataframe
plot.genedata$expression <- ifelse(test = plot.genedata$ACTB <1, yes = plot.genedata$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plot.genedata$label <- paste(rownames(plot.genedata)," - ", plot.genedata$ACTB, sep="")

# ==============================================================================
# Create Interactive 3D Gene Expression Plot
# ==============================================================================

# Plot your data
gene_umap <- plot_ly(data = plot.genedata, 
                     x = ~umap_1, y = ~umap_2, z = ~umap_3,
                     color = ~expression,  # Expression values determine color scaling
                     opacity = 0.5, 
                     colors = c("#050C9C", "#E90074"),
                     type = "scatter3d",
                     mode = "markers", 
                     marker = list(size = 2, width = 2), 
                     text = ~label,  # Cell ID or gene label
                     hoverinfo = "text") %>%
  layout(
    scene = list(
      xaxis = list(title = "UMAP 1", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      yaxis = list(title = "UMAP 2", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      zaxis = list(title = "UMAP 3", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      bgcolor = "black"  # Black scene background
    ),
    paper_bgcolor = "black",  # Black plot background
    plot_bgcolor = "black",   # Black figure background
    font = list(color = "white") # White text
  )

gene_umap
# Save html
saveWidget(gene_umap, "ACTB-3D_umap.html", selfcontained = TRUE)
# Convert HTML to PNG (high resolution)
webshot("ACTB-3D_umap.html", file = "ACTB-3D_umap.png", vwidth = 1200, vheight = 800)

# ==============================================================================
# 3D T-SNE VISUALIZATION
# ==============================================================================

# Generate 3D t-SNE Coordinates

# Re-run tSNE so that you have accurate calculations for all tSNE(s)
data_tsne <- RunTSNE(data, reduction.use = "pca",
                     dims.use = 1:20, dim.embed = 3)

# Extract tSNE information from Seurat Object
tsne_1 <- data_tsne[["tsne"]]@cell.embeddings[,1]
tsne_2 <- data_tsne[["tsne"]]@cell.embeddings[,2]
tsne_3 <- data_tsne[["tsne"]]@cell.embeddings[,3]

# ==============================================================================
# Prepare t-SNE Data for Visualization
# ==============================================================================

# Prepare a dataframe for cell plotting and change the "RNA_snn_res.0.6" value according to your set resolution during clustering
plot.data_tsne <- FetchData(object = data_tsne, vars = c("tSNE_1",
                                                         "tSNE_2",
                                                         "tSNE_3",
                                                         "RNA_snn_res.0.6"))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data_tsne$label <- paste(rownames(plot.data_tsne))

# ==============================================================================
# Create Interactive 3D t-SNE Plot by Cluster
# ==============================================================================

# Plot your data
tsne_3d <- plot_ly(data = plot.data_tsne, 
                   x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3,
                   color = ~RNA_snn_res.0.6, 
                   colors = viridis::turbo(25), 
                   type = "scatter3d",
                   mode = "markers", 
                   marker = list(size = 2, width=2), # Controls size of points
                   text = ~label, # Extra column for cell ID in hover
                   hoverinfo = "text") %>%
  layout(
    scene = list(
      xaxis = list(title = "tSNE 1", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      yaxis = list(title = "tSNE 2", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      zaxis = list(title = "tSNE 3", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
      bgcolor = "black" # Black scene background
    ),
    paper_bgcolor = "black",  # Black plot background
    plot_bgcolor = "black",   # Black figure background
    font = list(color = "white") # White text
  )
tsne_3d
# Save html
saveWidget(tsne_3d, "3D_tsne.html", selfcontained = TRUE)
# Convert HTML to PNG (high resolution)
webshot("3D_tsne.html", file = "3D_tsne.png", vwidth = 1200, vheight = 800)

# ==============================================================================
# 3D GENE EXPRESSION VISUALIZATION ON T-SNE
# ==============================================================================

# Prepare Gene Expression Data for 3D t-SNE Visualization
# Create a dataframe for gene expression on t-SNE
plotting.genedata <- FetchData(object = data_tsne, vars = c("tSNE_1",
                                                            "tSNE_2",
                                                            "tSNE_3",
                                                            "ACTB"),
                                                  layer = "scale.data")

# ==============================================================================
# Adjust Expression Scale for t-SNE Visualization
# ==============================================================================

# Adjust the scale so that every cell having an expression >1 will be one color
# Basically, you are re-adjusting the scale here, so that any cell having a certain expression 
# will light up on your 3D plot

# First make another column in your dataframe, where all values above 1 are re-assigned a value of 1
# This information is stored in the 'expression' column of your dataframe
plotting.genedata$expression <- ifelse(test = plotting.genedata$ACTB <1, yes = plotting.genedata$ACTB, no = 1)

# Add the label column, so that now the column has 'cellname-its expression value'
plotting.genedata$label <- paste(rownames(plotting.genedata)," - ", plotting.genedata$ACTB, sep="")

# ==============================================================================
# Create Interactive 3D Gene Expression Plot on t-SNE
# ==============================================================================

# Plot your data
gene_tsne <- plot_ly(data = plotting.genedata, x = ~tSNE_1, y = ~tSNE_2, z = ~tSNE_3,
                     color = ~expression, # you can just run this against the column for the gene as well using ~ACTB, the algorith will automatically scale in that case based on maximal and minimal values
                     opacity = .5, colors = c("#050C9C", "#E90074"), type = "scatter3d",
                     mode = "markers", marker = list(size = 2, width=2),
                     text=~label, hoverinfo="text") %>%
  layout(scene = list(xaxis = list(title = "tSNE 1", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
                      yaxis = list(title = "tSNE 2", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
                      zaxis = list(title = "tSNE 3", showgrid = TRUE, gridcolor = "white", zeroline = FALSE, color = "white"),
                      bgcolor = "black" # Black scene background
  ),
  paper_bgcolor = "black",  # Black plot background
  plot_bgcolor = "black",   # Black figure background
  font = list(color = "white") # White text
  )
gene_tsne
# Save html
saveWidget(gene_tsne, "ACTB-3d-tsne.html", selfcontained = TRUE)
# Convert HTML to PNG (high resolution)
webshot("ACTB-3d-tsne.html", file = "ACTB-3D_tsne.png", vwidth = 1200, vheight = 800)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================