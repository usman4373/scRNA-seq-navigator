# ==============================================================================
# Load Required Libraries
# ==============================================================================
library(EnhancedVolcano)
library(ggplot2)
library(readr)

# ==============================================================================
# DATA LOADING AND PREPARATION
# ==============================================================================
file_path <- "Cancer cells_DEGs.csv"
degs <- read_csv(file_path, show_col_types = FALSE)
degs <- as.data.frame(degs)              
rownames(degs) <- degs$gene              
degs$gene <- NULL                        

# ==============================================================================
# DATA CHECK
# ==============================================================================
print(head(degs[, c("avg_log2FC", "fdr")]))

# ==============================================================================
# VOLCANO PLOT GENERATION
# ==============================================================================
volc <- EnhancedVolcano(
  degs,
  lab              = NA,
  x                = "avg_log2FC",
  y                = "fdr",
  pCutoff          = 0.05,
  FCcutoff         = 2,  # Adjust log FC value cutoff
  title            = "Cancer cells-LM vs All Other Cell Types-Normal",   # Change according to your celltype degs file
  pointSize        = 2.0,
  labSize          = 3.0,
  legendLabSize    = 10,
  legendPosition    = "right",
  legendIconSize   = 4.0,
  drawConnectors   = FALSE,
  widthConnectors  = 0.5,
  gridlines.major  = TRUE,
  gridlines.minor  = TRUE,
  col = c(
    'grey50',             # Insignificant p-value and log2FC
    '#1DCD9F',            # Insignificant log2FC
    'royalblue',          # Insignificant p-value
    '#FF0B55'             # Significant p-value and log2FC
  ))

# ==============================================================================
# Display and Save Volcano Plot
# ==============================================================================
volc2 <- volc +
  theme_bw(base_size = 14) +
  theme(
    # ---- REMOVE PANEL BORDER (this is the real culprit) ----
    panel.border = element_blank(),
    
    # ---- RE-DRAW ONLY LEFT + BOTTOM SPINES ----
    axis.line = element_line(color = "black"),
    axis.line.x.top    = element_blank(),
    axis.line.y.right  = element_blank(),
    
    # ---- REMOVE TOP/RIGHT TICKS TOO ----
    axis.ticks.x.top   = element_blank(),
    axis.ticks.y.right = element_blank(),
    
    # ---- GRIDLINES ----
    panel.grid.major = element_line(color = "grey80", size = 0.6),
    #panel.grid.minor = element_line(color = "grey90", size = 0.6), # Uncomment to add minor gridlines
    
    # ---- REMOVE LEGEND TITLE ----
    legend.title = element_blank(),
    
    # ---- BOLD MAIN TITLE ----
    plot.title = element_text(face = "bold", size = 16),
    
    # ---- CLEAN BACKGROUND ----
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA)
  )

volc2
ggsave("Cancer cells_Volcano.png", volc2, dpi = 600)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================