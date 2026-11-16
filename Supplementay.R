library(pheatmap)
scPlot1(scedata, group_by = "orig.ident", point_size = 0.3,
        label_size = 0,
        circle_size = 18,circle_alpha = 0.6,color=colors,
        legend_cols = 2, raster_dpi = 600, arrow_len = 0.15,
        arrow_frac = 0.15, axis_label_size = 5) +
  theme(legend.position = "right")
scPlot1(scedata, group_by = "disease", point_size = 0.3,
        label_size = 0,
        circle_size = 18,circle_alpha = 0.6,color=colors,
        legend_cols = 1, raster_dpi = 600, arrow_len = 0.15,
        arrow_frac = 0.15, axis_label_size = 5) +
  theme(legend.position = "right")
scPlot1(scedatam, group_by = "orig.ident", point_size = 0.3,
        label_size = 0,
        circle_size = 18,circle_alpha = 0.6,color=colors,
        legend_cols = 2, raster_dpi = 600, arrow_len = 0.15,
        arrow_frac = 0.15, axis_label_size = 5) +
  theme(legend.position = "right")
scPlot1(scedatam, group_by = "group", point_size = 0.3,
        label_size = 0,
        circle_size = 18,circle_alpha = 0.6,color=colors,
        legend_cols = 1, raster_dpi = 600, arrow_len = 0.15,
        arrow_frac = 0.15, axis_label_size = 5) +
  theme(legend.position = "right")
marker_list=c("SLC8A1","SLC12A3","EMCN","SLC4A1","SLC4A9","CD163","NOS1","RBFOX1","NPHS2","LRP2","PDGFRB","CD247","SLC12A1","SLC44A5")
p <- FeatureDimPlot(
  scedata,
  features   = marker_list,
  reduction  = "UMAP",
  theme_use  = "theme_blank",
  theme_args = list(base_size = 18),
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  ncol = 7,
  legend.position = "right"
)
  filename = "E:/OB/paper//FeatureDimPlot_UMAP.jpg",
  plot = p,
  width = 13.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  quality = 100
)
marker_list1=c("Atp10b","Slc8a1","Slc12a3","Emcn","Slc4a1","Slc4a9","Apobec1",
               "Aqp2","Nphs2","Lrp2","Pdgfrb","Cd247","Slc12a1","Cryab")
p <- FeatureDimPlot(
  scedatam,
  features   = marker_list1,
  reduction  = "UMAP",
  theme_use  = "theme_blank",
  theme_args = list(base_size = 18),
  xlab = "UMAP_1",
  ylab = "UMAP_2",
  ncol = 7,
  legend.position = "right"
)
  filename = "E:/OB/paper//MiceFeatureDimPlot_UMAP.jpg",
  plot = p,
  width = 13.5,
  height = 4.5,
  units = "in",
  dpi = 600,
  quality = 100
)
-------------
scDot(
    POD,
    markers=c("WT1","NPHS1","NPHS2","SYNPO","PODXL"),
    celltype_var = "celltype",
    layer = "data",
    assay = NULL,
    color = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0"),
    dot_size_range = c(4, 15),
    standardize = TRUE,
    scale_range = c(-1, 1),
    flip = TRUE,
    theme = 1
  )
scDot(
  PODm,
  markers=c("Wt1","Nphs1","Nphs2","Synpo","Podxl"),
  celltype_var = "celltype",
  layer = "data",
  assay = NULL,
  color = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0"),
  dot_size_range = c(4, 15),
  standardize = TRUE,
  scale_range = c(-1, 1),
  flip = TRUE,
  theme = 1
)
cellRatioPlot(
  POD,
  sample.name = "orig.ident",
  celltype.name = "celltype",
  fill.col= c("#b2df8a","#FB9A99","#E31A1C")
) +
  theme(
    axis.title = element_text(
      size = 22,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 16,
      color = "black"
    ),
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    legend.text = element_text(
      size = 16
    )
  )
cellRatioPlot(
  PODm,
  sample.name = "orig.ident",
  celltype.name = "celltype",
  fill.col= c("#b2df8a","#FB9A99","#E31A1C")
) +
  theme(
    axis.title = element_text(
      size = 22,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      size = 16,
      color = "black"
    ),
    legend.title = element_text(
      size = 18,
      face = "bold"
    ),
    legend.text = element_text(
      size = 16
    )
  )
scDot(
  POD,
  markers=c("CALR","HMGB1","HSP90AA1","HSPA1A","HSAPA1B","HSPA1L"),
  celltype_var = "disease",
  layer = "data",
  assay = NULL,
  color = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0"),
  dot_size_range = c(4, 15),
  standardize = TRUE,
  scale_range = c(-1, 1),
  flip = TRUE,
  theme = 1
)
scDot(
  POD,
  markers=c("GK"),
  celltype_var = "celltype",
  layer = "data",
  assay = NULL,
  color = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0"),
  dot_size_range = c(4, 15),
  standardize = TRUE,
  scale_range = c(-1, 1),
  flip = TRUE,
  theme = 1
)
scDot(
  POD,
  markers=c("GK"),
  celltype_var = "celltype",
  layer = "data",
  assay = NULL,
  color = c("#FFFFCC", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0"),
  dot_size_range = c(4, 15),
  standardize = TRUE,
  scale_range = c(-1, 1),
  flip = TRUE,
  theme = 1
)
---
rownames(meta) <- meta$Metabolite.Name
meta_M3 <- meta[grep("M\\+3$", meta$Metabolite.Name), ]
heatmap_M3 <- meta_M3[, -1]
heatmap_M3 <- as.matrix(heatmap_M3)
heatmap_M3 <- heatmap_M3[
  apply(heatmap_M3, 1, sd) > 0,
]
pheatmap(
  heatmap_M3,
  scale = "row",
  cluster_rows = T,
  cluster_cols = F,
  border_color = "white",
  fontsize_row = 15,
  fontsize_col = 15
)
