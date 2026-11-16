library(cowplot)

sce <- Reduce(function(x, y) merge(x, y), ifnb.list)
sce=JoinLayers(sce)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce,selection.method = "vst", nfeatures = 3000)
sce <- ScaleData(sce,verbose = T)
sce <- RunPCA(sce,npcs = 50, verbose = T)
scedatam <- RunHarmony(sce, "orig.ident")
GetPC(scedatam,reduction = "harmony")
scedatam <- RunUMAP(scedatam, reduction = "harmony",dims = 1:41)
scedatam <- FindNeighbors(scedatam, reduction = "harmony", dims = 1:41)
scedatam <- FindClusters(scedatam,resolution = 0.4)
DimPlot(scedatam, reduction = "umap", cols = colors,label = T)
DotPlot(scedatam,c("Slc44a5"))
scedatam@meta.data$group[scedatam@meta.data$orig.ident %in% c("ADR1","ADR2","ADR3")]<-"ADRN";
scedatam@meta.data$group[scedatam@meta.data$orig.ident %in% c("ORG1","ORG2","ORG3")]<-"ORG";
scedatam@meta.data$group[scedatam@meta.data$orig.ident %in% c("Con1_BBC","Con2_BBC","Con3_BBC")]<-"Con_BBC";
scedatam@meta.data$group[scedatam@meta.data$orig.ident %in% c("Con1_C57","Con2_C57","Con3_C57")]<-"Con_C57";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("14")]<-"ATL";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("13","22")]<-"CNT";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("7")]<-"DCT";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("1","25")]<-"EC"
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("20")]<-"ICA";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("18")]<-"ICB";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("19")]<-"MAC";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("15")]<-"PC";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("23")]<-"POD"
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("0","3","4","6","8","9","10","12","16","21","27")]<-"PT";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("5","11")]<-"SC";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("26")]<-"T";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("2","24","28")]<-"TAL";
scedatam@meta.data$celltype[scedatam@meta.data$seurat_clusters %in% c("17")]<-"TL";
DimPlot(scedatam, reduction = "umap", cols = colors,label = T,group.by = "celltype")
--------------------------------------------------------------------------------
Idents(scedatam)="celltype"
PODm=subset(scedatam,idents="POD")
ExportGeneExpr<-GetAssayData(PODm,assay="RNA",slot="counts")
--------------------------------------------------------------------------------
colors <- c(
  "#54B0E4",'#ead2bd',"#a4cde1", "#4DAF4A", "#F781BF",
    "#BC9DCC", "#A65628", "#f9b769","#e32427", "#D2EBC8",'#d9d9d9',"#F5D2A8",
    "#F5CFE4","#BBDD78")
scPlot(scedatam, group_by = "celltype", point_size = 0.3,
        label_size = 8,
        circle_size = 18,circle_alpha = 0.6,colors = colors,
        legend_cols = 2, raster_dpi = 600, arrow_len = 0.15,
        arrow_frac = 0.15, axis_label_size = 5) +
  theme(legend.position = "none")
Idents(scedatam)="celltype"
PODm=subset(scedatam,idents="POD")
GetPC(PODm,reduction = "harmony")
PODm <- RunTSNE(PODm, reduction = "harmony",dims = 1:36)
PODm <- FindNeighbors(PODm, reduction = "harmony", dims = 1:36)
PODm <- FindClusters(PODm,resolution = 0.2)
DimPlot(PODm, reduction = "tsne", cols = colors,label = T)
DotPlot(PODm,"Wt1")
markers_plot <- c("Wt1", "Nphs1", "Nphs2", "Podxl")
plist2 <- vector("list", length(markers_plot))
for (i in seq_along(markers_plot)) {
  p <- single_CM_density(
    object = PODm,
    pt.size=1.5,
    features = markers_plot[i],
    night = TRUE
  )
  plist2[[i]] <- p
}
plot_grid(plotlist = plist2, ncol = 2)
PODm@meta.data$celltype[PODm@meta.data$seurat_clusters %in% c("0")]<-"hPOD";
PODm@meta.data$celltype[PODm@meta.data$seurat_clusters %in% c("1")]<-"iPOD1";
PODm@meta.data$celltype[PODm@meta.data$seurat_clusters %in% c("2")]<-"iPOD2";
PODm <- RunSlingshot(
  PODm,
  group.by = "celltype",
  reduction = "tsne",
  start = "hPOD")
scplotter::CellDimPlot(
  PODm,
  group_by = "celltype",
  reduction = "tsne",
  pt_size = 2,
  palcolor = c("#b2df8a","#E31A1C","#FB9A99"),
  lineages = paste0("Lineage", 1),
  theme = "theme_blank",
  lineages_span = 0.5,
  lineages_palcolor = "#B383B9",
  lineages_linewidth = 1.5,
)
FeatureDimPlot(
  PODm,
  features = paste0("Lineage", 1),
  pt.size = 2,
  reduction = "tsne",
  theme_use = "theme_blank",
  xlab = "tSNE_1",
  ylab = "tSNE_2"
)
DimPlot(scereal,reduction = "umap")
