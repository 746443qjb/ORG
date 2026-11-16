library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggwordcloud)
library(patchwork)
library(stringr)
library(simplifyEnrichment)
library(openxlsx)

sce <- Reduce(function(x, y) merge(x, y), ifnb.list)
sce=JoinLayers(sce)
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce,selection.method = "vst", nfeatures = 3000)
sce <- ScaleData(sce,verbose = T)
sce <- RunPCA(sce,npcs = 50, verbose = T)
scedat
scedata <- subset(scedata, idents = 11, invert = TRUE)
scedata <- RunHarmony(sce, "orig.ident")
scedata <- RunUMAP(scedata, reduction = "harmony",dims = 1:37)
scedata <- FindNeighbors(scedata, reduction = "harmony", dims = 1:37)
scedata <- FindClusters(scedata,resolution = 0.3)
Idents(scedata)="seurat_clusters"
DimPlot(scedata, reduction = "umap",cols= colors,label = T, group.by = "celltype")
Idents(scedata)="celltype"
colors <- c(
  "#FB9A99", "#377EB8", "#4DAF4A", "#F781BF", "#BC9DCC", "#A65628", "#54B0E4", "#F5D2A8",
  "#222F75","#B2DF8A", "#E3BE00", "#E7298A", "#E41A1C", "#D2EBC8", "#00CDD1",
  "#5E4FA2", "#8CA77B", "#1B9E77","#8DD3C7", "#7DBFC7", "#B3DE69", "#999999", "#FCED82",
  "#F5CFE4", "#F5D2A8", "#B383B9", "#EE934E", "#BBDD78"
)
----
scedata@meta.data$disease[scedata@meta.data$orig.ident %in% c("Con1","Con2","Con3")]<-"CON";
scedata@meta.data$disease[scedata@meta.data$orig.ident %in% c("OBN1","OBN2","OBN3")]<-"ORG";
scedata@meta.data$disease[scedata@meta.data$orig.ident %in% c("MCD1","MCD2","MCD3")]<-"MCD";
scedata@meta.data$disease[scedata@meta.data$orig.ident %in% c("FSGS1","FSGS2","FSGS3")]<-"FSGS";
--------------------------------------------------------------------------------
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("17")]<-"POD";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("0","1","8","10","18")]<-"PT";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("22")]<-"TL";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("12")]<-"PEC";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("14")]<-"MAC";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("16")]<-"T";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("3","6","11")]<-"TAL";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("20")]<-"MD";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("2")]<-"DCT";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("4","21")]<-"CT";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("9")]<-"SC";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("15","5")]<-"EC";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("13")]<-"ICB";
scedata@meta.data$celltype[scedata@meta.data$seurat_clusters %in% c("7","19")]<-"ICA";
marker=FindAllMarkers(scedata,logfc.threshold = 1,
                      test.use = "wilcox",
                      only.pos = T,
                      slot = "data",
                      min.pct = 0.3)
Idents(scedata)="celltype"
DimPlot(scedata_, reduction = "umap",label = T)
--------------------------------------------------------------------------------
Idents(scedata)="celltype"
POD=subset(scedata,idents="POD")
ExportGeneExpr<-GetAssayData(POD,assay="RNA",slot="counts")
--------------------------------------------------------------------------------
colors <- c(
  '#ead2bd',"#a4cde1", "#4DAF4A", "#F781BF",
  "#BC9DCC", "#A65628", "#54B0E4", "#f9b769","#e32427", "#D2EBC8",'#d9d9d9',"#F5D2A8",
"#F5CFE4","#BBDD78")
colors <- c('#4b6aa8','#3ca0cf','#c376a7','#ad98c3','#cea5c7',
            '#53738c','#a78982','#92699e',
            '#d69971','#df5734','#6c408e','#ac6894','#d4c2db',
            '#537eb7','#83ab8e','#ece399','#405993','#cc7f73',
            '#b95055','#d5bb72','#bc9a7f','#e0cfda','#d8a0c0',
            '#e6b884','#b05545','#d69a55','#64a776','#cbdaa9',
            '#efd2c9','#da6f6d','#ebb1a4','#a44e89','#a9c2cb',
            '#b85292','#6d6fa0','#8d689d','#c8c7e1','#d25774',
            '#c49abc','#927c9a','#3674a2','#9f8d89','#72567a',
            '#63a3b8','#c4daec','#61bada','#b7deea','#e29eaf',
            '#4490c4','#e6e2a3','#de8b36','#c4612f','#9a70a8',
            '#76a2be','#408444','#c6adb0','#9d3b62','#2d3462')
colorlist <- c("#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43",
               "#b79973","#f38989","#ec5051",,"#ef6a45","#f9a341","#f48521",
               "#ee8e46","#d4a6a8","#af93c4","#8660a8","#815e99","#c6b598","#f6f28f","#d4a55b",
               "#b05a28")
colors1<-c('#85d2e1','#e0d200','#eb5b5e',,'#749f87',
           '#e38501','#879da2','#02401b','#0073c2','#b40f1f',
           '#c5cdf7','#7092d2','#e69fc4','#b6854d')
colors2<-c('#8dd1c6',,'#bcb9d8','#80b2d4','#feb662',
           '#b4df65',,'#be7fbf','#cbecc4','#ffee72','#e6bdcb')
colors3<-c('#306311','#df772b','#e184ca',,'#e8e684',
           '#b2dfe8','#667570','#5ec9b8','#314e9b','#616e88',
           '#7623f1','#eaeaea')
scPlot1(scedata, group_by = "orig.ident", point_size = 0.3,
        label_size = 0,
        circle_size = 18,circle_alpha = 0.6,color=colors,
        legend_cols = 2, raster_dpi = 600, arrow_len = 0.15,
        arrow_frac = 0.15, axis_label_size = 5) +
  theme(legend.position = "right")
--------------------------------------------------------------------------------
POD=subset(scedata,idents = "POD")
POD <- NormalizeData(POD)
POD <- FindVariableFeatures(POD,selection.method = "vst", nfeatures = 3000)
POD <- ScaleData(POD,verbose = T)
POD <- RunPCA(POD,npcs = 50, verbose = T)
GetPC(POD,reduction = "harmony")
table(scedata$celltype)
POD <- RunTSNE(POD, reduction = "harmony",dims = 1:42)
POD <- FindNeighbors(POD, reduction = "harmony", dims = 1:42)
POD <- FindClusters(POD,resolution = 0.2)
DimPlot(POD, reduction = "tsne", cols = colors,label = T)
POD@meta.data$celltype[POD@meta.data$seurat_clusters %in% c("0")]<-"hPOD";
POD@meta.data$celltype[POD@meta.data$seurat_clusters %in% c("1")]<-"iPOD1";
POD@meta.data$celltype[POD@meta.data$seurat_clusters %in% c("2")]<-"iPOD2";
DimPlot(POD,group.by = "celltype",reduction = "tsne")
scPlot1(
  POD,
  group_by = "celltype",
  reduction = "tsne",
  colors=c("#b2df8a","#FB9A99","#E31A1C"),
  point_size = 3,
  label_size = 7,
  legend_cols = 2,
  raster_dpi = 600
) +
  theme(legend.position = "none")
Idents(POD)="celltype"
markers_plot <- c("WT1", "NPHS1", "NPHS2", "PODXL")
plist2 <- vector("list", length(markers_plot))
for (i in seq_along(markers_plot)) {
  p <- single_CM_density(
    object = POD,
    pt.size=1.5,
    features = markers_plot[i],
    night = TRUE
  )
  plist2[[i]] <- p
}
plot_grid(plotlist = plist2, ncol = 2)
scDot(PODm,markers=c("Wt1", "Nphs1"))
POD <- FindClusters(POD,resolution = seq(from = 0,to = 0.1,by = 0.01))
ClustreePlot(POD, prefix = "RNA_snn_res.",aspect.ratio=1.5,theme_args = list(
  text = element_text(size = 20),
  plot.tag = element_text(size = 50, face = "bold"),
  axis.text.y = element_text(size = 30),
  axis.title = element_text(size = 20),
  plot.title = element_text(size = 16),
  legend.text = element_text(size = 30),
  legend.title = element_text(size = 30)
))
cellRatioPlot(POD,sample.name = "disease",celltype.name = "celltype")
--------------------------------------------------------------------------------
POD <- RunSlingshot(
  POD,
  group.by = "celltype",
  reduction = "tsne",
  start = "hPOD")
CellDimPlot(
  POD,
  group_by = "celltype",
  reduction = "tsne",
  pt_size = 2,
  palcolor = c("#FB9A99","#b2df8a","#E31A1C"),
  lineages = paste0("Lineage", 1),
  theme = "theme_blank",
  lineages_span = 0.5,
  lineages_palcolor = "#B383B9",
  lineages_linewidth = 1.5,
)
FeatureDimPlot(
  POD,
  features = paste0("Lineage", 1),
  pt.size = 2,
  reduction = "tsne",
  theme_use = "theme_blank",
  xlab = "tSNE_1",
  ylab = "tSNE_2"
)
--------------------------------------------------------------------------------
Idents(POD)="celltype"
markers=FindAllMarkers(POD,
                    logfc.threshold = 0,
                    test.use = "MAST",
                    only.pos = T,
                    min.pct=0)
Idents(PODm)="celltype"
markers_m=FindAllMarkers(PODm,
                       logfc.threshold = 0,
                       test.use = "MAST",
                       only.pos = T,
                       min.pct=0)
head(markers,n=5)
head(markers_m,n=5)
----
human_hPOD <- markers %>%
  filter(cluster == "hPOD", p_val_adj < 0.05) %>%
  mutate(gene = toupper(gene)) %>%
  group_by(gene) %>%
  summarise(
    avg_log2FC_human = mean(avg_log2FC),
    p_val_adj_human = min(p_val_adj),
    .groups = "drop"
  )
mouse_hPOD <- markers_m %>%
  filter(cluster == "hPOD", p_val_adj < 0.05) %>%
  mutate(gene = toupper(gene)) %>%
  group_by(gene) %>%
  summarise(
    avg_log2FC_mouse = mean(avg_log2FC),
    p_val_adj_mouse = min(p_val_adj),
    .groups = "drop"
  )
shared_hPOD <- inner_join(human_hPOD, mouse_hPOD, by = "gene") %>%
  arrange(desc(avg_log2FC_human + avg_log2FC_mouse))
head(shared_iPOD2)
n_shared <- nrow(shared_iPOD2)
n_shared
model <- lm(avg_log2FC_mouse ~ avg_log2FC_human, data = shared_iPOD2)
r2 <- summary(model)$r.squared
pval <- summary(model)$coefficients[2, 4]
label_text <- sprintf(
  "R² = %.3f\np = %.2e\nShared genes = %d",
  r2, pval, n_shared
)
shared_iPOD2 <- shared_iPOD2 %>%
  mutate(
    avg_log2FC_mean = (avg_log2FC_human + avg_log2FC_mouse) / 2
  )
ggplot(
  shared_iPOD2,
  aes(x = avg_log2FC_human, y = avg_log2FC_mouse)
) +
  geom_point(
    aes(color = avg_log2FC_mean),
    alpha = 0.9,
    size = 3.2
  ) +
  scale_color_distiller(
    palette = "Spectral",
    direction = -1,
    name = "Mean log2FC"
  ) +
  geom_smooth(
    method = "lm",
    color = "#E31A1C",
    fill = "#E31A1C",
    se = TRUE,
    alpha = 0.3,
    linewidth = 1.5
  ) +
  geom_text_repel(
    data = label_data,
    aes(
      x = avg_log2FC_human,
      y = avg_log2FC_mouse,
      label = gene
    ),
    inherit.aes = FALSE,
    size = 4.2,
    fontface = "bold",
    color = "black",
    box.padding = 0.5,
    point.padding = 0.35,
    segment.color = "grey35",
    segment.linewidth = 0.4,
    max.overlaps = Inf
  ) +
  annotate(
    "text",
    x = Inf,
    y = -Inf,
    hjust = 1.1,
    vjust = -0.8,
    label = label_text,
    size = 9,
    fontface = "bold",
    color = "black"
  ) +
  labs(
    x = "Human log2FC (iPOD2)",
    y = "Mouse log2FC (iPOD2)"
  ) +
  theme_bw(base_size = 26) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    legend.position = "none"
  )+
  scale_x_continuous(
    labels = function(x) sprintf("%.1f", x)
  ) +
  scale_y_continuous(
    labels = function(x) sprintf("%.1f", x)
  )
--------------------------------------------------------------------------------
gene_list <- list(
  hPOD  = unique(shared_hPOD$gene),
  iPOD1 = unique(shared_iPOD1$gene),
  iPOD2 = unique(shared_iPOD2$gene)
)
max_len <- max(lengths(gene_list))
gene_df <- as.data.frame(
  lapply(gene_list, function(x) {
    length(x) <- max_len
    return(x)
  })
)
run_go_bp <- function(genes,
                      group_name,
                      OrgDb = org.Hs.eg.db,
                      pvalueCutoff = 1,
                      qvalueCutoff = 1) {
  genes <- unique(toupper(genes))
  genes <- genes[!is.na(genes)]
  gene_df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = OrgDb
  )
  if (nrow(gene_df) == 0) {
    warning(group_name, ": no gene can be converted to ENTREZID.")
    return(NULL)
  }
  ego <- enrichGO(
    gene          = gene_df$ENTREZID,
    OrgDb         = OrgDb,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = pvalueCutoff,
    qvalueCutoff  = qvalueCutoff,
    readable      = TRUE
  )
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    warning(group_name, ": no enriched GO BP terms.")
    return(NULL)
  }
  ego_df <- as.data.frame(ego)
  ego_df$Groups <- group_name
  ego_df$Database <- "GO_BP"
  return(
    list(
      enrich_result = ego,
      enrich_df = ego_df,
      gene_map = gene_df
    )
  )
}
go_hPOD <- run_go_bp(
  genes = gene_list$hPOD,
  group_name = "hPOD"
)
go_iPOD1 <- run_go_bp(
  genes = gene_list$iPOD1,
  group_name = "iPOD1"
)
go_iPOD2 <- run_go_bp(
  genes = gene_list$iPOD2,
  group_name = "iPOD2"
)
go_df_all <- bind_rows(
  go_hPOD$enrich_df,
  go_iPOD1$enrich_df,
  go_iPOD2$enrich_df
)

p_go_wordcloud_all <- ggplot(
  word_df_all,
  aes(
    label = keyword,
    size = count,
    color = score,
    angle = angle
  )
) +
  ggwordcloud::geom_text_wordcloud(
    rm_outside = TRUE,
    eccentricity = 1,
    shape = "square",
    show.legend = TRUE,
    grid_margin = 3
  ) +
  scale_color_distiller(
    palette = "Spectral",
    direction = -1,
    name = "Score\n-log10(adj.P)"
  ) +
  scale_size(
    name = "Count",
    range = c(3, 10)
  ) +
  facet_wrap(
    ~ Groups,
    nrow = 1
  ) +
  coord_flip() +
  theme_bw(base_size = 18) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    strip.background = element_rect(
      fill = "#E6F5E6",
      color = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 18
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right"
  )
p_go_wordcloud_all

word_df_all %>%
  dplyr::filter(Groups == "hPOD") %>%
  dplyr::arrange(dplyr::desc(score)) %>%
  dplyr::select(Groups, keyword, score, count)

ggplot(
  word_hPOD,
  aes(
    label = keyword,
    size = score,
    color = score,
    angle = angle
  )
) +
  ggwordcloud::geom_text_wordcloud(
    rm_outside = TRUE,
    eccentricity = 1,
    shape = "square",
    show.legend = TRUE,
    grid_margin = 3
  ) +
  scale_color_distiller(
    palette = "Spectral",
    direction = -1,
    name = "Score\n-log10(adj.P)"
  ) +
  scale_size(
    name = "Score",
    range = c(4, 12)
  ) +
  coord_flip() +
  labs(
    caption = "hPOD-GO-BP "
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0.5,
      size = 20,
      margin = margin(t = 8)
    )
  )

ggplot(
  word_iPOD1,
  aes(
    label = keyword,
    size = score,
    color = score,
    angle = angle
  )
) +
  ggwordcloud::geom_text_wordcloud(
    rm_outside = TRUE,
    eccentricity = 1,
    shape = "square",
    show.legend = TRUE,
    grid_margin = 3
  ) +
  scale_color_distiller(
    palette = "Spectral",
    direction = -1,
    name = "Score\n-log10(adj.P)"
  ) +
  scale_size(
    name = "Score",
    range = c(4, 12)
  ) +
  coord_flip() +
  labs(
    caption = "iPOD1-GO-BP "
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0.5,
      size = 20,
      margin = margin(t = 8)
    )
  )

ggplot(
  word_iPOD2,
  aes(
    label = keyword,
    size = score,
    color = score,
    angle = angle
  )
) +
  ggwordcloud::geom_text_wordcloud(
    rm_outside = TRUE,
    eccentricity = 1,
    shape = "square",
    show.legend = TRUE,
    grid_margin = 3
  ) +
  scale_color_distiller(
    palette = "Spectral",
    direction = -1,
    name = "Score\n-log10(adj.P)"
  ) +
  scale_size(
    name = "Score",
    range = c(4, 12)
  ) +
  coord_flip() +
  labs(
    caption = "iPOD2-GO-BP "
  ) +
  theme_bw(base_size = 20) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    legend.position = "right",
    plot.caption = element_text(
      hjust = 0.5,
      size = 20,
      margin = margin(t = 8)
    )
  )

