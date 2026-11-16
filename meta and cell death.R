library(clusterProfiler)
library(biomaRt)
library(babelgene)
POD_daixie=sc.metabolism.Seurat.2(POD, method = "VISION", imputation = F, ncores = 2,
  geneList = "KEGG")
data=POD_daixie@assays[["METABOLISM"]][["score"]]
colnames(data) <- gsub("\\.1_", "-1_", colnames(data))
head(colnames(data))
head(Cells(POD))
length(intersect(colnames(data), Cells(POD)))
pathways_use <- c(
  "Fatty acid degradation",
  "Citrate cycle (TCA cycle)",
  "Oxidative phosphorylation"
)
pathways_use %in% rownames(data)
score_mat <- data[pathways_use, , drop = FALSE]
score_df <- as.data.frame(t(score_mat))
colnames(score_df) <- c(
  "FA_degradation",
  "TCA_cycle",
  "Oxidative_phosphorylation"
)
all(rownames(score_df) %in% Cells(POD))
all(Cells(POD) %in% rownames(score_df))
score_df <- score_df[Cells(POD), , drop = FALSE]
all(rownames(score_df) == Cells(POD))
POD <- AddMetaData(POD, metadata = score_df)
POD$FAO_TCA_OXPHOS_score <- POD$FA_degradation +
  POD$TCA_cycle +
  POD$Oxidative_phosphorylation
FeaturePlot(POD, features = "FAO_TCA_OXPHOS_score",reduction="tsne")
score <- POD@meta.data$FAO_TCA_OXPHOS_score
POD@meta.data$FAO_TCA_OXPHOS_score_01 <-
  (score - min(score, na.rm = TRUE)) /
  (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
summary(POD$FAO_TCA_OXPHOS_score_01)
FeatureDimPlot(
  POD,
  features = "FAO_TCA_OXPHOS_score_01",
  pt.size = 2,
  reduction = "tsne",
  theme_use = "theme_blank",
  xlab = "tSNE_1",
  ylab = "tSNE_2",
  palcolor = rev(c(
    "#9E0142",
    "#D53E4F",
    "#F46D43",
    "#FDAE61",
    "#FEE08B",
    "#FFFFBF",
    "#fef9ef"
  ))
)
table(POD$celltype)
plot_df <- PODm@meta.data %>%
  dplyr::select(
    celltype,
    FA_degradation,
    TCA_cycle,
    Oxidative_phosphorylation
  ) %>%
  tidyr::pivot_longer(
    cols = c(
      FA_degradation,
      TCA_cycle,
      Oxidative_phosphorylation
    ),
    names_to = "Pathway",
    values_to = "Activity"
  )
plot_df$celltype <- factor(
  plot_df$celltype,
  levels = c("hPOD", "iPOD1", "iPOD2")
)
plot_df$Pathway <- factor(
  plot_df$Pathway,
  levels = c(
    "FA_degradation",
    "TCA_cycle",
    "Oxidative_phosphorylation"
  ),
  labels = c(
    "FA degradation",
    "TCA cycle",
    "OXPHOS"
  )
)
anova_res <- plot_df %>%
  group_by(Pathway) %>%
  anova_test(Activity ~ celltype)
anova_res
get_anova_table(anova_res)
tukey_res <- plot_df %>%
  group_by(Pathway) %>%
  tukey_hsd(Activity ~ celltype)
tukey_res
stat.test <- tukey_res %>%
  mutate(
    p_label = case_when(
      p.adj < 0.001 ~ paste0("p = ", formatC(p.adj, format = "e", digits = 2)),
      TRUE ~ paste0("p = ", signif(p.adj, 3))
    )
  )
stat.test <- stat.test %>%
  add_xy_position(
    x = "celltype",
    data = plot_df,
    formula = Activity ~ celltype,
    scales = "free_y",
    step.increase = 0.12
  )
ggplot(
  plot_df,
  aes(x = celltype, y = Activity, fill = celltype)
) +
  geom_violin(
    trim = FALSE,
    scale = "width",
    color = "black",
    linewidth = 0.4,
    alpha = 0.85
  ) +
  geom_boxplot(
    width = 0.25,
    outlier.shape = NA,
    color = "white",
    linewidth = 0.8,
    alpha = 0.9
  ) +
  geom_jitter(
    aes(color = celltype),
    width = 0.12,
    size = 0.8,
    alpha = 0.45,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ Pathway,
    nrow = 1,
    scales = "free_y"
  ) +
  stat_pvalue_manual(
    stat.test,
    label = "p_label",
    tip.length = 0.01,
    bracket.size = 0.4,
    size = 5
  ) +
  scale_fill_manual(
    values = c(
      "hPOD" = "#b2df8a",
      "iPOD1" = "#FB9A99",
      "iPOD2" = "#E31A1C"
    )
  ) +
  scale_color_manual(
    values = c(
      "hPOD" = "#b2df8a",
      "iPOD1" = "#FB9A99",
      "iPOD2" = "#E31A1C"
    )
  ) +
  labs(
    x = NULL,
    y = "Pathway activity (mouse)",
    fill = NULL
  ) +
  theme_classic(base_size = 20) +
  theme(
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.4
    ),
    strip.text = element_text(
      size = 21,
      face = "bold"
    ),
    axis.text.x = element_text(
      size = 20,
      color = "black",
      angle = 0,
      hjust = 0.5
    ),
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 20,
      face = "bold"
    ),
    legend.position = "none",
    panel.spacing = unit(1.2, "lines")
  )
--------------------------------------------------------------------------------
kegg_ids <- c(
  "Fatty_acid_degradation" = "00071",
  "Citrate_cycle_TCA_cycle" = "00020",
  "Oxidative_phosphorylation" = "00190"
)
path2eg <- as.list(org.Mm.egPATH2EG)
mouse_kegg_gene_sets <- lapply(kegg_ids, function(pid) {
  entrez_ids <- path2eg[[pid]]
  symbols <- mapIds(
    org.Mm.eg.db,
    keys = entrez_ids,
    keytype = "ENTREZID",
    column = "SYMBOL",
    multiVals = "first"
  )
  symbols <- unique(na.omit(as.character(symbols)))
  return(symbols)
})
writeGmtPathways(
  pathways = mouse_kegg_gene_sets,
  gmt.file  = "E:/OB/OB/Mice/mouse_KEGG_FA_TCA_OXPHOS.gmt"
)
PODm_daixie=sc.metabolism.Seurat.2(PODm, method = "VISION", imputation = F, ncores = 2,
                                  geneList = "E:/OB/OB/Mice/mouse_KEGG_FA_TCA_OXPHOS.gmt")
datam=PODm_daixie@assays[["METABOLISM"]][["score"]]
colnames(datam) <- gsub("\\.1_", "-1_", colnames(datam))
head(colnames(datam))
head(Cells(PODm))
length(intersect(colnames(datam), Cells(PODm)))
pathways_use <- c(
  "Fatty_acid_degradation",
  "Citrate_cycle_TCA_cycle",
  "Oxidative_phosphorylation"
)
pathways_use %in% rownames(datam)
score_mat <- datam[pathways_use, , drop = FALSE]
score_df <- as.data.frame(t(score_mat))
colnames(score_df) <- c(
  "FA_degradation",
  "TCA_cycle",
  "Oxidative_phosphorylation"
)
all(rownames(score_df) %in% Cells(PODm))
all(Cells(PODm) %in% rownames(score_df))
score_df <- score_df[Cells(PODm), , drop = FALSE]
all(rownames(score_df) == Cells(PODm))
PODm <- AddMetaData(PODm, metadata = score_df)
PODm$FAO_TCA_OXPHOS_score <- PODm$FA_degradation +
  PODm$TCA_cycle +
  PODm$Oxidative_phosphorylation
score <- PODm@meta.data$FAO_TCA_OXPHOS_score
PODm@meta.data$FAO_TCA_OXPHOS_score_01 <-
  (score - min(score, na.rm = TRUE)) /
  (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
summary(POD$FAO_TCA_OXPHOS_score_01)
FeatureDimPlot(
  PODm,
  features = "FAO_TCA_OXPHOS_score_01",
  pt.size = 2,
  reduction = "tsne",
  theme_use = "theme_blank",
  xlab = "tSNE_1",
  ylab = "tSNE_2",
  palcolor = rev(c(
    "#9E0142",
    "#D53E4F",
    "#F46D43",
    "#FDAE61",
    "#FEE08B",
    "#FFFFBF",
    "#fef9ef"
  ))
)
POD_death=sc.metabolism.Seurat.2(POD, method = "VISION", imputation = F, ncores = 2,
                                  geneList = "E:/OB///death1.gmt")
data=POD_death@assays[["METABOLISM"]][["score"]]
colnames(data) <- gsub("\\.1_", "-1_", colnames(data))
head(colnames(data))
head(Cells(POD))
length(intersect(colnames(data), Cells(POD)))
common_cells <- intersect(colnames(data), Cells(POD))
common_cells <- colnames(data)[colnames(data) %in% common_cells]
length(common_cells)
meta <- POD@meta.data[common_cells, c("celltype", "disease"), drop = FALSE]
score_df <- as.data.frame(t(data[, common_cells, drop = FALSE]))
new_data <- cbind(
  cell = rownames(score_df),
  meta,
  score_df
)
head(new_data)
head(human_gmt)
human_genes <- unique(human_gmt$gene)
orth_map <- orthologs(
  genes = human_genes,
  species = "mouse"
)
head(orth_map)
orth_map2 <- orth_map[, c("human_symbol", "symbol", "support_n")]
colnames(orth_map2) <- c("human_symbol", "mouse_symbol", "support_n")
head(orth_map2)
mouse_gmt <- merge(
  human_gmt,
  orth_map2,
  by.x = "gene",
  by.y = "human_symbol"
)
mouse_gmt <- mouse_gmt[, c("term", "mouse_symbol")]
colnames(mouse_gmt) <- c("term", "gene")
mouse_gmt <- unique(mouse_gmt)
head(mouse_gmt)
write_gmt <- function(gmt_df, file) {
  gene_sets <- split(gmt_df$gene, gmt_df$term)
  con <- file(file, open = "w")
  for (term in names(gene_sets)) {
    genes <- unique(gene_sets[[term]])
    line <- paste(c(term, "na", genes), collapse = "\t")
    writeLines(line, con)
  }
  close(con)
}
out_file <- "E:/OB///death1_mouse.gmt"
write_gmt(mouse_gmt, out_file)
PODm_death=sc.metabolism.Seurat.2(PODm, method = "VISION", imputation = F, ncores = 2,
                                   geneList = "E:/OB///death1_mouse.gmt")
datam=PODm_death@assays[["METABOLISM"]][["score"]]
colnames(datam) <- gsub("\\.1_", "-1_", colnames(datam))
head(colnames(datam))
head(Cells(PODm))
length(intersect(colnames(datam), Cells(PODm)))
common_cells <- intersect(colnames(datam), Cells(PODm))
common_cells <- colnames(datam)[colnames(datam) %in% common_cells]
length(common_cells)
metam <- PODm@meta.data[common_cells, c("celltype", "disease"), drop = FALSE]
score_dfM <- as.data.frame(t(datam[, common_cells, drop = FALSE]))
new_datam <- cbind(
  cell = rownames(score_dfM),
  metam,
  score_dfM
)
head(new_datam)
celltype_mean <- aggregate(
  new_data[, grep("Apoptosis|Pyroptosis|Ferroptosis|Autophagy|Necroptosis|Cuproptosis|Parthanatos|Entotic|Netotic|Lysosome|Alkaliptosis|Oxeiptosis|NETosis|Immunogenic|Anoikis|Paraptosis|Methuosis|Entosis", names(new_data))],
  by = list(Celltype = new_data$celltype),
  FUN = mean,
  na.rm = TRUE
)
disease_mean <- aggregate(
  new_data[, grep("Apoptosis|Pyroptosis|Ferroptosis|Autophagy|Necroptosis|Cuproptosis|Parthanatos|Entotic|Netotic|Lysosome|Alkaliptosis|Oxeiptosis|NETosis|Immunogenic|Anoikis|Paraptosis|Methuosis|Entosis", names(new_data))],
  by = list(Disease = new_data$disease),
  FUN = mean,
  na.rm = TRUE
)
celltypem_mean <- aggregate(
  new_datam[, grep("Apoptosis|Pyroptosis|Ferroptosis|Autophagy|Necroptosis|Cuproptosis|Parthanatos|Entotic|Netotic|Lysosome|Alkaliptosis|Oxeiptosis|NETosis|Immunogenic|Anoikis|Paraptosis|Methuosis|Entosis", names(new_data))],
  by = list(Celltype = new_datam$celltype),
  FUN = mean,
  na.rm = TRUE
)
diseasem_mean <- aggregate(
  new_datam[, grep("Apoptosis|Pyroptosis|Ferroptosis|Autophagy|Necroptosis|Cuproptosis|Parthanatos|Entotic|Netotic|Lysosome|Alkaliptosis|Oxeiptosis|NETosis|Immunogenic|Anoikis|Paraptosis|Methuosis|Entosis", names(new_data))],
  by = list(Disease = new_datam$disease),
  FUN = mean,
  na.rm = TRUE
)
head(disease_mean)
rownames(disease_mean)=disease_mean$Disease
disease_mat=disease_mean[, -1]
disease_mat=as.matrix(disease_mat)
color=rev(c("#D53E4F","#F46D43","#FDAE61","#FEE08B",
            "#FFFFBF","#E6F598","#ABDDA4","#66C2A5"))
jianbian <- colorRampPalette(color)
pheatmap(
  t(disease_mat),
  scale = "row",
  cluster_rows = F,
  cluster_cols = F,
  color = jianbian(100),
  border_color = "white",
  number_color = "white",
  number_format = "%.1f",
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 20,
  fontsize_number =20 ,
  fontsize_col = 20,
  display_numbers = T,
  angle_col = 90)
head(new_data)
head(new_datam)
human_df <- new_data %>%
  mutate(Species = "Human")
mouse_df <- new_datam %>%
  mutate(Species = "Mouse")
plot_df <- bind_rows(human_df, mouse_df)
plot_df$celltype <- factor(
  plot_df$celltype,
  levels = c("hPOD", "iPOD1", "iPOD2")
)
plot_df$Species <- factor(
  plot_df$Species,
  levels = c("Human", "Mouse")
)
anova_res <- plot_df %>%
  group_by(Species) %>%
  anova_test(Immunogenic ~ celltype)
anova_res
tukey_res <- plot_df %>%
  group_by(Species) %>%
  tukey_hsd(Immunogenic ~ celltype)
tukey_res
tukey_res <- tukey_res %>%
  add_xy_position(
    x = "celltype",
    scales = "free_y"
  )
