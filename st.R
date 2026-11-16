library(Seurat)
library(SCNT2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dbscan)
library(scMetabolism)
segmented_objects <- list(Con_HD, ORG1_2, ORG3)
segmented_objects <- lapply(segmented_objects, function(object) {
  object <- JoinLayers(object)
  object <- NormalizeData(object)
  object
})
Con_HD <- segmented_objects[[1]]
ORG1_2 <- segmented_objects[[2]]
ORG3 <- segmented_objects[[3]]
Con_HD <- easyRCTD(
  object = Con_HD,
  ref = reference,
  mode = "seg",
  max_cores = 8,
  UMI_min = 20,
  doublet_mode = "doublet"
)
ORG1_2 <- easyRCTD(
  object = ORG1_2,
  ref = reference,
  mode = "seg",
  max_cores = 8,
  UMI_min = 20,
  doublet_mode = "doublet"
)
ORG3 <- easyRCTD(
  object = ORG3,
  ref = reference,
  mode = "seg",
  max_cores = 8,
  UMI_min = 20,
  doublet_mode = "doublet"
)
assign_celltype <- function(object) {
  object$celltype <- ifelse(
    object$RCTD_spot_class %in% c("doublet_certain", "doublet_uncertain"),
    as.character(object$RCTD_second_type),
    as.character(object$RCTD_first_type)
  )
  subset(object, subset = RCTD_spot_class != "reject" & !is.na(celltype))
}
Con_HD <- assign_celltype(Con_HD)
ORG1_2 <- assign_celltype(ORG1_2)
ORG3 <- assign_celltype(ORG3)
spatial_colors <- c("#ead2bd", "#a4cde1", "#4DAF4A", "#F781BF", "#BC9DCC", "#A65628", "#f9b769", "#e32427", "#D2EBC8", "#d9d9d9", "#BBDD78")
stPlot(Con_HD, group = "celltype", color_list = spatial_colors, point_size = 0.5, scale = "hires", image_scale = "slice1", show_image = FALSE, draw = "seg", polygon.border = "NA", polygon.alpha = 0.7, theme = 1)
stPlot(ORG1_2, group = "celltype", color_list = spatial_colors, point_size = 0.5, scale = "hires", image_scale = "slice1", show_image = FALSE, draw = "seg", polygon.border = "NA", polygon.alpha = 0.7, theme = 1)
stPlot(ORG3, group = "celltype", color_list = spatial_colors, point_size = 0.5, scale = "hires", image_scale = "slice1", show_image = FALSE, draw = "seg", polygon.border = "NA", polygon.alpha = 0.7, theme = 1)
stFeature(Con_HD, feature = "NPHS2", assay = NULL, layer = "data", color = c("grey95", "gold", "red"), scale = "hires", image_scale = "slice1", show_image = FALSE, draw = "seg", polygon.border = "NA", polygon.alpha = 0.7, theme = 1)
stFeature(ORG1_2, feature = "NPHS2", assay = NULL, layer = "data", color = c("grey95", "gold", "red"), scale = "hires", image_scale = "slice1", show_image = FALSE, draw = "seg", polygon.border = "NA", polygon.alpha = 0.7, theme = 1)
stFeature(ORG3, feature = "NPHS2", assay = NULL, layer = "data", color = c("grey95", "gold", "red"), scale = "hires", image_scale = "slice1", show_image = FALSE, draw = "seg", polygon.border = "NA", polygon.alpha = 0.7, theme = 1)
Con_HD$group <- "0h"
ORG1_2$group <- "18h"
ORG3$group <- "36h"
POD_HD <- subset(merge(Con_HD, y = list(ORG1_2, ORG3)), subset = celltype == "POD")
POD_HD <- JoinLayers(POD_HD)
POD_HD <- NormalizeData(POD_HD)
POD_metabolism <- sc.metabolism.Seurat(POD_HD, method = "VISION", imputation = FALSE, ncores = 2, geneList = "KEGG")
metabolism_score <- POD_metabolism[["METABOLISM"]][["score"]]
metabolism_score <- metabolism_score[, Cells(POD_HD), drop = FALSE]
pathways_use <- c("Fatty acid degradation", "Citrate cycle (TCA cycle)", "Oxidative phosphorylation")
pathways_use <- intersect(pathways_use, rownames(metabolism_score))
pathway_df <- as.data.frame(t(metabolism_score[pathways_use, , drop = FALSE]))
colnames(pathway_df) <- c("FA_degradation", "TCA_cycle", "Oxidative_phosphorylation")
POD_HD <- AddMetaData(POD_HD, pathway_df[Cells(POD_HD), , drop = FALSE])
POD_HD$FAO_TCA_OXPHOS_score <- rowSums(POD_HD@meta.data[, c("FA_degradation", "TCA_cycle", "Oxidative_phosphorylation"), drop = FALSE], na.rm = TRUE)
plot_df <- POD_HD@meta.data %>%
  select(group, FA_degradation, TCA_cycle, Oxidative_phosphorylation) %>%
  pivot_longer(-group, names_to = "Pathway", values_to = "Activity")
plot_df$group <- factor(plot_df$group, levels = c("0h", "18h", "36h"))
ggplot(plot_df, aes(x = group, y = Activity, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", color = "black") +
  geom_boxplot(width = 0.15, outlier.shape = NA, color = "white") +
  geom_jitter(aes(color = group), width = 0.12, size = 0.5, alpha = 0.35, show.legend = FALSE) +
  facet_wrap(~ Pathway, scales = "free_y") +
  scale_fill_manual(values = c("0h" = "#b2df8a", "18h" = "#FB9A99", "36h" = "#E31A1C")) +
  scale_color_manual(values = c("0h" = "#b2df8a", "18h" = "#FB9A99", "36h" = "#E31A1C")) +
  theme_classic()
neighborhood_0h <- stNeighborhood(Con_HD, group.by = "celltype", target = "POD", neighbors = c("MAC", "T"), radius_um = c(20, 50, 100), coord_unit = "pixel", microns_per_pixel = 0.24631478859116646)
neighborhood_18h <- stNeighborhood(ORG1_2, group.by = "celltype", target = "POD", neighbors = c("MAC", "T"), radius_um = c(20, 50, 100), coord_unit = "pixel", microns_per_pixel = 0.24615418764136793)
neighborhood_36h <- stNeighborhood(ORG3, group.by = "celltype", target = "POD", neighbors = c("MAC", "T"), radius_um = c(20, 50, 100), coord_unit = "pixel", microns_per_pixel = 0.24616090770735294)
neighborhood_df <- bind_rows(
  mutate(neighborhood_0h, group = "0h"),
  mutate(neighborhood_18h, group = "18h"),
  mutate(neighborhood_36h, group = "36h")
)
neighborhood_long <- neighborhood_df %>%
  pivot_longer(cols = matches("^(MAC|T)_"), names_to = "Cell_distance", values_to = "Number") %>%
  separate(Cell_distance, into = c("Cell", "Distance"), sep = "_", extra = "merge")
ggplot(neighborhood_long, aes(x = Cell, y = Number, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  facet_wrap(~ Distance, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("0h" = "#b2df8a", "18h" = "#FB9A99", "36h" = "#E31A1C")) +
  theme_classic()
