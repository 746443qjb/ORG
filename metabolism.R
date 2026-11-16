library(tidyverse)
library(readr)
library(FactoMineR)
library(factoextra)
library(ggrepel)
library(ggfun)
library(grid)
rownames(data) <- data[[1]]
data <- data[ , -1]
data_t <- t(data)
data_t <- as.data.frame(data_t)
group <- c(rep("0h", 6), rep("18h", 6), rep("36h", 6))
data_t$Group <- group
pca_res <- PCA(data_t[, -ncol(data_t)], graph = FALSE)
explained_var <- pca_res$eig[1:2, 2]
pca_data <- as.data.frame(pca_res$ind$coord)
pca_data$Group <- data_t$Group
data1 <- pca_data[, c("Dim.1", "Dim.2", "Group")]
colnames(data1) <- c("PC1", "PC2", "Group")
pca <- data1
p1=ggplot(data = pca, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(color = Group)) +
  scale_color_manual(values = c("#b2df8a", "#FB9A99", "#E31A1C")) +
  labs(color = "Group")+
  stat_ellipse(aes(color = Group),fill = NA, geom = 'polygon')+scale_shape_manual(values=c(16,15)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p2=ggplot(data = pca) +
geom_density(aes(x = PC1, fill=Group),
color = 'black', alpha = 0.66, position = 'identity') +
scale_linetype_manual(values = c('solid','dashed'))+
  scale_fill_manual(values=c("#b2df8a", "#FB9A99", "#E31A1C")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
P3=ggplot(data = pca) +
  geom_density(aes(x = PC1, fill=Group),
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = c('solid','dashed'))+
  scale_fill_manual(values=c("#b2df8a", "#FB9A99", "#E31A1C")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  coord_flip()
pp <- p2 + plot_spacer() + p1 + P3 +
  plot_layout(guides = 'collect') &
  theme(legend.position='right')
pp
pp + plot_layout(widths = c(3, 1), heights = c(1,3))
-----------------------------------------
  X <- data_t[, -ncol(data_t)]
X <- as.matrix(X)
if (!all(sapply(X, is.numeric))) {
  stop("X contains non-numeric data. Please check your data.")
}
Y <- as.factor(data_t$Group)
if (any(is.na(X)) || any(is.na(Y))) {
  stop("Data contains NA values. Please handle missing data before proceeding.")
}
plsda_res <- plsda(X, Y, ncomp = 2)
plsda_scores <- plsda_res$variates$X
explained_var <- plsda_res$prop_expl_var$X[1:2] * 100
plsda_data <- data.frame(PC1 = plsda_scores[, 1],
                         PC2 = plsda_scores[, 2],
                         Group = Y)
p1=ggplot(data = plsda_data, aes(x = PC1, y = PC2)) +
  geom_point(size = 4, aes(color = Group)) +
  scale_color_manual(values = c("#b2df8a", "#FB9A99", "#E31A1C")) +
  labs(color = "Group")+
  stat_ellipse(aes(color = Group),fill = NA, geom = 'polygon')+scale_shape_manual(values=c(16,15)) +
  theme_bw() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank())
p2=ggplot(data = plsda_data) +
  geom_density(aes(x = PC1, fill=Group),
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = c('solid','dashed'))+
  scale_fill_manual(values=c("#b2df8a", "#FB9A99", "#E31A1C")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
P3=ggplot(data = plsda_data) +
  geom_density(aes(x = PC1, fill=Group),
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = c('solid','dashed'))+
  scale_fill_manual(values=c("#b2df8a", "#FB9A99", "#E31A1C")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  coord_flip()
pp <- p2 + plot_spacer() + p1 + P3 +
  plot_layout(guides = 'collect') &
  theme(legend.position='right')
pp
pp + plot_layout(widths = c(3, 1), heights = c(1,3))
------------------------------------------
X <- as.matrix(data_t[, -ncol(data_t)])
Y <- as.factor(data_t$Group)
oplsda_res <- opls(X, Y, predI = 1, orthoI = 1)
----------------------------------------- limma
str(data)
group <- factor(c(rep("0h", 6), rep("18h", 6), rep("36h", 6)))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
fit <- lmFit(data, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
results_0h_vs_18h <- topTable(fit2, coef = "0h_vs_18h", number = Inf)
results_18h_vs_36h <- topTable(fit2, coef = "18h_vs_36h", number = Inf)
results_0h_vs_36h <- topTable(fit2, coef = "0h_vs_36h", number = Inf)
DotPlot(PODm,c("Tpi1"),assay="RNA")
DotPlot(POD,c("TPI1","GK"),assay="RNA")
--------------------------------------------
df1 <- df %>%
  mutate(change = case_when(
    VIP > 1 & padj < 0.05 & log2FoldChange > 1 ~ "Up",
    VIP > 1 & padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Normal"
  ))
p1=ggplot(data = df1) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj),
                 color = log2FoldChange,
                 size = -log10(padj))) +
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) +
  scale_size(range = c(1,7)) +
  xlim(c(-8, 8)) +
  ylim(c(-1, 14)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_blank(),
        plot.subtitle = element_blank()
  )
df2 <- df2 %>%
  mutate(change = case_when(
    VIP > 1 & padj < 0.05 & log2FoldChange > 1 ~ "Up",
    VIP > 1 & padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE ~ "Normal"
  ))
p2=ggplot(data = df2) +
  geom_point(aes(x = log2FoldChange, y = -log10(padj),
                 color = log2FoldChange,
                 size = -log10(padj))) +
  scale_color_gradientn(colours = c("#3288bd", "#66c2a5","#ffffbf", "#f46d43", "#9e0142"),
                        values = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = 4) +
  scale_size(range = c(1,7)) +
  xlim(c(-9, 9)) +
  ylim(c(-1, 14)) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        axis.text = element_text(size = 13, color = "#000000"),
        axis.title = element_text(size = 15),
        plot.title = element_blank(),
        plot.subtitle = element_blank()
  )
grid.newpage()
combined_plot <- grid.arrange(p1, p2, ncol = 2)
--------------------------------------------------
df_list <- as.list(df)
df_list <- lapply(df_list, function(x) Filter(nzchar, x))
ggVennDiagram(df_list,edge_size=1,set_size = 6,size=5,label_alpha=0,) + scale_fill_gradient(low="white",high = "red")
