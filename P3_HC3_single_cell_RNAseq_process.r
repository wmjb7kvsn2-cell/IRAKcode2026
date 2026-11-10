####Wujunfeng####
getwd()
setwd('/home/yc/R_analysis/Wujunfeng/')
library(SCP)
data("pancreas_sub")
library('Seurat')
setRepositories(ind = 1:3, addURLs = c('https://satijalab.r-universe.dev', 'https://bnprks.r-universe.dev/'))
install.packages(c("BPCells", "presto", "glmGamPoi"))
remotes::install_github("satijalab/seurat-data", quiet = TRUE)
remotes::install_github("satijalab/azimuth", quiet = TRUE)
remotes::install_github("satijalab/seurat-wrappers", quiet = TRUE)
BiocManager::install("TFBSTools", type = "source", force = TRUE)
BiocManager::install("DirichletMultinomial", type="source")
install.packages('XML')
remotes::install_github(
  "satijalab/azimuth",
  upgrade = "never",
  dependencies = TRUE
)
remotes::install_local(
  "~/azimuth",
  dependencies = TRUE,
  upgrade = "never"
)
library(Seurat)
HC1 <- Read10X("/home/yc/R_analysis/Wujunfeng/XYC_HC/XYC_HC_S49_L1/Outs/filtered_feature_bc_matrix/")
HC1 <- CreateSeuratObject(counts = HC1, project = "HC1", min.cells = 3, min.features = 200)
HC2 <- Read10X("/home/yc/R_analysis/Wujunfeng/GSM6833297/")
HC2 <- CreateSeuratObject(counts = HC2, project = "HC2", min.cells = 3, min.features = 200)
HC3=HC2
HC2=HC1
P2<-Read10X("/home/yc/R_analysis/Wujunfeng/LYX/LYX1027_S19_L003/Outs/filtered_feature_bc_matrix/")
P2 <- CreateSeuratObject(counts = P2, project = "P2", min.cells = 3, min.features = 200)
P3<-Read10X("/home/yc/R_analysis/Wujunfeng/CYM/CYM_S32_L2/Outs/filtered_feature_bc_matrix/")
P3 <- CreateSeuratObject(counts = P3, project = "P3", min.cells = 3, min.features = 200)

#线粒体比例
HC2[["percent.mt"]] <- PercentageFeatureSet(HC2, pattern = "^mt-")
HC3[["percent.mt"]] <- PercentageFeatureSet(HC3, pattern = "^mt-")

P2[["percent.mt"]] <- PercentageFeatureSet(P2, pattern = "^mt-")
P3[["percent.mt"]] <- PercentageFeatureSet(P3, pattern = "^mt-")

#血红蛋白基因
HC2[["percent.hb"]] <- PercentageFeatureSet(HC2, pattern = "^Hb[^(p)]")
HC3[["percent.hb"]] <- PercentageFeatureSet(HC3, pattern = "^Hb[^(p)]")

P2[["percent.hb"]] <- PercentageFeatureSet(P2, pattern = "^Hb[^(p)]")
P3[["percent.hb"]] <- PercentageFeatureSet(P3, pattern = "^Hb[^(p)]")

#核糖体基因
HC2[["percent.rb"]] <- PercentageFeatureSet(HC2, pattern = "^Rbs|Rpl")
HC3[["percent.rb"]] <- PercentageFeatureSet(HC3, pattern = "^Rbs|Rpl")

P2[["percent.rb"]] <- PercentageFeatureSet(P2, pattern = "^Rbs|Rpl")
P3[["percent.rb"]] <- PercentageFeatureSet(P3, pattern = "^Rbs|Rpl")

p1= VlnPlot(HC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),pt.size = 0.1, ncol = 5)
p2 = VlnPlot(HC3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),pt.size = 0.1, ncol = 5)
p3 = VlnPlot(P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),pt.size = 0.1, ncol = 5)
p4 = VlnPlot(P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.hb","percent.rb"),pt.size = 0.1, ncol = 5)

p1/p2/p3/p4

#QC质控
HC2 <- subset(HC2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000 &  percent.mt < 5)
HC3 <- subset(HC3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000 &  percent.mt < 5)
P2 <- subset(P2, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000 &  percent.mt < 5)
P3 <- subset(P3, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 30000 &  percent.mt < 5)
#merge data
sce <- merge(HC2, y=c(HC3,P2,P3),add.cell.ids = c("HC2","HC3", "P2","P3"))
#NormalizeData & ScaleData
sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce)
sce <- ScaleData(sce, vars.to.regress = c("percent.mt"))
sce <- RunPCA(sce, verbose=F)
sce_cca <- IntegrateLayers(object = sce, 
                           method = HarmonyIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.harmony",verbose = FALSE)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(Azimuth)
library(ggplot2)
library(patchwork)
sce_cca[["RNA"]] <- JoinLayers(sce_cca[["RNA"]])
Seurat::ElbowPlot(sce_cca, ndims = 50)
sce_cca <- FindNeighbors(sce_cca, reduction = "integrated.harmony", dims = 1:30)
sce_cca <- FindClusters(sce_cca, resolution = 1.0)
sce_cca <- RunUMAP(sce_cca, dims = 1:30, reduction = "integrated.harmony")
library(ggsci)
getwd()
setwd('/home/yc/R_analysis/Wujunfeng/')
pdf(file = './harmony.merged.res1.umap.pdf',width = 6,height = 5)
DimPlot(sce_cca, reduction = "umap", group.by = "orig.ident")+
  ggtitle("Harmony")+scale_color_bmj()
dev.off()
library(RColorBrewer)

cols <- c(brewer.pal(12,"Paired"),brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),
          brewer.pal(12,"Set3"),brewer.pal(8,"Pastel2"),brewer.pal(9,"Pastel1"),brewer.pal(8,"Accent"))

pdf(file = './harmony.merged.res1.umap.pdf',width = 6,height = 5)
DimPlot(sce_cca, reduction = "umap",cols = cols)+
  ggtitle("Harmony")
dev.off()
Allmarkers <- FindAllMarkers(sce_cca, logfc.threshold = 0.3, min.pct = 0.3, only.pos = T)
write.csv(Allmarkers, file = 'Allmarkers.res1.dim30.csv')
BiocManager::install('SingleR')
install.packages('SingleR')
install_github("dviraran/SingleR")
BiocManager::install("SingleR" ,version = "3.20", force = TRUE)
devtools::install_github("chloelulu/scMayoMap")
library(scMayoMap)
obj<-scMayoMap(data = Allmarkers,tissue = 'blood')
res<-obj$res
head(res)
markers<-obj$markers
head(markers)
plt<-scMayoMap.plot(scMayoMap.object = obj,directory = './',width = 12,height = 8)
plt
sce_cca=subset(sce_cca,seurat_clusters!=12)
sce_cca=subset(sce_cca,seurat_clusters!=14)
cell_map <- c(
  `0`="Neutrophil", `1`="Neutrophil", `24`="Neutrophil",
  `2`="CD4_Naive", `9`="CD4_Naive", `15`="CD4_Naive",
  `3`="B_Naive", `8`="gdT", `18`="B_Naive", `21`="B_Naive",
  `6`="CD14_mono", `16`="CD14_mono",
  `23`="DC", `26`="DC",
  `11`="CD8_Naive",
  `4`="CD8_Temra",
  `5`="CD8_Tcm",
  `7`="NKT",
  `10`="CD4_Tcm",
  `13`="B_helper",
  `17`="CD16_mono",
  `19`="Megakaryocyte",
  `20`="NK",
  `22`="Plasma",
  `25`="Basophils"
)

sce_cca$cell_type <- unname(
  cell_map[as.character(sce_cca$seurat_clusters)]
)

sce_cca$cell_type[is.na(sce_cca$cell_type)] <- "Proliferating"
table(sce_cca$seurat_clusters, sce_cca$cell_type)                         
Idents(sce_cca)=sce_cca$cell_type
save(sce_cca,file = 'pbmc.wujunfeng.res1.renamed.20260928.rds')
library(SCP)
dim.use=1:30
pdf("./CellCluster-DimPlot_umap_res1..renamed.label.pdf"),width = 6,height = 5)
CellDimPlot(sce_cca, group.by = "cell_type", reduction = "UMAP", theme_use = "theme_blank", label = FALSE)  +
  guides(color = guide_legend(override.aes = list(size = 2)))  # 让legend中的点变大
dev.off()
pdf("./CellCluster-DimPlot_umap_res1.1.renamed.label.pdf",width = 6,height = 5)
CellDimPlot(srt = sce_cca, group.by = "cell_type", reduction = "UMAP", label = TRUE)
dev.off()
sce_cca$sample=ifelse(sce_cca$orig.ident%in%c('HC2','HC3'),'HC','Patient')
sce_cca$sample=factor(sce_cca$sample,levels = c('HC','Patient'))
pdf(file = 'sce_cca.remove.counts.sample.pdf',width = 6,height = 3)
CellStatPlot(sce_cca, stat.by = "sample", group.by = "cell_type", bg.by = "cell_type", palette = "Paired", stat_type = "count", position = "dodge")
dev.off()
table(sce_cca$orig.ident)
load('./pbmc.wujunfeng.res1.renamed.20260928.rds')
CellStatPlot(sce_cca, stat.by = "cell_type", group.by = "sample", plot_type = "trend")
sce_cca$orig.ident=ifelse(sce_cca$orig.ident=='HC1','HC2',
                          ifelse(sce_cca$orig.ident=='HC2','HC3',
                                 ifelse(sce_cca$orig.ident=='P2','P2','P3')))
sce_cca$orig.ident=factor(sce_cca$orig.ident,levels = c('HC2','P2','HC3','P3'))
CellStatPlot(sce_cca, stat.by = "cell_type", group.by = "orig.ident", plot_type = "trend")
library(Seurat)
library(msigdbr)
library(UCell)
library(dplyr)
library(ggplot2)
library(patchwork)

hallmark <- msigdbr(
  species = "Homo sapiens",
  collection = "H"
)
library(GSEABase)
gmt <- getGmt("/home/yc/h.all.v2026.1.Hs.symbols.gmt")
gene_sets <- geneIds(gmt)
head(names(gene_sets))
length(gene_sets)
pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)

gene_sets4 <- gene_sets[pathways]
sapply(gene_sets4, length)
library(UCell)

sce_cca <- AddModuleScore_UCell(
  sce_cca,
  features = gene_sets4,
  name = NULL
)
grep("HALLMARK", colnames(sce_cca@meta.data), value = TRUE)
table(sce_cca$orig.ident)
plot_split_violin <- function(
    object,
    score,
    sample1,
    sample2,
    title = NULL,
    color1 = "#2878B5",
    color2 = "#D73027"
) {
  
  df <- object@meta.data %>%
    filter(sample %in% c(sample1, sample2)) %>%
    filter(!is.na(cell_type)) %>%
    mutate(sample = factor(sample, levels = c(sample1, sample2)))
  
  ggplot(df, aes(x = cell_type, y = .data[[score]])) +
    
    geom_half_violin(
      data = df %>% filter(sample == sample1),
      aes(fill = sample),
      side = "l",
      trim = TRUE,
      scale = "width",
      color = NA
    ) +
    
    geom_half_violin(
      data = df %>% filter(sample == sample2),
      aes(fill = sample),
      side = "r",
      trim = TRUE,
      scale = "width",
      color = NA
    ) +
    
    scale_fill_manual(
      values = setNames(
        c(color1, color2),
        c(sample1, sample2)
      )
    ) +
    
    labs(
      x = NULL,
      y = "UCell score",
      title = title
    ) +
    
    theme_classic(base_size = 11) +
    
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = 9
      ),
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(
        size = 11,
        face = "bold"
      )
    )
}
grep(
  "HALLMARK",
  colnames(sce_cca@meta.data),
  value = TRUE
)
score1 <- "HALLMARK_INFLAMMATORY_RESPONSE"
score2 <- "HALLMARK_IL6_JAK_STAT3_SIGNALING"
score3 <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
score4 <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"
VlnPlot(
  sce_cca,
  features = score1,
  group.by = "cell_type",
  split.by = "orig.ident",
  idents = NULL,
  pt.size = 0
)
obj_HC2P2 <- subset(
  sce_cca,
  subset = orig.ident %in% c("HC2", "P2")
)

p1 <- VlnPlot(
  obj_HC2P2,
  features = score1,
  group.by = "cell_type",
  split.by = "orig.ident",
  pt.size = 0,
  cols = c("#2878B5", "#D73027")
) +
  ggtitle("Hallmark inflammatory response") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p1
plot_half_violin <- function(
    object,
    score,
    sample1,
    sample2,
    title,
    col1 = "#2878B5",
    col2 = "#D73027",
    show_x = FALSE,
    show_legend = FALSE
) {
  
  df <- object@meta.data %>%
    filter(orig.ident %in% c(sample1, sample2)) %>%
    filter(!is.na(cell_type)) %>%
    mutate(
      orig.ident = factor(
        orig.ident,
        levels = c(sample1, sample2)
      )
    )
  
  p <- ggplot(
    df,
    aes(
      x = cell_type,
      y = .data[[score]]
    )
  ) +
    
    # 左半边：HC
    geom_half_violin(
      data = df %>% filter(orig.ident == sample1),
      aes(fill = orig.ident),
      side = "l",
      scale = "width",
      trim = TRUE,
      width = 0.9,
      color = "grey25",
      linewidth = 0.25
    ) +
    
    # 右半边：P
    geom_half_violin(
      data = df %>% filter(orig.ident == sample2),
      aes(fill = orig.ident),
      side = "r",
      scale = "width",
      trim = TRUE,
      width = 0.9,
      color = "grey25",
      linewidth = 0.25
    ) +
    
    scale_fill_manual(
      values = setNames(
        c(col1, col2),
        c(sample1, sample2)
      )
    ) +
    
    labs(
      x = NULL,
      y = "UCell Score",
      title = title,
      fill = NULL
    ) +
    
    theme_classic(base_size = 11) +
    
    theme(
      plot.title = element_text(
        size = 11,
        face = "plain",
        hjust = 0,
        margin = margin(b = 2)
      ),
      
      axis.title.y = element_text(size = 9),
      
      axis.text.y = element_text(size = 8),
      
      axis.line = element_line(
        linewidth = 0.4
      ),
      
      axis.ticks = element_line(
        linewidth = 0.4
      ),
      
      plot.margin = margin(
        t = 1,
        r = 2,
        b = 1,
        l = 2
      ),
      
      legend.position = if (show_legend) "top" else "none",
      
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 11)
    )
  
  # 只有最下面一幅图显示 cell type
  if (!show_x) {
    
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    
  } else {
    
    p <- p +
      theme(
        axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          vjust = 1,
          size = 8
        )
      )
  }
  
  return(p)
}
library(gghalves)
p1 <- plot_half_violin(
  sce_cca,
  score1,
  "HC2", "P2",
  "Hallmark_Inflammatory_response",
  col1 = "#2878B5",
  col2 = "#E53935",
  show_x = FALSE,
  show_legend = TRUE
)

p2 <- plot_half_violin(
  sce_cca,
  score2,
  "HC2", "P2",
  "Hallmark_IL6_JAK_STAT3_signaling",
  col1 = "#2878B5",
  col2 = "#E53935",
  show_x = FALSE
)

p3 <- plot_half_violin(
  sce_cca,
  score3,
  "HC2", "P2",
  "Hallmark_Interferon_alpha_response",
  col1 = "#2878B5",
  col2 = "#E53935",
  show_x = FALSE
)

p4 <- plot_half_violin(
  sce_cca,
  score4,
  "HC2", "P2",
  "Hallmark_Interferon_gamma_response",
  col1 = "#2878B5",
  col2 = "#E53935",
  show_x = TRUE
)
p_HC2P2 <- p1 / p2 / p3 / p4 +
  plot_layout(heights = c(1, 1, 1, 1))

p_HC2P2

# HC3: Light cyan
# P3: Deep orange-brown

q1 <- plot_half_violin(
  sce_cca, score1,
  "HC3", "P3",
  "Hallmark_Inflammatory_response",
  col1 = "#56B4E9",
  col2 = "#D55E00",
  show_x = FALSE,
  show_legend = TRUE
)

q2 <- plot_half_violin(
  sce_cca, score2,
  "HC3", "P3",
  "Hallmark_IL6_JAK_STAT3_signaling",
  col1 = "#56B4E9",
  col2 = "#D55E00",
  show_x = FALSE,
  show_legend = FALSE
)

q3 <- plot_half_violin(
  sce_cca, score3,
  "HC3", "P3",
  "Hallmark_Interferon_alpha_response",
  col1 = "#56B4E9",
  col2 = "#D55E00",
  show_x = FALSE,
  show_legend = FALSE
)

q4 <- plot_half_violin(
  sce_cca, score4,
  "HC3", "P3",
  "Hallmark_Interferon_gamma_response",
  col1 = "#56B4E9",
  col2 = "#D55E00",
  show_x = TRUE,
  show_legend = FALSE
)

# Combine plots
p_HC3P3 <- q1 / q2 / q3 / q4

p_HC3P3

HC3 = "#00A087"
P3  = "#3C5488"
HC3 = "#4DBBD5"
P3  = "#E64B35"
save(sce_cca,file = 'sce_acc.res1.renamed.hallmark.score.rds')
