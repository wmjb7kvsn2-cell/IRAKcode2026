#20200930  测试用reselution=1.5

rm(list=ls())

##载入seurat???
library(dplyr)
library(Seurat)

##读入pbmc数据
sp1 <- Read10X(data.dir = "YJZHC/filtered_feature_bc_matrix/")
youngliver <- CreateSeuratObject(counts = sp1, project = "HC", min.cells = 3, min.features = 200)
youngliver

##读入pbmc1数据
sp2 <- Read10X(data.dir = "YJZ/filtered_feature_bc_matrix/")
oldliver <- CreateSeuratObject(counts = sp2, project = "PT", min.cells = 3, min.features = 200)
oldliver

merged<-merge(youngliver,oldliver)
table(merged$orig.ident)

merged[["percent.mt"]] <- PercentageFeatureSet(merged, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB.genes <- CaseMatch(HB.genes, rownames(merged))
merged[["percent.HB"]]<-PercentageFeatureSet(merged, features=HB.genes) 
merged[["percent.rb"]] <- PercentageFeatureSet(merged, pattern = "^RP[LS]")
VlnPlot(merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.HB"), ncol = 3)

plot1 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

str(merged)

merged <- subset(merged, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 10 & percent.HB < 1)

table(merged$orig.ident)

##A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 )
merged <- NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000)
merged

merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 3000)

##细胞周期归类
g2m_genes <- cc.genes$g2m.genes
g2m_genes <- CaseMatch(search=g2m_genes, match=rownames(merged))
s_genes <- cc.genes$s.genes    
s_genes <- CaseMatch(search=s_genes, match=rownames(merged))
merged <- CellCycleScoring(merged, g2m.features=g2m_genes, s.features=s_genes)
DimPlot(merged, reduction = "tsne",group.by = "Phase",pt.size = 1.5)

##而对所有基因进行标准化的方法如下：
all.genes <- rownames(merged)
merged <- ScaleData(merged, features = all.genes, vars.to.regress = c("percent.mt","percent.rb","S.Score","G2M.Score"))

merged <- RunPCA(merged, features = all.genes)

DimPlot(merged, reduction = "pca",group.by = "orig.ident",pt.size = 1.5)

ElbowPlot(merged, ndims = 50)

merged <- FindNeighbors(merged, dims = 1:30)
merged <- FindClusters(merged, resolution = 0.3)

merged <- RunTSNE(merged, dims = 1:40)#耗时
merged <- RunUMAP(merged, dims = 1:40)

DimPlot(merged,reduction = "tsne",label = TRUE,pt.size = 1.5)
DimPlot(merged,reduction = "umap",label = TRUE,pt.size = 1.5)

plot1<-DimPlot(subset(merged, subset = orig.ident=='HC'),reduction = "tsne",label = TRUE,pt.size = 1.5)
plot2<-DimPlot(subset(merged, subset = orig.ident=='PT'),reduction = "tsne",label = TRUE,pt.size = 1.5)
CombinePlots(plots = list(plot1, plot2))

DimPlot(merged,reduction = "tsne",label = TRUE,group.by="orig.ident",pt.size = 1.5)
DimPlot(merged,reduction = "umap",label = TRUE,group.by="celltype.main",pt.size = 1)


DefaultAssay(merged) <- "RNA"
options(future.globals.maxSize = 20 * 1024^3)

Idents(merged) <- "celltype.main"
table(Idents(merged))
ClusterMarker <- FindAllMarkers(merged, assay = "RNA", slot = "data", only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
ClusterMarker <- ClusterMarker[,c(7,1:6)]
#dir.create("process_data")
write.csv(ClusterMarker,'process_data/CelltypeMarkers.csv', row.names=F)

###### 提取差异显著的marker genes
top = 30   #可根据需要调整
TopMarkers1 <- ClusterMarker %>% filter(p_val_adj == 0) %>% group_by(cluster) %>% 
  top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>%
  top_n(n = top, wt = avg_log2FC)
TopMarkers <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(TopMarkers,'process_data/TopMarkers.csv', row.names=F)

###### 提取没有核糖体的Markers
ClusterMarker_filter <- ClusterMarker[!grepl("^RP[SL]", ClusterMarker$gene, ignore.case = F),]
ClusterMarker_filter <- ClusterMarker_filter[!grepl("^MT-", ClusterMarker_filter$gene, ignore.case = F),]
top = 30   #可根据需要调整
TopMarkers1 <- ClusterMarker_filter %>% filter(p_val_adj == 0) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
TopMarkers2 <- ClusterMarker_filter %>% filter(p_val_adj < 0.01) %>% 
  group_by(cluster) %>% top_n(n = top, wt = avg_log2FC)
ClusterMarker_filter <- rbind(TopMarkers1, TopMarkers2) %>% unique.matrix() %>% arrange(cluster)
write.csv(ClusterMarker_filter,'process_data/TopMarker_noRibo_noMito.csv', row.names=F)

#导出一个方便对比的表格(每行一个cluster的marker gene)
TopMarkers2Lines(ClusterMarker_filter, output = "process_data/")

#细胞群进行富集分析--------------------------------------
library(readxl)
library(DOSE)
options(connectionObserver = NULL)
#BiocManager::install("org.Dr.eg.db")
#library(org.Dr.eg.db)  ###瀵瑰簲鍚庨潰閮介渶瑕佹敼
#library(org.Mm.eg.db)  ###
library(org.Hs.eg.db)  ###
#BiocManager::install("topGO")
library(topGO)
library(clusterProfiler)
#BiocManager::install("pathview")
library(pathview)
#keytypes(org.Mm.eg.db) 
library(ggplot2)
library(openxlsx)

data <- ClusterMarker_filter
head(data)
type <- names(table(data$cluster))
type

# 3. Annotated and Enrichment
result1 <- "Celltype-Enrichment2"
if ( !dir.exists(result1)) {
  print("0. Create specified directory...")
  dir.create(result1)
}

go <- list()
kegg <- list()

for (i in type) {
  df <- data[data$cluster == i & data$p_val < 0.05,]
  name <- i
  df$Gene_symbol <- df$gene
  df2 <- df
  df2$avg_log2FC <- 2^df2$avg_log2FC
  colnames(df2)[2] <- "Fold_change"
  
  ## 02. Enrichment 
  result2 <- paste(result1,"/",name,sep = "")
  if ( !dir.exists(result2)) {
    print("0. Create specified directory...")
    dir.create(result2)
  }
  genelits = bitr(df$gene, 
                  fromType="SYMBOL", 
                  toType="ENTREZID", 
                  OrgDb="org.Hs.eg.db")
  write.table(genelits,paste(result2,"/01.",name,"_genelist.txt",sep = ""),sep = "\t",row.names = FALSE,quote = FALSE)
  id <- as.vector(genelits[,2])
  
  ## GO
  ego.all <- enrichGO(
    gene          = id,
    keyType = "ENTREZID",
    OrgDb         = org.Hs.eg.db,
    ont           = "ALL",    # optional : "CC" "BP" "MF" "ALL"
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,
    qvalueCutoff  = 1,
    readable      = TRUE)
  go[[name]] <- as.data.frame(ego.all)
  barplot(ego.all, showCategory=20,title="EnrichmentGO_ALL")
  write.table(ego.all,paste(result2,"/02.",name,"_GO.txt",sep = ""),row.names = FALSE,sep = "\t")
  
  ## KEGG
  ekk.all <- enrichKEGG(gene         = id,
                        organism     = "hsa",###dre=fish,mmu=mouse,hsa=human
                        pvalueCutoff = 1)
  y <- setReadable(ekk.all, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  kegg[[name]] <- as.data.frame(y)
  write.table(ekk.all, paste(result2,"/03.",name,"_KEGG.txt",sep = ""),row.names = FALSE,sep="\t")
  
  ## 03. GO & KEGG Plot 
  # go
  df1 <- as.data.frame(ego.all)
  top1 <- head(df1[,c(3,6)],15)
  top1[,2] <- -log10(top1[,2])
  p1 <- ggplot(top1,aes(reorder(Description,pvalue),pvalue))+
    geom_bar(stat="identity",position="dodge",fill="#B8002D",width = 0.5)+
    theme(axis.ticks.length=unit(0.2,'cm'),
          axis.text.x=element_text(size=10,color = "black"),
          axis.text.y=element_text(size=10,color = "black"))+
    guides(fill=guide_legend(title=NULL))+
    ggtitle("Top 15 enriched GO terms")+ coord_flip() +
    theme(axis.title = element_blank(),legend.position='none') +
    theme(panel.grid =element_blank()) +   
    theme(axis.line = element_line(size=0.5, colour = "black")) +
    theme(strip.background=element_blank(),legend.position="none",panel.background = element_rect(fill='transparent'))+ 
    theme(axis.text.x = element_text(size = 10, color = "black", face = "plain")) + 
    theme(axis.text.y = element_text(size = 10, color = "black", face = "plain"))
  pdf(paste(result2,"/04.",name,"_GO_plot.pdf",sep = ""),width = 7,height = 5)
  print(p1)
  dev.off()
  
  # kegg
  df2 <- as.data.frame(ekk.all)
  top2 <- head(df2[,c(2,5)],15)
  top2[,2] <- -log10(top2[,2])
  p2 <- ggplot(top2,aes(reorder(Description,pvalue),pvalue))+
    geom_bar(stat="identity",position="dodge",fill="#B8002D",width = 0.5)+
    theme(axis.ticks.length=unit(0.2,'cm'),
          axis.text.x=element_text(size=10,color = "black"),
          axis.text.y=element_text(size=10,color = "black"))+
    guides(fill=guide_legend(title=NULL))+
    ggtitle("Top 15 enriched KEGG Pathways")+ coord_flip() +
    theme(axis.title = element_blank(),legend.position='none') +
    theme(panel.grid =element_blank()) +   
    theme(axis.line = element_line(size=0.5, colour = "black")) +
    theme(strip.background=element_blank(),legend.position="none",panel.background = element_rect(fill='transparent'))+ 
    theme(axis.text.x = element_text(size = 10,  color = "black", face = "plain")) + 
    theme(axis.text.y = element_text(size = 10,  color = "black", face = "plain"))
  
  pdf(paste(result2,"/05.",name,"_KEGG_plot.pdf",sep = ""),width = 7,height =5)
  print(p2)
  dev.off()
}

write.xlsx(go,"Celltype-Enrichment2/GO_result.xlsx") 
write.xlsx(kegg,"Celltype-Enrichment2/KEGG_result.xlsx") 

table(merged$celltype.main)
marker<-FindMarkers(merged, group.by="celltype.main",ident.1 = "CD4 TCM", ident.2 = "CD4 Naive", logfc.threshold = 0.25, min.pct = 0.25)
write.xlsx(marker,'Celltype-Enrichment/CD4 TCM vs CD4 Naive.xlsx',quote=F,row.names =T)

#细胞鉴定-----------------------------
library(SingleR)
library(future)
library(tidyverse)
library(patchwork)

load("D:\\单细胞练习\\07 生信会客厅\\20210915-1500\\m09_course\\Resource\\SingleR_ref\\ref_BE_259s.RData")#读取软件自带参考数据
source("D:\\单细胞练习\\07 生信会客厅\\20210915-1500\\m09_course\\Resource\\sc_function.R")
merged <- cell_identify(merged, ref.BE)
p1 <- DimPlot(merged, label = T, label.size = 5,reduction = "tsne") + NoLegend()
p2 <- DimPlot(merged, group.by = "SingleR", label = T, label.size = 5,reduction = "tsne")
p2
p <- p1|p2
p
ggsave("process_data/Celltype_SingleR.png", p, width = 12, height = 5)

library(ggplot2) 
genes_to_check = c('PTPRC', 'CD3D', 'CD3E', 'CD4','CD8A','FOXP3','IL2RA',
                   'CD19', 'CD79A', 'MS4A1' ,
                   'IGHG1', 'MZB1', 'SDC1',
                   'CD68', 'CD163', 'CD14', 
                   'TPSAB1' , 'TPSB2',  # mast cells,
                   'RCVRN','FPR1' , 'ITGAM' ,
                   'C1QA',  'C1QB','MRC1','MSR1' ,  # mac
                   'S100A9', 'S100A8', 'MMP19',  # monocyte
                    'CD86', 'CCL22','S100A4', #DC
                   'LAMP3', 'IDO1','IDO2',## DC3 
                   'CD1E','CD1C', # DC2
                   'KLRB1','NCR1', # NK 
                   'MKI67','TOP2A')
library(stringr)   
DotPlot(merged, features = genes_to_check,assay='RNA',group.by = "seurat_clusters"  ) + coord_flip()

###### 
merged$celltype.main <- recode(merged$seurat_clusters,
                             "0" = "CD4 TCM",
                             "1" = "B cell1",
                             "2" = "CD4 Naive",
                             "3" = "CD8 effector",
                             "4" = "CD8 TCM",
                             "5" = "NK",
                             "6" = "Inflammatory_Monocyte",
                             "7" = "CD8 TEM",
                             "8" = "Classical_Monocyte",
                             "9" = "B cell2",
                             "10" = "NKT",
                             "11" = "Platelet1",
                             "12" = "Plasma cell1",
                             "13" = "Neutrophil",
                             "14" = "Platelet2",
                             "15" = "Plasma cell2",
                             "16" = "pDC",
                             "17" = "CD8 effector")

table(merged$celltype.main)
Idents(merged) <- "celltype.main"
p1 <- DimPlot(merged, label = T, label.size = 5,group.by = "seurat_clusters",reduction = "tsne") + NoLegend()
p2 <- DimPlot(merged, group.by = "celltype.main", label = T, label.size = 4,reduction = "tsne")
p2
p <- p1|p2
p
ggsave("process_data/Celltype.pdf", p, width = 12, height = 5)


table(merged$orig.ident)
p1<-DimPlot(subset(merged, subset = orig.ident=='HC'),reduction = "tsne",pt.size = 1,label = TRUE)+ NoLegend()
p2<-DimPlot(subset(merged, subset = orig.ident=='PT'),reduction = "tsne",label = TRUE,pt.size = 1)
p <- p1|p2
ggsave("process_data/Celltype_group.pdf", p, width = 12, height = 5)

merged@meta.data$sample <- paste(merged@meta.data$orig.ident, merged@meta.data$celltype.main, sep = "_")
table(merged$orig.ident,merged$celltype.main)

VlnPlot(merged, features = c('CD68','FCN1','CD1C','TPSB2','CD14','MARCO','CXCR2',
                             'CLEC9A','IL3RA','KLRD1','MS4A1','IGHG4',
                             'CD3D','CD8A','KLRF1',
                             'CD79A','IGHG4','MS4A1'),group.by = "sample",pt.size = 0)

cuscolors <-c('#E5D2DD','#53A85F','#F1BB72', '#F3B1A0','#D6E7A3','#57C3F3', '#8C549C',
'#476D87','#E59CC4','#AB3282', '#23452F','#BD956A','#8C549C', '#585658','#E95C59','#53A85F','#F1BB72', '#F3B1A0')

v <- VlnPlot(merged, features = c("GZMK","GZMB","IFNG","PRF1","TIGIT","HAVCR2","PDCD1","LAG3"), 
             stack=TRUE, group.by = "celltype.main", cols=cuscolors ) + NoLegend()##

merged$orig.ident=factor(merged$orig.ident,levels = c("Young","Old"))   #20210120
#T cell
VlnPlot(merged, features = c("GZMA","GZMK","GZMB","GZMH","IFNG","PRF1"),split.plot = T,
                   split.by = "orig.ident",pt.size = 0,cols = c("#5858FA","#FE2EC8"),ncol = 3)
#Mono
VlnPlot(pbmc, features = c("IL1B","IL6","TNF","CCL2","CCL3","CCL4"),split.plot = T,
        split.by = "type",pt.size = 0,cols = c("#5858FA","#FE2EC8"),ncol = 3)
#Neuteophils
VlnPlot(pbmc, features = c("S100A4","S100A6","S100A8","S100A9","S100A10","S100A11","S100A12","CSF2","CXCL10"),split.plot = T,
        split.by = "type", group.by = "celltype", pt.size = 0,cols = c("#5858FA","#FE2EC8"),ncol = 3)

VlnPlot(merged, features = c("S100A10"),split.plot = T,split.by = "type",pt.size = 0,cols = c("#5858FA","#FE2EC8"))
ggsave("process_data/Violin_marker3.pdf", width = 12, height = 12)

##存储结果
saveRDS(merged, file = "process_data\\PBMC_all.rds")

#20220222开始-----------------------------
rm(list=ls())
library(dplyr)
library(Seurat)

merged <- readRDS("process_data/PBMC_all.rds")
head(merged@meta.data)
table(merged$seurat_clusters)
median(merged$nFeature_RNA)

table(Idents(merged))

#画比例图-----------------------------------
library(IOBR)
library(stringr)
library(readr)
library(RColorBrewer)
library(scales)
library(ggplot2)
library(cowplot)

use_colors = c("#e5192c","#3a77b7","#3cac4c","#813c93","#f36c24",
               "#37b8c3","#a54922","#6b7627","#28996b",
               "#965b6a","#e9148f","#595b5e","#76c3ad",
               "#80d08a","#d29099","#f2e010","#e9148f","#595b5e","#76c3ad",
               "#80d08a","#d29099","#f2e010")
show_col(use_colors)

table(pbmc$group)
pbmc$celltype <- factor(pbmc$celltype,levels = c("CD4_TCM","CD4_TEM","CD4_CTL","Treg","CD8_TCM","CD8_TEM",
                                                 "CD8_CTL","NKT","proliferating_T","Other_T"))
pbmc$celltype <- factor(pbmc$celltype,levels = c("Naive_B","CD83_B","Memory_B","Activated_memoryB","Plasma_cell","Transitional_B"))
pbmc$celltype <- factor(pbmc$celltype,levels = c("Classical_CD56d","Inflamed_CD56d","HLA_CD56d","CD56b"))
pbmc$Group <- factor(pbmc$Group,levels = c("Cd","Tg","Ad","Ed","Sc"))
pbmc$celltype <- factor(pbmc$celltype,levels = c("Naive_CD8T","IFN_CD8T","Effector_CD8T","Exhausted_CD8T","Dim"))
table(pbmc$celltype)

pbmc$celltype <- factor(pbmc$celltype,levels = c("CD14_Mono1", "CD14_Mono2", "CD16_Mono", "cDC",
                                                 "pDC","Classical_CD56d","Inflamed_CD56d","HLA_CD56d","CD56b",
                                                 "Naive_B","CD83_B","Memory_B","Activated_memoryB","Plasma_cell","Transitional_B",
                                                 "CD4_TCM","CD4_TEM","CD4_CTL","Treg","CD8_TCM","CD8_TEM", "CD8_CTL","NKT","proliferating_T","Other_T"))

scRNA <- pbmc
bar.df=scRNA@meta.data
text.df=as.data.frame(table(bar.df$group))
p1=bar.df%>%ggplot(aes(x=group))+geom_bar(aes(fill=celltype))+
  scale_x_discrete("")+
  scale_y_continuous("cell number",expand = c(0.02,0))+
  scale_fill_manual(values = use_colors)+
  geom_text(data = text.df, aes(x=Var1,y=Freq,label=Freq),size=4)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank())
p1

bar.df$celltype <- factor(bar.df$celltype,levels = c("Naive_CD8T","IFN_CD8T","Effector_CD8T","Exhausted_CD8T","Dim"))
p2 = bar.df %>% ggplot(aes(x = Group)) + geom_bar(aes(fill = celltype), position = "fill")+
  scale_x_discrete("")+
  scale_y_continuous("cell ratio", expand = c(0.02,0), labels = scales::label_percent())+
  scale_fill_manual(values = use_colors)+
  theme_bw()+
  theme(
    panel.grid = element_blank(),
    axis.ticks.x = element_blank())
p2
p1|p2
ggsave2("percentage.pdf", path = "subT", width = 8, height = 5)

library(ggplot2)
epi_clusters <- FetchData(pbmc, vars = c("celltype", "type"))
head(epi_clusters)
count_tumor <- epi_clusters %>% filter(type == "PT") %>% dplyr::count() %>% as.numeric() 
count_normal <- epi_clusters %>% filter(type == "HC") %>% dplyr::count() %>% as.numeric() 


epi_counts <- epi_clusters %>% group_by(type) %>% dplyr::count(celltype)

proportion_tumor <- epi_counts %>% filter(type == "PT") %>% mutate(proportion = n/count_tumor)
proportion_normal <- epi_counts %>% filter(type == "HC") %>% mutate(proportion = n/count_normal)

proportion_epi <- full_join(proportion_normal, proportion_tumor, by = "celltype") %>% 
  mutate(proportion.x = ifelse(is.na(proportion.x), 0,  proportion.x)) %>%  
  mutate(proportion.y = ifelse(is.na(proportion.y), 0,  proportion.y)) %>%
  mutate(tissue_type.x = "HC") %>%
  mutate(tissue_type.y = "PT") %>%
  mutate(cluster_type = ifelse(proportion.x > proportion.y, "HC", "PT"))

cluster_type_data <- left_join(x = epi_clusters, y = proportion_epi, by = "celltype")
rownames(cluster_type_data) <- rownames(epi_clusters)

epi <- AddMetaData(merged, dplyr::select(cluster_type_data, cluster_type))
dim(cluster_type_data)
head(epi@meta.data)

### Bar plot for figure 2
n1 <-  dplyr::select(proportion_epi, c(tissue_type.x, celltype, proportion.x)) %>%
  mutate(tissue_type = tissue_type.x) %>% 
  mutate(proportion = proportion.x) %>%
  mutate(tissue_type.x = NULL) %>%
  mutate(proportion.x = NULL)

t1 <-  dplyr::select(proportion_epi, c(tissue_type.y, celltype, proportion.y)) %>%
  mutate(tissue_type = tissue_type.y) %>% 
  mutate(proportion = proportion.y) %>%
  mutate(tissue_type.y = NULL) %>%
  mutate(proportion.y = NULL)

proportion_epi2 <- rbind(n1, t1)

proportion_epi2$tissue_type = factor(proportion_epi2$tissue_type,levels = c("HC", "PT"))
ggplot(proportion_epi2, aes(fill = tissue_type, y = proportion, x = celltype)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  scale_fill_manual(values =  c("#5858FA","#FE2EC8"))

ggsave("figs/10x_proportion.pdf", width = 12, height = 5)
table(merged$orig.ident,merged$celltype.main)
p <- as.data.frame(proportion_epi2)
write.csv(p,"process_data/liver/propations.csv",row.names = F)

