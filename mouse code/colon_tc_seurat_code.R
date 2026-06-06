setwd("E:/AOM_tumor/filter_data/T_cell")
rm(list=ls())


##Load R packages
{
  library(densityClust)
  library(scran)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(DoubletFinder)
  library(BiocSingular)
  library(sctransform)
  library(glmGamPoi)
  library(scater)
  library(scDblFinder)
  library(gridExtra)
  library(SingleR)
  library(cowplot)
  library(harmony)
  library(devtools)
  library(tidyverse)
  library(scCATCH)
  library(mindr)
  
}

set.seed(123)  

##Import data
exp <- readRDS("TC5_exp.rds")

meta <- read.csv("TC5_meta.csv", row.names = 1)


##Initialize Seurat object
tc <- CreateSeuratObject(counts =exp , project = "T_cell", min.cells = 3, min.features = 200) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = tc@var.genes, npcs = 20, verbose = FALSE)


##Add group information
#Group information
orig.ident <- meta$orig.ident
names(orig.ident) <- colnames(x = tc)
tc <- AddMetaData(
  object = tc,
  metadata = orig.ident,
  col.name = "orig.ident"
)
table(tc@meta.data$orig.ident)

#Sample information
group <- meta$group
names(group) <- colnames(x = tc)
tc <- AddMetaData(
  object = tc,
  metadata = group,
  col.name = "group"
)
table(tc@meta.data$group)

save(tc, file = "5TC.RData")

#load("5TC.RData")

table(tc@meta.data$orig.ident)


p1 <- DimPlot(object = tc, reduction = "pca", pt.size = .1,group.by = "orig.ident")
p2 <- VlnPlot(object = tc, features = "PC_1", pt.size = .1,group.by = "orig.ident")
plot_grid(p1,p2)
p = plot_grid(p1,p2)
ggsave("tc_plot.png", p ,width=9 ,height=6)

##Run Harmony
tc <- tc %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(tc, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = tc, reduction = "harmony", pt.size = .1,group.by = "orig.ident")
p2 <- VlnPlot(object = tc, features = "harmony_1", pt.size = .1,group.by = "orig.ident")
plot_grid(p1,p2)


##Downstream analysis
tc <- FindNeighbors(tc, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.5)

tc <- RunUMAP(tc, reduction = "harmony", dims = 1:20)

DimPlot(tc, reduction = "umap",label = T,pt.size = 0.7) 

##Find highly variable genes
markers <- FindAllMarkers(object = tc, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

write.csv(top30, "5tc_markers.csv")


##marker Visualization
genes_to_check = c("Il7r","Icos","Gzmk","Gzmb", #Memory T cells 0
                   "Chn2","Prkch","Bcl2", #Cd8+T_Bcl2 2/6/11
                   "Ms4a4b","Tcf7", "Cd69","Cd28","Cd44",#Cd4-Activated T 1
                   "Foxp3","Ctla4",#Tregs 3/8/13/15
                   "Il23r","Cd163l1","Trdc",#γδT cell 4/9
                   "Cd8a","Cd7","Nkg7","Klrd1","Gzma",#Cd8+_Cytotoxic T cell 5
                   "Pid1","Mrc1","C1qc","Adgre1",#Mac 
                   "Sell","S100a8", "S100a9",#Neutrophils 10
                   "Gata3","Il1rl1") #ILC 11/12

p <- DotPlot(tc, features = genes_to_check,assay='RNA', col.min = 0, dot.min = 0)  + RotatedAxis() + coord_flip()

p

ggsave("colon_dotplot.png", p, width=8,height=8)


new.cluster.ids <- c("0"="Memory_like T",
                     "1"="Cd4 Activated T",
                     "2"="Cd8 T_Bcl2",
                     "3"="Tregs",
                     "4"="γδT cell",
                     "5"="Cd8 Cytotoxic T",
                     "6"="Cd8 T_Bcl2",
                     "7"="other",
                     "8"="Tregs",
                     "9"="γδT cell",
                     "10"="other",
                     "11"="Cd8 T_Bcl2",
                     "12"="other",
                     "13"="Tregs",
                     "14"="other",
                     "15"="Tregs",
                     "16"="other")
                     
                     
names(new.cluster.ids) <- levels(tc)
tc <- RenameIdents(tc, new.cluster.ids)

table(tc@meta.data$orig.ident,tc@active.ident)

#Check cell counts
table(tc@meta.data$orig.ident,tc@active.ident)

ident <- data.frame(tc@active.ident)

head(ident)

tc@meta.data$cell_type <- ident$tc.active.ident

table(tc@meta.data$cell_type)
table(tc@meta.data$orig.ident, tc@meta.data$cell_type)

##Remove non-immune cells
tc_scRNA <- subset(x=tc, idents =c("Memory_like T","Cd8 T_Bcl2","Cd4 Activated T","Tregs","γδT cell",
                                    "Cd8 Cytotoxic T"))

tc_scRNA

write.csv(tc_scRNA@meta.data, "5tc_meta_harmony.csv")

saveRDS(tc_scRNA, file = "5tc_seurat_harmony.rds")

TC <- readRDS("tc_seurat_harmony.rds")


##ggplot plotting
##Extract the first two principal component coordinates
# extact PC ranges
pc12 <- Embeddings(object = tc_scRNA,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)

##Build labels and positions needed for coordinate axes
# get botomn-left coord
lower <- floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2))) - 2

# get relative line length
linelen <- abs(0.3*lower) + lower

# mid point
mid <- abs(0.3*lower)/2 + lower

# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('UMAP_2','UMAP_1'),each = 2))

# axies label
label <- data.frame(lab = c('UMAP_2','UMAP_1'),angle = c(90,0),
                    x = c(lower - 3,mid),y = c(mid,lower - 2.5))

##Visualization
# plot
color <- c("#f39c90", "#4DBBD5B2","#A9D179","#00A087B2","#F5AE6B","#CC79A7","#4387B5")

plot1 <- DimPlot(tc_scRNA, reduction = 'umap', label = T,repel = TRUE,
                 pt.size = .1, cols = color,label.size = 4) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot1 <- DimPlot(tc_scRNA, reduction = 'umap', label = F,repel = TRUE,
                 pt.size = .1, cols = color,split.by  = "group", label.size = 4,ncol = 2) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot1
ggsave('tc_group_umap_v2.png',plot1,width = 8,height = 8)


##Heatmap marker visualization
##DotPlot
genes_to_check = c("Il7r","Icos","Gzmk","Gzmb", 
                   "Chn2","Prkch","Bcl2", 
                   "Ms4a4b","Tcf7", "Cd69","Cd28","Cd44",
                   "Foxp3","Ctla4",
                   "Il23r","Cd163l1","Trdc",
                   "Cd8a","Cd7","Nkg7","Klrd1","Gzma")

p <- DotPlot(tc_scRNA, features = genes_to_check,assay='RNA',cols = c(low="white",high="darkred"),
             col.min = 0, dot.min = 0)  + RotatedAxis() + coord_flip()
p
ggsave("5tc_marker_dotplot.png", p, width=7,height=7)


##Violin plot
genes_to_check = c("Il7r","Icos","Gzmk","Gzmb", 
                   "Chn2","Prkch","Bcl2", 
                   "Ms4a4b","Tcf7", "Cd69","Cd28","Cd44",
                   "Foxp3","Ctla4",
                   "Il23r","Cd163l1","Trdc",
                   "Cd8a","Cd7","Nkg7","Klrd1","Gzma")



markers <- CaseMatch(genes_to_check, rownames(tc_scRNA))
markers <- as.character(markers)
VlnPlot(tc_scRNA, features = markers, pt.size = 0)

VlnPlot(tc_scRNA, features = markers, pt.size = 0, group.by = 'cell_type', stack = T)+NoLegend()


#FeaturePlot
library(viridis)
pal <- viridis(n = 15, option = "C", direction = -1)

p1 = FeaturePlot(tc_scRNA,features = c("Il7r","Icos","Gzmk","Gzmb", 
                                    "Chn2","Prkch","Bcl2", 
                                    "Ms4a4b","Tcf7", "Cd69","Cd28","Cd44",
                                    "Foxp3","Ctla4",
                                    "Il23r","Cd163l1","Trdc",
                                    "Cd8a","Cd7","Nkg7","Klrd1","Gzma"),cols = pal, order = T, ncol =6)

ggsave("5tc_VlnPlot.png", p1, width=14,height=8)



####Cell proportion plots across groups

##Load R packages
options(stringsAsFactors = F)
{
  library(densityClust)
  library(scran)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(DoubletFinder)
  library(BiocSingular)
  library(sctransform)
  library(glmGamPoi)
  library(scater)
  library(scDblFinder)
  library(gridExtra)
  library(SingleR)
  library(cowplot)
  library(harmony)
  library(devtools)
  library(tidyverse)
  library(scCATCH)
  library(mindr)
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(tidyr)# use gather & spread
  library(reshape2) # use melt & dcast 
  library (gplots) 
  
}

TC <- readRDS("TC_seurat_harmony.rds")
TC

phe=tc_scRNA@meta.data
colnames(phe)
tb=table(phe$cell_type,
         phe$orig.ident)

head(tb)

#Count cells
balloonplot(tb, main ="tumor-T cells", xlab ="celltype", ylab="sample",
            label = T, show.margins = T)


bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
write.csv(bar_per,file = "celltype_by_group_percent.csv")
col =c("#BD6263","#8EA325","#3FA116","#CE2820","#9265C1","#885649","#2EBEBE")
ggplot(bar_per, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "right",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = "% Relative cell source", fill = NULL)+labs(x = NULL)+
  scale_fill_manual(values=col)

ggsave("5tc_celltype.png",width = 8,height = 4)



