rm(list=ls())
setwd("E:/AOM_tumor/filter_data/cellchat")


##Load R packages
{
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(openxlsx)
}

# Part1 data CellChat object initialize
## Import data
imm <- readRDS("colon_imm_seurat.rds")

bc <- readRDS("5bc_seurat.rds")

tc <- readRDS("TC5_seurat.rds")


rm(NT_BT)

#merge data
all <- merge(imm, y = c(bc, tc), add.cell.ids = c("4K", "8K"), project = "PBMC12K",merge.data = TRUE)

table(NT_BT$orig.ident)
table(NT_BT$group)
table(NT_BT$hash.ID)
table(NT_BT$cell_type)

saveRDS(NT_BT, file = "NT_BT_seurat.rds")

#Extract normal samples
N_BT <- subset(x = NT_BT, subset = group == "normal")

saveRDS(N_BT, file = "N_BT_seurat.rds")

#Extract tumor samples
T_BT <- subset(x = NT_BT, subset = group == "tumor")

saveRDS(T_BT, file = "T_BT_seurat.rds")


####Import data analysis
rm(list=ls())
setwd("E:/spleen_alldata/cellchat/tumor_BT")

T_BT <-  readRDS("T_BT_seurat.rds")
T_BT
table(T_BT$orig.ident)
table(T_BT$cell_type)

#analysis wt
tb_ko <- T_BT[,T_BT@meta.data[["orig.ident"]]=='T-ko']
exp <- tb_ko@assays$RNA@data
meta <- tb_ko@meta.data

index1 <- which(colnames(exp) %in% rownames(meta))
data.input <- exp[,index1]
dim(data.input)

## create CellChat object
cellchat <- createCellChat(object = data.input,
                           meta = meta,
                           group.by = "cell_type")

## import data
CellChatDB <- CellChatDB.mouse  #CellChatDB.human,CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)# Show the structure of the database
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB

## data preprocessing
cellchat <- subsetData(cellchat)# subset the expression data of signaling genes for saving computation cost
#future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)# project gene expression data onto PPI network (optional)
#- between cells pathways stored in "net" "netP".

# Part2 between cells
## calculate
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

## extract
df.net <- subsetCommunication(cellchat)
#df.net <- subsetCommunication(all_cellchat, signaling = 'WNT')#netp pathways data /sources.use /signaling pathways

## pathways between cells
cellchat <- computeCommunProbPathway(cellchat)

## calculate integrate
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# save CellChat object
cellchat@netP$pathways

saveRDS(cellchat, file = "cellchat_T_BT_wt.rds")

rm(list=ls())


## Import data
cellchat_wt <- readRDS(file = "cellchat_T_BT_wt.rds")
cellchat_ko <- readRDS(file = "cellchat_T_BT_ko.rds")


#merge cellchat object
cco.list <- list(T_wt=cellchat_wt, T_ko=cellchat_ko)

cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)


#all cells: count
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p

ggsave("Overview_number_strength.png", p, width = 6, height = 4)


#count
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#count heatmap
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2

#count
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}


#specific pathways identify Visualization
## pathways analysis
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
p
ggsave("Compare_pathway_strengh.png", p, width = 10, height = 6)


#display
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),  comparison = c(1, 2), angle.x = 45)
p
ggsave("Compare_bubble.png", p, width = 12, height = 8)

#display upregulated downregulated
p1 <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in T_ko", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Increased signaling in T_ko.png", p1, width = 12, height = 7)


#####pathways display
## inspect pathways
cellchat_ko@netP$pathways
cellchat_wt@netP$pathways

color <- c("#c04327","#dd6f56","#e79a88","#c491a5","#8c8ec1","#2a4986",
           "#8189b0")

color <- c('#800020','#9A3671','#E1C392','#25BCCD','#95A238','#FADA5E',
           '#B05923',"#c04327","#dd6f56","#e79a88","#c491a5","#8c8ec1","#2a4986",
           "#8189b0") 
source <- c('Cd8+ NKT',
            'Effector Memory Cd4+ TC',
            'Follicular BC',
            'GC BC',
            'IFI+ BC',
            'MZ BC',
            'Naive Cd4+ TC',
            'Naive Cd8+ TC',
            'Plasma cell',
            'Plasmblast',
            'Pro/Pre BC',
            'Tfh',
            'Transitional BC',
            'Tregs')
target <- c('Cd8+ NKT',
            'Effector Memory Cd4+ TC',
            'Follicular BC',
            'GC BC',
            'IFI+ BC',
            'MZ BC',
            'Naive Cd4+ TC',
            'Naive Cd8+ TC',
            'Plasma cell',
            'Plasmblast',
            'Pro/Pre BC',
            'Tfh',
            'Transitional BC',
            'Tregs')


## Visualization pathways
CD40 <- c("CD40")
MHCI <- c("MHC-I")

CD22 <- c("CD22")
CD45 <- c("CD45")
GALECTIN <- c("GALECTIN")
PECAM1 <- c("PECAM1")
SELPLG <- c("SELPLG")
LAMININ <- c("LAMININ")
SEMA4 <- c("SEMA4")
TGFb <- c("TGFb")
CD80 <- c("CD80")
SEMA7 <- c("SEMA7")

## 
#CD52
netVisual_chord_gene(cellchat_ko, big.gap = 20,small.gap = 5,
                     targets.use = target,transparency = 0.5,
                     signaling = MHCI,show.legend = T,color.use = color,
                     title.name = 'MHCI signaling pathway interaction in ko')

netVisual_chord_gene(cellchat_wt, big.gap = 20,small.gap = 5,
                     targets.use = target,transparency = 0.5,
                     signaling = MHCI,show.legend = T,color.use = color,
                     title.name = 'MHCI signaling pathway interaction in wt')



## heatmap
#CD52
netVisual_heatmap(cellchat_ko, color.heatmap = 'Reds',color.use = color,
                  signaling = MHCI,targets.use = target,
                  sources.use = source,
                  title.name = 'MHC-I signaling pathway interaction in ko')

p2 = netVisual_heatmap(cellchat_wt, color.heatmap = 'Reds',color.use = color,
                  signaling = MHCI,targets.use = target,
                  sources.use = source,
                  title.name = 'MHC-I signaling pathway interaction in wt')
p2

ggsave("MHCI_wt_heatmap.pdf",p2,width = 12, height = 7)

## 

netVisual_bubble(cellchat_ko,
                 signaling = "MHC-I",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MHC-I signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "MHC-I",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MHC-I signaling pathway in wt')

#MIF
netVisual_bubble(cellchat_ko,
                 signaling = "MIF",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MIF signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "MIF",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MIF signaling pathway in ko')

#CD22
netVisual_bubble(cellchat_ko,
                 signaling = "CD22",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'CD22 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "CD22",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'wt-CD22 signaling pathway in wt')
#CD45
netVisual_bubble(cellchat_ko,
                 signaling = "CD45",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'CD45 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "CD45",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'wt-CD45 signaling pathway in wt')
#GALECTIN
netVisual_bubble(cellchat_ko,
                 signaling = "GALECTIN",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'GALECTIN signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "GALECTIN",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'GALECTIN signaling pathway in wt')
#PECAM1
netVisual_bubble(cellchat_ko,
                 signaling = "PECAM1",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'PECAM1 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "PECAM1",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'PECAM1 signaling pathway in wt')

##SELPLG
netVisual_bubble(cellchat_ko,
                 signaling = "SELPLG",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SELPLG signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "SELPLG",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SELPLG signaling pathway in wt')

#LAMININ
netVisual_bubble(cellchat_ko,
                 signaling = "LAMININ",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'LAMININ signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "LAMININ",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'LAMININ signaling pathway in wt')
##SEMA4
netVisual_bubble(cellchat_ko,
                 signaling = "SEMA4",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SEMA4 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "SEMA4",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SEMA4 signaling pathway in wt')



#prob data label heatmap
KO_MHCI <- cellchat_ko@netP[["prob"]][,,'MHC-I'] 
WT_MHCI <- cellchat_wt@netP[["prob"]][,,'MHC-I'] 

KO_MIF <- cellchat_ko@netP[["prob"]][,,'MIF'] 
WT_MIF <- cellchat_wt@netP[["prob"]][,,'MIF'] 

KO_CD22 <- cellchat_ko@netP[["prob"]][,,'CD22'] 
WT_CD22 <- cellchat_wt@netP[["prob"]][,,'CD22'] 

KO_CD45 <- cellchat_ko@netP[["prob"]][,,'CD45'] 
WT_CD45 <- cellchat_wt@netP[["prob"]][,,'CD45'] 

KO_GALECTIN <- cellchat_ko@netP[["prob"]][,,'GALECTIN'] 
WT_GALECTIN <- cellchat_wt@netP[["prob"]][,,'GALECTIN']

KO_PECAM1 <- cellchat_ko@netP[["prob"]][,,'PECAM1'] 
WT_PECAM1 <- cellchat_wt@netP[["prob"]][,,'PECAM1']

KO_SELPLG <- cellchat_ko@netP[["prob"]][,,'SELPLG'] 
WT_SELPLG <- cellchat_wt@netP[["prob"]][,,'SELPLG']

list_of_datasets <- list('ko_MHCI' = KO_MHCI, 
                         'wt_MHCI' = WT_MHCI)

list_of_datasets <- list('ko_MHCI' = KO_CD52, 
                         'wt_MHCI' = WT_CD52,
                         'ko_MIF' = KO_MIF,
                         'wt_MIF' = WT_MIF,
                         'ko_CD22' = KO_CD22,
                         'wt_CD22' = WT_CD22,
                         'ko_CD45' = KO_CD45,
                         'wt_CD45' = WT_CD45,
                         'ko_GALECTIN' = KO_GALECTIN,
                         'wt_GALECTIN' = WT_GALECTIN,
                         'ko_PECAM1' = KO_PECAM1,
                         'wt_PECAM1' = WT_PECAM1,
                         'ko_SELPLG' = KO_SELPLG,
                         'wt_SELPLG' = WT_SELPLG)

write.xlsx(list_of_datasets, "MHCI_cellchat_pathways_prob.xlsx",rowNames=TRUE)
