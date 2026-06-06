setwd("D:/Desktop/RData/mice_scRNA/GEX+ADT/spleen_identify/B_cell")

##Clear environment
rm(list=ls())

##Load R packages
library(Seurat)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(readxl)

##Load data
bc <- readRDS("bc_seurat_harmony.rds")

table(bc@meta.data$orig.ident)

meta <- bc@meta.data[,c("orig.ident","group")]#Group information for later plotting
bc <- as.matrix(bc@assays$RNA@counts)#Extract count matrix

library(msigdbr)

msigdbr_species() 

mouse <- msigdbr(species = "Mus musculus")

mouse[1:5,1:5]

table(mouse$gs_cat) 

#Extract GO file
table(mouse$gs_subcat)

mouse_GO_bp = msigdbr(species = "Mus musculus",
                      category = "C5", #GO is in C5
                      subcategory = "GO:BP") %>% 
  dplyr::select(gs_name,gene_symbol)#Here you can choose gene symbols or IDs depending on your data needs, mainly for convenience
mouse_GO_bp_Set = mouse_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)#GSVA requires a list later, so convert it to a list

#After the expression matrix and pathway information are ready, GSVA analysis can be run

T_gsva <- gsva(expr = bc, 
               gset.idx.list = mouse_GO_bp_Set,
               kcdf="Poisson", #Check the help function to choose an appropriate kcdf method 
               parallel.sz = 5)

write.table(T_gsva, '5bc_gsva.xlsx', row.names=T, col.names=NA, sep="\t")

T_gsva <- read_excel("D:/Desktop/RData/mice_scRNA/GEX+ADT/spleen_identify/B_cell/Fob/Fob_gsva.xlsx")

#Use the limma package for differential analysis
library(edgeR)
library(limma)

table(bc@meta.data$group)

#Set groups
group <- c(rep("WT", 2167), rep("CKO", 189)) %>% as.factor()#Set groups with the control first
desigN <- model.matrix(~ 0 + group) #Build comparison matrix
colnames(desigN) <- levels(group)

#Set thresholds
#logFCcutoff <- log2(1.5)
logFCcutoff <- log2(1.2)
PvalueCutoff <- 0.05      #Use the unadjusted p value here

#limma pathway differential analysis: mcao vs sham
fit = lmFit(T_gsva, desigN)
fit2 <- eBayes(fit)
diff=topTable(fit2,adjust='fdr',coef=2,number=Inf)#Adjust as needed; topTable

DEgeneSets <- diff[(diff$P.Value < PvalueCutoff & (diff$logFC>logFCcutoff | diff$logFC < (-logFCcutoff))),]

write.csv(DEgeneSets, file = "5bc_GSVA_pathway_1.2.xlsx")

DEgeneSets <- read_excel("Naive_Cd4_GSVA_DEG.xlsx")

#Show heatmap
pheatmap::pheatmap(T_gsva[rownames(DEgeneSets,),],
                   show_rownames = T,
                   show_colnames = T)
p = pheatmap::pheatmap(T_gsva[rownames(DEgeneSets,),],
                       show_rownames = T,
                       show_colnames = T)
ggsave('pro_gsva.png', p,width = 25,height = 15)


#Finally visualize the differential pathways of interest
ko_up <- c("GOBP_PHAGOCYTOSIS_RECOGNITION",
           "GOBP_COMPLEMENT_ACTIVATION",
           "GOBP_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN")
wt_up <- c("GOBP_REGULATION_OF_GAP_JUNCTION_ASSEMBLY",
           "GOBP_MESENCHYMAL_STEM_CELL_PROLIFERATION",
           "GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_STEM_CELL_PROLIFERATION",
           "GOBP_REGULATION_OF_LYMPHANGIOGENESIS",
           "GOBP_MAST_CELL_PROLIFERATION",
           "GOBP_INTERCELLULAR_TRANSPORT",
           "GOBP_NEGATIVE_REGULATION_OF_CELL_FATE_SPECIFICATION",
           "GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT",
           "GOBP_POSITIVE_REGULATION_OF_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY")
TEST <- c(ko_up, wt_up)
diff$ID <- rownames(diff) 
Q <- diff[TEST,]
group1 <- c(rep("ko_up", 3), rep("wt_up", 9)) 
df <- data.frame(ID = Q$ID, score = Q$t,group=group1 )
# Sort by t score
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#Add the pathway ID column

ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #Remove background color
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="t value of GSVA score")+
  scale_fill_manual(values = c("red","blue"))#Set colors

p = ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #Remove background color
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="t value of GSVA score")+
  scale_fill_manual(values = c("#008020","#08519C"))#Set colors

ggsave('PC_gsva.png', p, width = 20,height = 8)
