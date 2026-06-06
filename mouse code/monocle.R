setwd("E:/spleen_alldata/NT_spleen_BC")
rm(list=ls())


##Load R packages
{
  library(Seurat)
  library(cowplot)
  library(harmony)
  library(dplyr)
  library(tidyverse)
  library(patchwork)
  library(BiocManager)
  library(monocle)
  library(patchwork)
}

###Import data
bc <- readRDS("bc_seurat_harmony.rds")

expr_matrix <- as(as.matrix(bc@assays$RNA@counts), 'sparseMatrix')

p_data <- bc@meta.data

p_data$celltype <- bc@active.ident 

f_data <- data.frame(gene_short_name = row.names(bc),row.names = row.names(bc))

#Build CDS object
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

#Convert p_data and f_data from data.frame objects to AnnotatedDataFrame objects.
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 3)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)


###Step 5 trajectory-defining genes select Visualization build trajectory
#1. Use highly variable genes selected by Seurat
express_genes <- VariableFeatures(bc)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

#2. Use cluster differentially expressed genes
deg.cluster <- FindAllMarkers(bc)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

#3. Use highly variable genes selected by Monocle
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)


#4. dpFeature method
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1)
head(diff)

deg <- subset(diff, qval < 0.01) #Select 2,829 genes
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)

#Around 2,000 ordering genes is usually appropriate; if there are too many genes, select the top genes
ordergene <- row.names(deg)[order(deg$qval)][1:400]



##Step 6 (select for cell ordering, use reversed graph embedding (DDRTree) data perform)
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)

#Use the root_state parameter to set the pseudotime root. In the pseudotime-colored plot below, the left side is the root. The state plot shows the root is State1; to set the other end as the root, do the following
cds <- orderCells(cds, root_state = 5) #Set State5 as the starting point of the pseudotime axis


##Save metadata and pseudotime object
saveRDS(cds,file = "bc_monocle_WT_DEG.rds")
write.csv(pData(cds),"bc_monocleMeta_WT_DEG.csv")

##Plotting
pData(cds)$State <- as.factor(pData(cds)$State)
pData(cds)$pseudotime <- as.factor(pData(cds)$pseudotime)

p1 <- plot_cell_trajectory(cds, color_by="Pseudotime")
p1
ggsave("Pseudotime.pdf",p1,width = 8,height = 6)


plot_cell_trajectory(cds, color_by="State", show_tree=F, 
                     show_branch_points = F,show_cell_names =F,
                     show_state_number = F,show_backbone = T)

p2 <- plot_cell_trajectory(cds, color_by="cell_type", show_tree=F, 
                           show_branch_points = F,show_cell_names =F,
                           show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#f39c90","#4d62ae","#b37eb7"))
p2
ggsave("Pseudotime_celltype.pdf",p2,width = 8,height = 6)


p3 <- plot_cell_trajectory(cds, color_by="pseudotime", show_tree=F, 
                           show_branch_points = F,show_cell_names =F,
                           show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#c04327","#dd6f56","#e79a88",
                                 "#c491a5","#8c8ec1","#2a4986"))
p3
ggsave("Pseudotime_after.pdf",p3,width = 8,height = 6)

