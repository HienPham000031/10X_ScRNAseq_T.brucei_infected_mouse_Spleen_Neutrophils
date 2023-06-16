### REMOVE EMPTY DROPLET/AMBIENT RNA FOR SEPERATED CELLRANGER COUNT DATA (NAIVE and 14DPI)
library(SoupX)
# Typical SoupX workflow for 10X data which would produce a new matrix that has had
#the contaminating reads removed
# Load data and estimate soup profile
sc_naive = load10X("./Naive")
sc_14dpi = load10X("./14dpi")
# Estimate rho
sc_naive = autoEstCont(sc_naive)
sc_14dpi = autoEstCont(sc_14dpi)
# Clean the data
Naive_out = adjustCounts(sc_naive)
D14_out = adjustCounts(sc_14dpi)

### CREATE SEURAT OBJECT AND PERFORM QCs FILTERING
library(Seurat)
library(dplyr)
library(tidyverse)
# Loading output of SoupX into Seurat
naive <- CreateSeuratObject(counts = Naive_out, project = "Naive", min.cells = 3, min.features = 200)
D14 <- CreateSeuratObject(counts = D14_out, project = "14dpi", min.cells = 3, min.features = 200)
# How big is the matrix: report number of genes (rows) and number of cells (columns)
dim(naive)
dim(D14)
# Exclude ribosomal RNA/non-coding RNA contamination
#Naive
counts_naive <- GetAssayData(naive, assay = "RNA")
counts_naive <- counts_naive[-(which(rownames(counts_naive) %in% c("Gm42418","Gm26917","AY036118"))),]
naive <- subset(naive, features = rownames(counts_naive))
#14dpi
counts_14dpi <- GetAssayData(D14, assay = "RNA")
counts_14dpi <- counts_14dpi[-(which(rownames(counts_14dpi) %in% c("Gm42418","Gm26917","AY036118"))),]
D14 <- subset(D14, features = rownames(counts_14dpi))
## Standard pre-processing workflow, perform for 2 Seurat object seperately
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
naive[["percent.mt"]] <- PercentageFeatureSet(naive, pattern = "^mt-")
D14[["percent.mt"]] <- PercentageFeatureSet(D14, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(naive, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(D14, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
## We filter low quality cells 
# Naive: cells that have unique feature counts over 4,000 or less than 200 and >10% mitochondrial counts
naive <- subset(naive, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
# 14dpi: cells that have unique feature counts over 5,000 or less than 200 and >10% mitochondrial counts
D14 <- subset(D14, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

###THEN INTEGRATING THE 2 DATA SETS BY SEURAT
# Creating list of 2 data sets:
brucei.list <- list(naive, D14)
# Normalize and identify variable features for each dataset independently
brucei.list <- lapply(X = brucei.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
## Select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = brucei.list)
brucei.list <- lapply(X = brucei.list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
brucei.list
# Perform Integration
brucei.anchors <- FindIntegrationAnchors(object.list = brucei.list, 
                                         anchor.features = features, reduction = "rpca")
# Create list of common genes to keep
to_integrate <- lapply(brucei.list, rownames)
to_integrate
# Integrate data and keep full geneset and creates an 'integrated' data assay
brucei.integrated <- IntegrateData(anchorset = brucei.anchors)
# Perform Cell Cycle Scoring for combined data: PERFORM CELL CYCLE SCORING FOR INTEGRATED DATA
# Assigining list of cell cycle markers
s.genes <- c("Mcm5","Pcna","Tyms","Fen1","Mcm2","Mcm4","Rpm1","Ung","Gins2","Mcm6",
             "Cdca7","Dtl","Prim1","Uhrf1","Mlf1ip","Hells","Rfc2","Rpa2","Nasp",
             "Rad51ap1","Gmnn","Wdr76","Slbp","Ccne2","Ubr7","Pold3","Msh2","Atad2",
             "Rad51","Rrm2","Cdc45","Cdc6","Ex01","Tipin","Dscc1","Blm","Casp8ap2",
             "Usp1","Clspn","Pola1","Chaf1b","Brip1","E2f8")
g2m.genes<-c("Hmgb2","Cdk1","Nusa1","Ube2c","Birc5","Tpx2","Top2a","Ndc80","Cks2",
             "Nuf2","Cks1b","Mki67","Tmpo","Cenpf","Tacc3","Fam64a","Smc4","Ccnb2",
             "Ckap2l","Ckap2","Aurkb","Bub1","Kif11","Anp32e","Tubb4b","Gtse1","Kif20b",
             "Hjurp","Cdca3","Hn1","Cdc20","Ttk","Cdc25c","Kif2c","Rangap1","Ncapd2",
             "Dlgap5","Cdca2","Cdca8","Ect2","Kif23","Hmmr","Aurka","Psrc1","Anln",
             "Lbr","Ckap5","Cenpe","Ctcf","Nek2","G2e3","Gas2l3","Cbx5","Cenpa")
#Assign CellCycle Score
DefaultAssay(brucei.integrated) <- "RNA"
CellCycle.integrated<-CellCycleScoring(brucei.integrated, s.features = s.genes, 
                                       g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(CellCycle.integrated[[]])

# Visualize the distribution of cell cycle markers across
RidgePlot(CellCycle.integrated, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2)
# Identify the most variable genes
CellCycle.integrated <- FindVariableFeatures(CellCycle.integrated, 
                                           selection.method = "vst",
                                           nfeatures = 2000, 
                                           verbose = FALSE)
# Scale the counts:
DefaultAssay(CellCycle.integrated) <- "integrated"
CellCycle.integrated<- ScaleData(CellCycle.integrated)
# Perform PCA
CellCycle.integrated <- RunPCA(CellCycle.integrated)
# Plot the PCA colored by cell cycle phase
DimPlot(CellCycle.integrated,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
DimPlot(CellCycle.integrated)
##  Here we do see some differences amongs G1, G2M and S phase
#Based on this plot, we would regress out the variation due to cell cycle during
#the scaling step of integrated data
#It is suggested to regress out the difference between the G2M and S phase scores. 
#This means that signals separating non-cycling cells and cycling cells will be maintained, 
#but differences in cell cycle phase among proliferating cells (which are often uninteresting), 
#will be regressed out of the data
brucei.integrated$CC.Difference <- CellCycle.integrated$S.Score - CellCycle.integrated$G2M.Score
head(brucei.integrated@meta.data)

### RUN THE STANDARD WORKFLOW FOR SCALING AND CLUSTERING
# specify that we will perform downstream analysis on the corrected data note that the original
# unmodified data still resides in the 'RNA' assay
DefaultAssay(brucei.integrated) <- "integrated"
# Scaling the data
brucei.integrated <- ScaleData(brucei.integrated, vars.to.regress = "CC.Difference", verbose = FALSE)
# Perform linear dimensional reduction and Run non-linear dimensional reduction (UMAP)
brucei.integrated <- RunPCA(brucei.integrated, npcs = 100, verbose = FALSE)
DimHeatmap(brucei.integrated, dims = 1:20, cells = 500, balanced = TRUE)
DimHeatmap(brucei.integrated, dims = 21:40, cells = 500, balanced = TRUE)
#Determine the 'dimensionality' of the dataset
JackStraw.integrated<-JackStraw(brucei.integrated, dims = 50)
JackStraw.integrated<- ScoreJackStraw(JackStraw.integrated, dims = 1:50)
JackStrawPlot(JackStraw.integrated, dims = 1:50)
ElbowPlot(JackStraw.integrated)

brucei.integrated <- RunUMAP(brucei.integrated, reduction = "pca", dims = 1:30)
brucei.integrated <- FindNeighbors(brucei.integrated, reduction = "pca", dims = 1:30)
brucei.integrated <- FindClusters(brucei.integrated, resolution = 0.5)

# Visualization
DimPlot(brucei.integrated, reduction = "umap", label = TRUE)
head(brucei.integrated@meta.data)
DefaultAssay(brucei.integrated) <- "RNA"
DotPlot(brucei.integrated, features = c("Ly6g","Ly6c2","Cd177","Itgam","Csf1r","Ms4a3","Cebpe","Mpo","Cxcr2","Cxcr4"), 
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text = element_text(size = 25, face = "italic"),
                                                        legend.title = element_text(size = 20),
                                                        legend.text = element_text(size = 15),
                                                        axis.title.x = element_text(size = 25),
                                                        axis.title.y.left = element_text(size = 25))

saveRDS(brucei.integrated,"./brucei.integrated.rds")

### CELLTYPE ANNOTATION BY SINGLE-R
library(SingleR)
library(celldex)
library(scran)
# Loading the reference library annotation
ref <- ImmGenData()
ref
# Perform Annotation
brucei.annotation <- SingleR(test = as.SingleCellExperiment(brucei.integrated), ref = ref,
                             assay.type.test = 1, labels = ref$label.main)
brucei.annotation
# Copy over the labels and pruned.labels/Align the annotation to Seurat object
brucei.integrated$SingleR.pruned.labels<-brucei.annotation$pruned.labels
brucei.integrated$SingleR.labels<-brucei.annotation$labels
DimPlot(brucei.integrated, reduction = "umap", group.by = "SingleR.labels", label = TRUE, 
        pt.size = 0.5, label.size = 9)
##1st SUBSET Neutrophils population
Idents(brucei.integrated) <- "SingleR.labels"
Neutrophils1 <- subset((brucei.integrated), idents = "Neutrophils")
head(Neutrophils1@meta.data)
##Annotation, but with 'label.fine'
Neutrophils1.annotation <- SingleR(test = as.SingleCellExperiment(Neutrophils1), ref = ref,
                                   assay.type.test = 1, labels = ref$label.fine)
Neutrophils1$SingleR.Neutrophils1.labels <- Neutrophils1.annotation$labels
Neutrophils1 <- RunUMAP(Neutrophils1, dims = 1:10)
DimPlot(Neutrophils1, reduction = "umap", group.by = "SingleR.Neutrophils1.labels", label = TRUE, pt.size = 0.5)
##2nd SUBSET of Neutrophils
head(Neutrophils1@meta.data)
Idents(Neutrophils1) <- "SingleR.Neutrophils1.labels"
Neutrophils2 <- subset((Neutrophils1), idents = c("Neutrophils (GN.URAC)",
                                                  "Neutrophils (GN.ARTH)",
                                                  "Neutrophils (GN)",
                                                 "Neutrophils (GN.Thio)"))
DefaultAssay(Neutrophils2) <- "integrated"
# Perform linear dimensional reduction and Run non-linear dimensional reduction (UMAP)
Neutrophils.integrated <- RunPCA(Neutrophils2, npcs = 50, verbose = FALSE)
DimHeatmap(Neutrophils.integrated, dims = 1:20, cells = 500, balanced = TRUE)
Neutrophils.integrated <- RunUMAP(Neutrophils.integrated, reduction = "pca", dims = 1:4)
Neutrophils.integrated <- FindNeighbors(Neutrophils.integrated, reduction = "pca", dims = 1:4)
Neutrophils.integrated <- FindClusters(Neutrophils.integrated, resolution = 0.1)
# Visualization
DimPlot(Neutrophils.integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(Neutrophils.integrated, reduction = "umap", split.by = "orig.ident", label = TRUE)
DimPlot(Neutrophils.integrated, reduction = "umap")
#Assign new cluster IDs
new.cluster.ids <- c("N3","N2","N4","N1")
names(new.cluster.ids) <- levels(Neutrophils.integrated)
new.cluster.ids
Neutrophils.integrated <- RenameIdents(Neutrophils.integrated, new.cluster.ids)
levels(Neutrophils.integrated) <- c("N1", "N2", "N3", "N4")
#Add new Identity as SeuratObject level
Neutrophils.integrated$Identity<-Idents(Neutrophils.integrated)
head(Neutrophils.integrated)
table(Idents(Neutrophils.integrated), Neutrophils.integrated$orig.ident)
#Perform Differential Expressed Genes (DEGs) test
DefaultAssay(Neutrophils.integrated) <- "RNA"
Neutrophil.markers <- FindAllMarkers(Neutrophils.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = log(1.5))
Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Neutrophil.markers<- subset(Neutrophil.markers, Neutrophil.markers$p_val<0.05)
Neutrophil.markers
nrow(Neutrophil.markers)
top10<- Neutrophil.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#Heatmap of top 10 DEGs per clusters
all.genes <- rownames(Neutrophils.integrated)
Neutrophils.scaled <- ScaleData(Neutrophils.integrated, features = all.genes)
DoHeatmap(Neutrophils.scaled, features = top10$gene, label = TRUE, size = 7) + NoLegend() + scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')),
                                                                                                                       mid = "white", 
                                                                                                                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                                                                                                                        midpoint = 0, guide = "colourbar",
                                                                                                                        aesthetics = "fill") + theme(axis.text = element_text(size = 21, face = "italic"),
                                                                                                                                                     legend.title = element_text(size = 20),
                                                                                                                                                     legend.text = element_text(size = 15))
# Assign genes lists
TFs <- c("Gata1","Gata2","Cebpa","Runx1","Gfi1","Pclaf","Top2a","Ube2c","Ube2s",
         "Cebpe","Cebpb","Cebpd","Spi1","Fos","Junb","Jund","Irf1",
         "Irf2","Irf5","Irf7","Irf9","Mapk13","Socs13","Pten","Ms4a6c","Ms4a3",
         "Ms4a4c","Pou2af1","Ms4a6b")
Aging <- c("Sell","Itgam","Itga4","Cxcr4","Cxcr2","Cd47","Cd24a","Tlr4","Icam1",
           "Itgax","Cebpd","Spi1")
AntiApoptosis <- c("Mcl1","Bcl2","Bcl2l1","Bclxl","Bcl2a1b","Bcl10","Xiap","Birc5",
                   "Csf2ra","Csf3r","Csf2rb","Xiap","S100a8","S100a9","Cflar",
                   "Gadd45","Traf1","Sdf1","Mki67")
ProApoptosis <- c("Bad","Bax","Bak1","Casp2","Casp3","Casp7","Casp8","Casp9","Casp10",
                  "Casp12","Fas","Tnfrsf6","Tnfsf10","Tnfrsf10a","Tnfrsf1a","Bbc3",
                  "Pmaip1","Diablo","Bcl2l11","Apaf1","Endog","Pidd1")
NADPH <- c("Cybb","Cyba","Rac2","Rac1","Ncf2","Ncf1","Ncf4")
Granules <- c("Mpo","Elane","Ctsg","Prtn3","Prss57","Ctsc","Camp","Ltf","Lcn2","Lyz2",
               "Cybb","Cyba","Cyb5r4","Hp","Ngp","Mmp1","Mmp8","Mmp9","Slpi","Itgam") 

#Calculate Module Score
Neutrophils.integrated <- AddModuleScore(Neutrophils.integrated, features = list(Aging),
                                         ctrl.size = 5, name = "Aging") 
Neutrophils.integrated <- AddModuleScore(Neutrophils.integrated, features = list(AntiApoptosis),
                                         ctrl.size = 5, name = "AntiApoptosis") 
Neutrophils.integrated <- AddModuleScore(Neutrophils.integrated, features = list(ProApoptosis),
                                         ctrl.size = 5, name = "ProApoptosis") 
Neutrophils.integrated <- AddModuleScore(Neutrophils.integrated, features = list(NADPH),
                                         ctrl.size = 5, name = "NADPHOxidase")

Neutrophils.integrated[['AgingModule']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophils.integrated,
                                                                                vars = 'Aging1')))
Neutrophils.integrated[['AntiApoptosisModule']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophils.integrated,
                                                                                    vars = 'AntiApoptosis1')))
Neutrophils.integrated[['ProApoptosisModule']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophils.integrated,
                                                                                    vars = 'ProApoptosis1')))
Neutrophils.integrated[['NADPHOxidaseModule']] <- CreateAssayObject(data = t(x = FetchData(object = Neutrophils.integrated,
                                                                                    vars = 'NADPHOxidase1')))
#Visualization with Violin Plots
VlnPlot(object = Neutrophils.integrated, features = "Aging1", 
        split.by = 'orig.ident')
VlnPlot(object = Neutrophils.integrated, features = "AntiApoptosis1", 
        split.by = 'orig.ident', idents = "N4") & theme(axis.title.x = element_blank(),
                                                        axis.text.x = element_blank()) & ggtitle("Anti-Apoptosis")
VlnPlot(object = Neutrophils.integrated, features = "ProApoptosis1", 
        split.by = 'orig.ident', idents = "N4") & theme(axis.title.x = element_blank(),
                                                        axis.text.x = element_blank()) & ggtitle("Pro-Apoptosis")
VlnPlot(object = Neutrophils.integrated, features = "NADPHOxidase1", 
        split.by = 'orig.ident', idents = "N4") & theme(axis.title.x = element_blank(),
                                                        axis.text.x = element_blank()) & ggtitle("NADPH Oxidase")
#Visualization by MultiBar Heatmaps
library(DoMultiBarHeatmap)
DoMultiBarHeatmap(object = Neutrophils.scaled, features = TFs,
                  group.by="Identity", additional.group.by="orig.ident") + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 20),legend.text = element_text(size = 15))
DoMultiBarHeatmap(object = Neutrophils.scaled, features = AntiApoptosis,
                  group.by="Identity", additional.group.by="orig.ident") + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 20),legend.text = element_text(size = 15))
DoMultiBarHeatmap(object = Neutrophils.scaled, features = AntiApoptosis,
                  group.by="Identity", additional.group.by="orig.ident") + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 20),legend.text = element_text(size = 15))
DoMultiBarHeatmap(object = Neutrophils.scaled, features = NADPH,
                  group.by="Identity", additional.group.by="orig.ident") + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 20),legend.text = element_text(size = 15))
DoMultiBarHeatmap(object = Neutrophils.scaled, features = Granules,
                  group.by="Identity", additional.group.by="orig.ident") + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", 
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 25, face = "italic"),
        legend.title = element_text(size = 20),legend.text = element_text(size = 15))

###COMPARISON SUB-CLUSTER N4-DEGs 14dpi vs Naive
N4 <- subset(Neutrophils.integrated, ident = "N4")
N4
Naivevs14dpi.deg.N4 <- FindMarkers(N4, ident.1 = "14dpi", ident.2 = "Naive", 
                                   verbose = TRUE, group.by="orig.ident", 
                                   logfc.threshold = log(1.5))
Naivevs14dpi.deg.N4
write.csv(Naivevs14dpi.deg.N4,"./Naivevs14dpi.deg.N4.csv")
#Violin Plots for genes belongs to N4-DEGs 14dpi vs Naive
VlnPlot(object = Neutrophils.integrated,
        features = c("Stat1","Irf7","Ifit3","Ifi204","Ifi47","Isg15","Ifitm3",
                     "Igtp","Gbp7","Junb","Jund","Il1rap"), split.by = 'orig.ident', 
        idents = "N4") & theme(axis.title.x = element_blank(), axis.text.x = element_blank())

### TRAJECTORIES ANALYSIS WITH MONOCLE
library(monocle3)
library(SeuratWrappers)
# Convert the Seurat object to a CellDataSet object
Neutrophils.cds <- as.cell_data_set(Neutrophils.integrated)
Neutrophils.cds <- cluster_cells(cds = Neutrophils.cds, reduction_method = "UMAP")
Neutrophils.cds <- learn_graph(Neutrophils.cds, use_partition = TRUE)
# Order cells
Neutrophils.cds <- order_cells(Neutrophils.cds, reduction_method = "UMAP")
# Plot trajectories colored by pseudotime
plot_cells(
  cds = Neutrophils.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)





