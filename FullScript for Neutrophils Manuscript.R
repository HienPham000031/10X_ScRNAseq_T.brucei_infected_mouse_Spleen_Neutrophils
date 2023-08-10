###T.brucei.Neutrophils - CREATED BY HIEN PHAM - EMAIL: ptthien3108@gmail.com
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
DimPlot(Neutrophils.integrated, reduction = "umap", group.by = "orig.ident",
        split.by = "orig.ident", cols = c("darkorange1","cyan4"))
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
DoMultiBarHeatmap(object = Neutrophils.scaled, features = ProApoptosis,
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
DoMultiBarHeatmap(object = Neutrophils.integrated, features = top10$gene,
                  group.by="Identity", additional.group.by="orig.ident") + 
  scale_fill_gradient2(low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white",
                       high = rev(c('#b2182b','#ef8a26','#fddbc7')),
                       midpoint = 0, guide = "colourbar",
                       aesthetics = "fill") + 
  theme(axis.text = element_text(size = 20, face = "italic"),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15))


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
DotPlot(Neutrophils.integrated, features = c("Timp1","Timp2"), 
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text = element_text(size = 15, face = "italic"),
                                                        legend.title = element_text(size = 20),
                                                        legend.text = element_text(size = 15),
                                                        axis.title.x = element_text(size = 25),
                                                        axis.title.y.left = element_text(size = 25))
VlnPlot(object = Neutrophils.integrated,
        features = c("Top2a","Pclaf","Zmpste24","Cebpe","Thbs1","Mmp8","Il1b",
                     "Ccl6")) & theme(axis.title.x = element_blank(), 
                                      axis.text.x = element_blank())

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


Full Splenocytes Annotation
### FULL ANNOTATION OF SPLENOCYTES
#This script operates on a previously saved Seurat object. 
#Please note that the number of Seurat clusters may vary if the script is executed 
#on a different computer or with an alternate version of R and Seurat
Idents(brucei.integrated) <- 'seurat_clusters'
##Visualizing Recognized Canonical Markers for Each Cell Type
#T cells markers
DotPlot(brucei.integrated, features = c("Trac","Trbc1","Cd3e","Cd3g","Cd4","Cd8a","Cd8b",
                                        "Foxp3","Il2ra","Ctla4","Il7r","Sell",
                                        "Ccr7","Ifng","Id2","Ccl5","Icos","Maf","Il21",
                                        "Tnfrsf18","Itga1","Itga2","Cd27","Klrk1","Ncr1",
                                        "Klrg1","Klrd1","Gzma","Gzmb","Cd44"), 
        cols = c("blue", "red"), dot.scale = 8)
#Macrophages, Monocytes, Erythroid and HSPC markers
DotPlot(brucei.integrated, features = c("Lst1","Ms4a6c","Plac8","Ifitm3",
                                        "Adgre1","C1qa","C1qb","C1qc","Fcgr1","Ly6c2",
                                        "Itgam","Ccr2","Kit","Gata1","Tal1","Gypa",
                                        "Hbb-bt","Tfrc","Cst3","Itgax","Cd8a","Clec9a","Klk1",
                                        "Siglech","Pacsin1","Pntt","Irf8","Ccr9","Ly6d",
                                        "Ly6i","Cebpb","Cybb"), 
        cols = c("blue", "red"), dot.scale = 8)
#DCs and pDCs markers
DotPlot(brucei.integrated, features = c("Cst3","Itgax","Cd8a","Clec9a","Klk1",
                                        "Siglech","Pacsin1","Pntt","Irf8","Ccr9","Ly6d"),
        cols = c("blue", "red"), dot.scale = 8)
#Annotate OF 26-27
DotPlot(brucei.integrated, features = c("Acta2","Vim","Fap","Pdgfra","Pdgfrb","S100a4",
                                        "Fn1","Col1a1","Postn","Dcn","Col1a2","Fbln2",
                                        "Des","Cdm11","Cd73","Cd90","Cd105"),
        cols = c("blue", "red"), dot.scale = 8)
DotPlot(brucei.integrated, features = c("Kit","Ly6a","Fcgr3","Fcgr2","Cd34","Fcer1a"),
        cols = c("blue", "red"), dot.scale = 8)
#B cells markers
DotPlot(brucei.integrated, features = c("Cd19","Cd79a","Cd79b","Pax5","Ly6d","Ighd",
                                        "Cd55","Fcer2a","Cd93","Cd24a","Ms4a1",
                                        "Vpreb3","Cr2","Cd1d1","Ebf1","Apoe",
                                        "Plac8","Zbtb20","Zbtb32","Ahnak","Vim",
                                        "S100a6","Fcrl5","Ly86","Bhlhe41","Aicda",
                                        "Mki67","Cdk1","Cdc20","Pcna","Pclaf","Ube2c",
                                        "Sdc1","Xbp1","Prdm1","Irf4","Jchain",
                                        "Iglc1","Ighm","Ighg1","Ighg2b","Ighg2c","Ighg3","Cd5"), 
        cols = c("blue", "red"), dot.scale = 8)

##RENAME and COLOR by CELLTYPES
new.cluster.ids <- c("FoB","Macrophage","Naive CD4 T","CD4 Tfh","Naive CD8 T","Neutrophils",
                     "MZB","Erythroid","CD8 Eff T","cDC","GC-like","Unknown B cells",
                     "Plasma cells","Erythroid","NK/NKT","cDC","HPSC","Neutrophils",
                     "TransB","Th1","Neutrophils","Tregs","Macrophage","Macrophage",
                     "Monocyte"," ","Monocyte","pDC")
names(new.cluster.ids) <- levels(brucei.integrated)
brucei.annotated <- RenameIdents(brucei.integrated, new.cluster.ids)

###SUBSET B CELLS (WITHOUT PLASMA CELLS)
#Due to the presence of an 'unknown B cells' population from previous steps, 
#which is challenging to annotate as a specific B cell subtype, 
#we will subset B cells for an additional round of clustering
B.cells.subset <- subset((brucei.annotated), idents = c("Unknown B cells","TransB","MZB",
                                                        "FoB","GC-like"))
DefaultAssay(B.cells.subset) <- "integrated"
B.cells.subset <- ScaleData(B.cells.subset,vars.to.regress = "CC.Difference", verbose = FALSE)

# Perform linear dimensional reduction and Run non-linear dimensional reduction (UMAP)
B.cells.subset <- RunPCA(B.cells.subset, npcs = 50, verbose = FALSE)
DimHeatmap(B.cells.subset, dims = 1:20, cells = 500, balanced = TRUE)
B.cells.subset <- RunUMAP(B.cells.subset, reduction = "pca", dims = 1:20)
B.cells.subset <- FindNeighbors(B.cells.subset, reduction = "pca", dims = 1:20)
B.cells.subset <- FindClusters(B.cells.subset, resolution = 0.5)

# Visualization
DimPlot(B.cells.subset, reduction = "umap", group.by = "orig.ident")
DimPlot(B.cells.subset, reduction = "umap", split.by = "orig.ident", label = TRUE, 
        repel = TRUE)
DimPlot(B.cells.subset, reduction = "umap")
DotPlot(B.cells.subset, features = c("Cd19","Cd79a","Cd79b","Pax5","Ly6d","Ighd",
                                     "Cd55","Fcer2a","Cd93","Cd24a","Ms4a1",
                                     "Vpreb3","Cr2","Cd1d1","Ebf1","Apoe",
                                     "Plac8","Zbtb20","Zbtb32","Ahnak","Vim",
                                     "S100a6","Fcrl5","Ly86","Bhlhe41","Aicda",
                                     "Mki67","Cdk1","Cdc20","Pcna","Pclaf","Ube2c",
                                     "Sdc1","Xbp1","Prdm1","Irf4","Jchain",
                                     "Iglc1","Ighm","Ighg1","Ighg2b","Ighg2c","Ighg3","Cd5"), 
        cols = c("blue", "red"), dot.scale = 8)
new.cluster.ids.B <- c("FoB","FoB","FoB","GC-like","TransB","FoB","GC-like","B1","MZB",
                       "FoB","FoB","FoB","FoB","CD4 Tfh")
names(new.cluster.ids.B) <- levels(B.cells.subset)
B.cells.subset.annotated <- RenameIdents(B.cells.subset, new.cluster.ids.B)
saveRDS(B.cells.subset.annotated,"./B.cell.subset.annotated")

##Copying the label of subset back to main Seurat object
Idents(brucei.annotated)
brucei.annotated[["annotation"]] <- Idents(brucei.annotated);
B.cells.subset.annotated[["annotation"]] <- Idents(B.cells.subset.annotated)
head(brucei.annotated)
rownames(brucei.annotated@meta.data)
brucei.annotated@meta.data$annotation
B.cells.subset.annotated@meta.data
brucei.annotated@meta.data$annotation
head(brucei.annotated)

ind = match(rownames(B.cells.subset.annotated@meta.data), rownames(brucei.annotated@meta.data))
ind
brucei.annotated@meta.data$annotation <- as.character(brucei.annotated@meta.data$annotation)
brucei.annotated@meta.data[ind, "annotation"] = as.character.factor(B.cells.subset.annotated@meta.data[ ,"annotation"])
brucei.annotated$annotation
brucei.annotated@meta.data$annotation <- as.factor(brucei.annotated@meta.data$annotation)
B.cells.subset.annotated@meta.data$annotation
Idents(brucei.annotated) <- brucei.annotated@meta.data$annotation
table(Idents(B.cells.subset.annotated), B.cells.subset.annotated$orig.ident)
DimPlot(brucei.annotated, reduction = "umap", split.by = "orig.ident", label = TRUE, 
        repel = TRUE)
levels(brucei.annotated) <- c("Plasma cells","GC-like","FoB","MZB","B1","TransB",
                              "NK/NKT","CD8 Eff T","Naive CD8 T","Tregs","Th1","CD4 Tfh",
                              "Naive CD4 T","Macrophage","Monocyte","HPSC","Erythroid","cDC",
                              "pDC","Neutrophils"," ")
saveRDS(brucei.annotated,"./brucei.annotated.rds")

##Copying also label N1 N2 N3 N4 to the big brucei object
neutrophils[["annotation"]] <- Idents(neutrophils)
head(neutrophils)
brucei.annotated2<- brucei.annotated
ind.neutrophils = match(rownames(neutrophils@meta.data), rownames(brucei.annotated2@meta.data))
ind.neutrophils
brucei.annotated2@meta.data$annotation <- as.character(brucei.annotated2@meta.data$annotation)
brucei.annotated2@meta.data[ind.neutrophils, "annotation"] = as.character.factor(neutrophils@meta.data[ ,"annotation"])
brucei.annotated2$annotation
Idents(brucei.annotated2) <- brucei.annotated2@meta.data$annotation
table(Idents(brucei.annotated2), brucei.annotated2$orig.ident)
table(Idents(neutrophils), neutrophils$orig.ident)
brucei.annotated2@meta.data$annotation <- as.factor(brucei.annotated2@meta.data$annotation)
DimPlot(brucei.annotated2, reduction = "umap", split.by = "orig.ident", label = TRUE, 
        repel = TRUE)
levels(brucei.annotated2)
new.cluster.ids.full2 <- c("Macrophage","Naive CD8 T","Erythroid","FoB","Naive CD4 T",
                           "CD4 Tfh","NK/NKT","MZB","B1"," ","Tregs","CD8 Eff T","cDC",
                           "Th1","N4","HPSC","GC-like","Plasma cells","TransB","Monocyte",
                           "pDC","N1","GMP","N2","N3")
names(new.cluster.ids.full2) <- levels(brucei.annotated2)
brucei.annotated2 <- RenameIdents(brucei.annotated2, new.cluster.ids.full2)
levels(brucei.annotated2) <- c("Plasma cells","GC-like","FoB","MZB","B1","TransB",
                               "NK/NKT","CD8 Eff T","Naive CD8 T","Tregs","Th1","CD4 Tfh",
                               "Naive CD4 T","Macrophage","Monocyte","HPSC","Erythroid","cDC",
                               "pDC","GMP","N1","N2","N3","N4"," ")
saveRDS(brucei.annotated2,"./brucei.annotated2.rds")

###VISUALIZATION OF GENES FOR ANNOTATION OF CELLTYPES
DefaultAssay(brucei.annotated2) <- "RNA"
#B cells
DotPlot(brucei.annotated2, idents = c("Plasma cells","GC-like","FoB","MZB","B1","TransB"),
        features = c("Cd19","Cd79a","Cd79b","Pax5","Ly6d","Ighd",
                     "Cd55","Fcer2a","Cd93","Cd24a","Ms4a1",
                     "Vpreb3","Cr2","Cd1d1","Ebf1","Apoe",
                     "Plac8","Zbtb20","Zbtb32","Itgam","Itgax",
                     "Tbx21","Ahnak","Vim",
                     "S100a6","Fcrl5","Ly86","Bhlhe41","Aicda",
                     "Mki67","Cdk1","Cdc20","Pcna","Pclaf","Ube2c",
                     "Sdc1","Xbp1","Prdm1","Irf4","Jchain",
                     "Iglc1","Ighm","Ighg1","Ighg2b","Ighg2c","Ighg3","Cd5"), 
        cols = c("blue", "red"), dot.scale = 8) + NoLegend() + theme(axis.text = element_text(size = 15, face = "italic"),
                                                                     axis.text.x = element_text(angle=45, hjust = 1),
                                                                     legend.title = element_text(size = 20),
                                                                     legend.text = element_text(size = 15),
                                                                     axis.title.x = element_text(size = 25),
                                                                     axis.title.y.left = element_text(size = 25))
#T cells
DotPlot(brucei.annotated2,idents = c("NK/NKT","CD8 Eff T","Naive CD8 T","Tregs","Th1","CD4 Tfh","Naive CD4 T"),
        features = c("Trac","Trbc1","Cd3e","Cd3g","Cd4","Cd8a","Cd8b",
                     "Foxp3","Il2ra","Ctla4","Il7r","Sell",
                     "Ccr7","Cd44","Ifng","Id2","Ccl5","Icos","Maf","Il21",
                     "Cxcr5","Cxcr3","Tbx21","Bcl6",
                     "Tnfrsf18","Itga1","Itga2","Klrk1","Ncr1",
                     "Klrg1","Klrd1","Nkg2d","Klrb1c","Gzma","Gzmb"), 
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text = element_text(size = 15, face = "italic"),
                                                        axis.text.x = element_text(angle=45, hjust = 1),
                                                        legend.title = element_text(size = 20),
                                                        legend.text = element_text(size = 15),
                                                        axis.title.x = element_text(size = 25),
                                                        axis.title.y.left = element_text(size = 25))
DotPlot(brucei.annotated, idents = c("Macrophage","Monocyte","HPSC","Erythroid","cDC","pDC"),
        features = c("Adgre1","C1qa","C1qb","C1qc","Fcgr1","Ly6c2","Plac8",
                     "Itgam","Ccr2","Kit","Gata1","Tal1","Gypa",
                     "Hbb-bt","Tfrc","Cst3","Itgax","Trbc1","Cadm1","Xcr1",
                     "Id2","Batf3",
                     "Cd8a","Clec9a","Klk1",
                     "Siglech","Pacsin1","Pntt","Irf8",
                     "Irf7","Tlr7","Ccr9","Ly6d"), 
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text = element_text(size = 15, face = "italic"),
                                                        axis.text.x = element_text(angle=45, hjust = 1),
                                                        legend.title = element_text(size = 20),
                                                        legend.text = element_text(size = 15),
                                                        axis.title.x = element_text(size = 25),
                                                        axis.title.y.left = element_text(size = 25))
DotPlot(brucei.annotated2, idents = c("GMP","N1","N2","N3","N4"),
        features = c("Cd34","Cd63","Ms4a3","Mpo","Elane","Irf8","Ly86","Gfi1",
                     "Vcam1","Rgcc","Cd81","Fcnb","Cebpe","Cd177","S100a8","Ly6g",
                     "Ly6c2","Itgam","Mmp8","Mmp9","Csf1r","Cxcr4","Ccl6","Cxcr2","Il1b"), 
        cols = c("blue", "red"), dot.scale = 8) + theme(axis.text = element_text(size = 15, face = "italic"),
                                                        axis.text.x = element_text(angle=45, hjust = 1),
                                                        legend.title = element_text(size = 20),
                                                        legend.text = element_text(size = 15),
                                                        axis.title.x = element_text(size = 25),
                                                        axis.title.y.left = element_text(size = 25))   
DimPlot(object = brucei.annotated2, reduction = 'umap',label = TRUE, split.by = "orig.ident",
        cols = c('Macrophage'='#F8766D','Naive CD8 T'='#EF7F46','Erythroid'='#E38900',
                 'FoB'='#D59100','Naive CD4 T'='#C49A00','CD4 Tfh'='#B0A100',
                 'NK/NKT'='#99A800','MZB'='#53B400','B1'='#00B823','Tregs'='#00BC56',
                 'CD8 Eff T'='#00BF77','cDC'='#00C094','Th1'='#00C1AD','GMP'='#D59100',
                 'N1'='#F166E9','N2'='#FB61D7','N3'='#FF62C1','N4'='#FF66A8','HPSC'='#06A4FF',
                 'GC-like'='#06A4FF','Plasma cells'='#7698FF','TransB'='#A58AFF',
                 'Monocyte'='#00BFC4','pDC'='#00B6EB',' '='#D59100')) + NoLegend()

###VISUALIZATION EXPRESSION OF CERTAIN GENES RELEVANT TO THE MANUSCRIPT
FeaturePlot(brucei.annotated2,features = c("Mmp8","Mmp9","Mmp25"), cols = c("grey","red"), keep.scale = "all")
VlnPlot(brucei.annotated2, idents = c("N1","N2","N3","N4"),
        features = "Mmp8", split.by = "orig.ident", y.max = 5)
VlnPlot(brucei.annotated2, idents = c("N1","N2","N3","N4"),
        features = "Mmp9", split.by = "orig.ident", y.max = 5)
VlnPlot(brucei.annotated2, idents = c("N1","N2","N3","N4"),
        features = "Mmp25", split.by = "orig.ident", y.max = 5)
FeaturePlot(brucei.annotated2, features = c("Tnfrsf13b"), split.by = "orig.ident", cols = c("grey","red"), keep.scale = "feature")                                                        

Nichenet PCs
###NICHENET FOR CELL-CELL INTERATION ANALYSIS
library(nichenetr)
library(circlize)
library(tidyverse)
#Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks:
ligand_target_matrix = readRDS("./ligand_target_matrix.rds")
ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns

head(ligand_target_matrix)
dim(ligand_target_matrix)

lr_network = readRDS("./lr_network.rds")
head(lr_network)
dim(lr_network)

weighted_networks = readRDS("./weighted_networks.rds")
head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
weighted_networks_lr


#Because the expression data is of mouse origin, we will convert the NicheNet 
#network gene symbols from human to mouse based on one-to-one orthology:
lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()

ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]

weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
weighted_networks_lr

###INFER THE INTERATION BETWEEN NEUTROPHILS (SENDER) AND PLASMA CELLS (RECEIVER)
##Define a “sender/niche” cell population and a “receiver/target” cell population present 
##in your expression data and determine which genes are expressed in both populations
## receiver
receiver = "Plasma cells"
expressed_genes_receiver = get_expressed_genes(receiver, brucei.annotated2, pct = 0.10, assay_oi = "RNA")
expressed_genes_receiver
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
background_expressed_genes
## sender
sender_celltypes = c("N1","N2","N3","N4")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Neutrophils.integrated, pct = 0.10, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
expressed_genes_sender

##Define a gene set of interest: these are the genes in the “receiver/target” cell population 
##that are potentially affected by ligands expressed by interacting cells 
##(e.g. genes differentially expressed upon cell-cell interaction)
brucei_receiver= subset(brucei.annotated2, idents = receiver)
brucei_receiver
brucei_receiver = SetIdent(brucei_receiver, value = brucei_receiver[["orig.ident"]])

condition_oi = "14dpi"
condition_reference = "Naive" 

DE_table_receiver = FindMarkers(object = brucei_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
DE_table_receiver
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.58) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
geneset_oi

##Define a set of potential ligands: these are ligands that are expressed by the
##“sender/niche” cell population and bind a (putative) receptor expressed by the 
##“receiver/target” population
ligands = lr_network %>% pull(from) %>% unique()
ligands
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expressed_ligands
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands

##Perform NicheNet ligand activity analysis: rank the potential ligands based on 
##the presence of their target genes in the gene set of interest 
##(compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

#The number of top-ranked ligands that are further used to predict active target 
#genes and construct an active ligand-receptor network is here 20.
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
best_upstream_ligands <- best_upstream_ligands[! best_upstream_ligands %in% c('Ptprc','Csf1')] ##Remove ligands which receptors are irrelevant to plasma cells
DotPlot(Neutrophils.integrated, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

###Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() 
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

###Circos plots to visualize ligand-target and ligand-receptor interactions
##Calculate average receptors expression in sender cells
# avg_expression_ligands = AverageExpression(seuratObj %>% subset(subset = aggregate == "LCMV"),features = nichenet_output$top_ligands) # if want to look specifically in LCMV-only cells
DefaultAssay(neutrophils) <- "RNA"
avg_expression_ligands = AverageExpression(Neutrophils.integrated, features = best_upstream_ligands, assays = "RNA")
avg_expression_ligands

#Assign ligands to sender cells
#To assign ligands to sender cell type, we can e.g. look for which sender cell types 
#show an expression that is higher than the average + SD.
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

N1_specific_ligands = sender_ligand_assignment$N1 %>% names() %>% setdiff(general_ligands)
N2_specific_ligands = sender_ligand_assignment$N2 %>% names() %>% setdiff(general_ligands)
N3_specific_ligands = sender_ligand_assignment$N3 %>% names() %>% setdiff(general_ligands)
N4_specific_ligands = sender_ligand_assignment$N4 %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("N1-specific", times = N1_specific_ligands %>% length()),
                  rep("N2-specific", times = N2_specific_ligands %>% length()),
                  rep("N3-specific", times = N3_specific_ligands %>% length()),
                  rep("N4-specific", times = N4_specific_ligands %>% length()),
                  rep("General", times = general_ligands %>% length())),
  ligand = c(N1_specific_ligands, N2_specific_ligands, N3_specific_ligands,N4_specific_ligands,general_ligands))

##Define the ligand-target links of interest
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Plasma Cells") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

###Prepare the circos visualization: give each segment of ligands and targets a specific color and order
grid_col_ligand =c("General" = "lawngreen",
                   "N1-specific" = "royalblue",
                   "N2-specific" = "darkgreen",
                   "N3-specific" = "violet",
                   "N4-specific" = "steelblue2")
grid_col_target =c(
  "Plasma Cells" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target,weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

###Prepare the circos visualization: order ligands and targets
target_order = circos_links$target %>% unique()
ligand_order = c(N1_specific_ligands, N2_specific_ligands, N3_specific_ligands, N4_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


###Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N2-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N3-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N4-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Plasma Cells") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)
table(circos_links)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

circos.clear()

#Render the circos plot (degree of transparancy determined by the regulatory 
#potential value of a ligand-target interaction)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

circos.clear()
#Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
weighted_networks_lr
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
lr_network_top_df

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_receptors
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network
colnames(vis_ligand_receptor_network)
###Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor 
###interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

###Visualize ligand-receptor interactions of the prioritized ligands in a circos plot
lr_network_top_df = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)
lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "Plasma cells_receptors") %>% inner_join(ligand_type_indication_df)
lr_network_top_df

grid_col_ligand =c("General" = "lawngreen",
                   "N1-specific" = "royalblue",
                   "N2-specific" = "darkgreen",
                   "N3-specific" = "violet",
                   "N4-specific" = "steelblue2")
grid_col_receptor =c(
  "Plasma cells_receptors" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# Prepare the circos visualization: order ligands and receptors
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(N1_specific_ligands, N2_specific_ligands, N3_specific_ligands,N4_specific_ligands, general_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)

#Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N2-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N3-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N4-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "Plasma cells_receptors") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)
circos.clear()
#Render the circos plot (all links same transparancy). Only the widths of the blocks 
#that indicate each receptor is proportional the ligand-receptor interaction weight 
#(~prior knowledge supporting the interaction).
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()

#Render the circos plot (degree of transparancy determined by the prior interaction 
#weight of the ligand-receptor interaction - just as the widths of the blocks indicating each receptor)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()

### VERIFICATION THE DEGS PER CELLTYPES
library(EnhancedVolcano)
Plasmacells<- subset(brucei.annotated2, idents = "Plasma cells")
#Finding DEGs of 14dpi compares to Naive
Naivevs14dpi.deg.Plasmacells <- FindMarkers(Plasmacells, ident.1 = "14dpi", ident.2 = "Naive", verbose = TRUE, group.by="orig.ident", logfc.threshold = log(1.5))
VolcanoPlasma_Nichenet_R <- subset(Naivevs14dpi.deg.Plasmacells, 
                                   rownames(Naivevs14dpi.deg.Plasmacells) %in% receptor_order)
VolcanoPlasma_Nichenet_R
EnhancedVolcano(VolcanoPlasma_Nichenet_R,
                lab = as.character(rownames(VolcanoPlasma_Nichenet_R)),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                drawConnectors = TRUE,
                widthConnectors = 0.75)
FeaturePlot(brucei.annotated2, features = c("Tnfsf13b"), cols = c("grey","red"),
            split.by = "orig.ident")
VlnPlot(object = brucei.annotated2,idents = "Plasma cells", 
        features = c("Tnfrsf13b","Tnfrsf17"), split.by = "orig.ident") + theme(plot.title = element_text(size = 20, face = "bold"),
                                                                               axis.text.x = element_text(angle=0, hjust = 0.5),
                                                                               legend.title=element_text(size=20), 
                                                                               legend.text=element_text(size=20),
                                                                               axis.text = element_text(size = 20),
                                                                               axis.title.y.left = element_text(size = 25))



Nichenet FoB
###NICHENETR FOR Fo B CELLS
##Define a “sender/niche” cell population and a “receiver/target” cell population present 
##in your expression data and determine which genes are expressed in both populations
## receiver
receiver = "FoB"
expressed_genes_receiver = get_expressed_genes(receiver, brucei.annotated2, pct = 0.10, assay_oi = "RNA")
expressed_genes_receiver
background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
background_expressed_genes
## sender
sender_celltypes = c("N1","N2","N3","N4")

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, Neutrophils.integrated, pct = 0.10, assay_oi = "RNA") # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()
expressed_genes_sender

##Define a gene set of interest: these are the genes in the “receiver/target” cell population 
##that are potentially affected by ligands expressed by interacting cells 
##(e.g. genes differentially expressed upon cell-cell interaction)
brucei_receiver= subset(brucei.annotated2, idents = receiver)
brucei_receiver
brucei_receiver = SetIdent(brucei_receiver, value = brucei_receiver[["orig.ident"]])

condition_oi = "14dpi"
condition_reference = "Naive" 

DE_table_receiver = FindMarkers(object = brucei_receiver, ident.1 = condition_oi, ident.2 = condition_reference, min.pct = 0.10) %>% rownames_to_column("gene")
DE_table_receiver
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.58) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
geneset_oi

##Define a set of potential ligands: these are ligands that are expressed by the
##“sender/niche” cell population and bind a (putative) receptor expressed by the 
##“receiver/target” population
ligands = lr_network %>% pull(from) %>% unique()
ligands
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)
expressed_ligands
potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()
potential_ligands

##Perform NicheNet ligand activity analysis: rank the potential ligands based on 
##the presence of their target genes in the gene set of interest 
##(compared to the background set of genes)
ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))
ligand_activities

#The number of top-ranked ligands that are further used to predict active target 
#genes and construct an active ligand-receptor network is here 20.
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()
DotPlot(Neutrophils.integrated, features = best_upstream_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()

###Infer receptors and top-predicted target genes of ligands that are top-ranked in the ligand activity analysis
active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 200) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

###Circos plots to visualize ligand-target and ligand-receptor interactions
##Calculate average receptors expression in sender cells
# avg_expression_ligands = AverageExpression(seuratObj %>% subset(subset = aggregate == "LCMV"),features = nichenet_output$top_ligands) # if want to look specifically in LCMV-only cells
DefaultAssay(neutrophils) <- "RNA"
avg_expression_ligands = AverageExpression(Neutrophils.integrated, features = best_upstream_ligands, assays = "RNA")
avg_expression_ligands

#Assign ligands to sender cells
#To assign ligands to sender cell type, we can e.g. look for which sender cell types 
#show an expression that is higher than the average + SD.
sender_ligand_assignment = avg_expression_ligands$RNA %>% apply(1, function(ligand_expression){
  ligand_expression > (ligand_expression %>% mean() + ligand_expression %>% sd())
}) %>% t()
sender_ligand_assignment = sender_ligand_assignment %>% apply(2, function(x){x[x == TRUE]}) %>% purrr::keep(function(x){length(x) > 0})
names(sender_ligand_assignment)

all_assigned_ligands = sender_ligand_assignment %>% lapply(function(x){names(x)}) %>% unlist()
unique_ligands = all_assigned_ligands %>% table() %>% .[. == 1] %>% names()
general_ligands = best_upstream_ligands %>% setdiff(unique_ligands)

N1_specific_ligands = sender_ligand_assignment$N1 %>% names() %>% setdiff(general_ligands)
N2_specific_ligands = sender_ligand_assignment$N2 %>% names() %>% setdiff(general_ligands)
N3_specific_ligands = sender_ligand_assignment$N3 %>% names() %>% setdiff(general_ligands)
N4_specific_ligands = sender_ligand_assignment$N4 %>% names() %>% setdiff(general_ligands)

ligand_type_indication_df = tibble(
  ligand_type = c(rep("N1-specific", times = N1_specific_ligands %>% length()),
                  rep("N2-specific", times = N2_specific_ligands %>% length()),
                  rep("N3-specific", times = N3_specific_ligands %>% length()),
                  rep("N4-specific", times = N4_specific_ligands %>% length())),
  ligand = c(N1_specific_ligands, N2_specific_ligands, N3_specific_ligands,N4_specific_ligands))

##Define the ligand-target links of interest
active_ligand_target_links_df = active_ligand_target_links_df %>% mutate(target_type = "Fo B Cells") %>% inner_join(ligand_type_indication_df) # if you want ot make circos plots for multiple gene sets, combine the different data frames and differentiate which target belongs to which gene set via the target type
cutoff_include_all_ligands = active_ligand_target_links_df$weight %>% quantile(0.40)

active_ligand_target_links_df_circos = active_ligand_target_links_df %>% filter(weight > cutoff_include_all_ligands)

ligands_to_remove = setdiff(active_ligand_target_links_df$ligand %>% unique(), active_ligand_target_links_df_circos$ligand %>% unique())
targets_to_remove = setdiff(active_ligand_target_links_df$target %>% unique(), active_ligand_target_links_df_circos$target %>% unique())

circos_links = active_ligand_target_links_df %>% filter(!target %in% targets_to_remove &!ligand %in% ligands_to_remove)

###Prepare the circos visualization: give each segment of ligands and targets a specific color and order
grid_col_ligand =c("N1-specific" = "royalblue",
                   "N2-specific" = "darkgreen",
                   "N3-specific" = "violet",
                   "N4-specific" = "steelblue2")
grid_col_target =c(
  "Fo B Cells" = "tomato")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_target = tibble(target_type = grid_col_target %>% names(), color_target_type = grid_col_target)

circos_links = circos_links %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as target!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_target)
links_circle = circos_links %>% select(ligand,target,weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
target_color = circos_links %>% distinct(target,color_target_type)
grid_target_color = target_color$color_target_type %>% set_names(target_color$target)

grid_col =c(grid_ligand_color,grid_target_color)

# give the option that links in the circos plot will be transparant ~ ligand-target potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

###Prepare the circos visualization: order ligands and targets
target_order = circos_links$target %>% unique()
###To remove 'general'
ligand_order = c(N1_specific_ligands, N2_specific_ligands, N3_specific_ligands, N4_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,target_order)


###Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_target = 15
width_same_cell_same_target_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N2-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N3-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N4-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_target_type, times = (circos_links %>% filter(target_type == "Fo B Cells") %>% distinct(target) %>% nrow() -1)),
  width_ligand_target
)
table(circos_links)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

circos.clear()

#Render the circos plot (degree of transparancy determined by the regulatory 
#potential value of a ligand-target interaction)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 1)
}, bg.border = NA) #

#Receptors of top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)
weighted_networks_lr
lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)
lr_network_top_df

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_receptors
order_receptors <- order_receptors[! order_receptors %in% c('Cd2')] ##Remove receptors irrelevant to B cells
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network

###Receptors of top-ranked ligands, but after considering only bona fide ligand-receptor 
###interactions documented in literature and publicly available databases
lr_network_strict = lr_network %>% filter(database != "ppi_prediction_go" & database != "ppi_prediction")
ligands_bona_fide = lr_network_strict %>% pull(from) %>% unique()
receptors_bona_fide = lr_network_strict %>% pull(to) %>% unique()

lr_network_top_df_large_strict = lr_network_top_df_large %>% distinct(from,to) %>% inner_join(lr_network_strict, by = c("from","to")) %>% distinct(from,to)
lr_network_top_df_large_strict = lr_network_top_df_large_strict %>% inner_join(lr_network_top_df_large, by = c("from","to"))

lr_network_top_df_strict = lr_network_top_df_large_strict %>% spread("from","weight",fill = 0)
lr_network_top_matrix_strict = lr_network_top_df_strict %>% select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df_strict$to)

dist_receptors = dist(lr_network_top_matrix_strict, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix_strict %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix_strict))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix_strict))

vis_ligand_receptor_network_strict = lr_network_top_matrix_strict[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network_strict) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network_strict) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network_strict = vis_ligand_receptor_network_strict %>% t() %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential\n(bona fide)")
p_ligand_receptor_network_strict

###Visualize ligand-receptor interactions of the prioritized ligands in a circos plot
lr_network_top_df = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors) %>% rename(ligand = from, receptor = to)
lr_network_top_df = lr_network_top_df %>% mutate(receptor_type = "Fo B cells_receptors") %>% inner_join(ligand_type_indication_df)
lr_network_top_df

grid_col_ligand =c("N1-specific" = "royalblue",
                   "N2-specific" = "darkgreen",
                   "N3-specific" = "violet",
                   "N4-specific" = "steelblue2")
grid_col_receptor =c(
  "Fo B cells_receptors" = "darkred")

grid_col_tbl_ligand = tibble(ligand_type = grid_col_ligand %>% names(), color_ligand_type = grid_col_ligand)
grid_col_tbl_receptor = tibble(receptor_type = grid_col_receptor %>% names(), color_receptor_type = grid_col_receptor)

circos_links = lr_network_top_df %>% mutate(ligand = paste(ligand," ")) # extra space: make a difference between a gene as ligand and a gene as receptor!
circos_links = circos_links %>% inner_join(grid_col_tbl_ligand) %>% inner_join(grid_col_tbl_receptor)
links_circle = circos_links %>% select(ligand,receptor, weight)

ligand_color = circos_links %>% distinct(ligand,color_ligand_type)
grid_ligand_color = ligand_color$color_ligand_type %>% set_names(ligand_color$ligand)
receptor_color = circos_links %>% distinct(receptor,color_receptor_type)
grid_receptor_color = receptor_color$color_receptor_type %>% set_names(receptor_color$receptor)

grid_col =c(grid_ligand_color,grid_receptor_color)

# give the option that links in the circos plot will be transparant ~ ligand-receptor potential score
transparency = circos_links %>% mutate(weight =(weight-min(weight))/(max(weight)-min(weight))) %>% mutate(transparency = 1-weight) %>% .$transparency 

# Prepare the circos visualization: order ligands and receptors
receptor_order = circos_links$receptor %>% unique()
ligand_order = c(N1_specific_ligands, N2_specific_ligands, N3_specific_ligands,N4_specific_ligands) %>% c(paste(.," ")) %>% intersect(circos_links$ligand)
order = c(ligand_order,receptor_order)

#Prepare the circos visualization: define the gaps between the different segments
width_same_cell_same_ligand_type = 0.5
width_different_cell = 6
width_ligand_receptor = 15
width_same_cell_same_receptor_type = 0.5

gaps = c(
  # width_ligand_target,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N1-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N2-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N3-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_ligand_type, times = (circos_links %>% filter(ligand_type == "N4-specific") %>% distinct(ligand) %>% nrow() -1)),
  width_different_cell,
  rep(width_same_cell_same_receptor_type, times = (circos_links %>% filter(receptor_type == "Fo B cells_receptors") %>% distinct(receptor) %>% nrow() -1)),
  width_ligand_receptor
)
circos.clear()
#Render the circos plot (all links same transparancy). Only the widths of the blocks 
#that indicate each receptor is proportional the ligand-receptor interaction weight 
#(~prior knowledge supporting the interaction).
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = 0, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()

#Render the circos plot (degree of transparancy determined by the prior interaction 
#weight of the ligand-receptor interaction - just as the widths of the blocks indicating each receptor)
circos.par(gap.degree = gaps)
chordDiagram(links_circle, directional = 1,order=order,link.sort = TRUE, link.decreasing = FALSE, grid.col = grid_col,transparency = transparency, diffHeight = 0.005, direction.type = c("diffHeight", "arrows"),link.arr.type = "big.arrow", link.visible = links_circle$weight >= cutoff_include_all_ligands,annotationTrack = "grid", 
             preAllocateTracks = list(track.height = 0.075))
# we go back to the first track and customize sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.55), cex = 0.8)
}, bg.border = NA) #
circos.clear()

### VERIFICATION THE DEGS PER CELLTYPES
FoBcells<- subset(brucei.annotated2, idents = "FoB")
#Finding DEGs of 14dpi compares to Naive
Naivevs14dpi.deg.FoBcells <- FindMarkers(FoBcells, ident.1 = "14dpi", ident.2 = "Naive", verbose = TRUE, group.by="orig.ident", logfc.threshold = log(1.5))
head(Naivevs14dpi.deg.FoBcells)
Naivevs14dpi.deg.FoBcells
VolcanoFoB_Nichenet_R <- subset(Naivevs14dpi.deg.FoBcells, 
                                rownames(Naivevs14dpi.deg.FoBcells) %in% receptor_order)
EnhancedVolcano(VolcanoFoB_Nichenet_R,
                lab = as.character(rownames(VolcanoFoB_Nichenet_R)),
                x = 'avg_log2FC',
                y = 'p_val',
                pCutoff = 0.05,
                FCcutoff = 0.1,
                drawConnectors = TRUE,
                widthConnectors = 0.75)



