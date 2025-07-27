suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
  library(ggplot2)
  library(ggsci)
  library(dplyr)
  library(psych)
  library(pheatmap)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(DOSE)
  library(GOSemSim)
  library(enrichplot)
  library(patchwork)
  library(gridExtra)
  library(data.table)
  library(tidyverse)
  library(AnnoProbe)
  library(VENcelltypes)
  library(feather)
  library(matrixStats)
  library(edgeR)
  library(GenomicFeatures)
  library(parallel)
  library(COSG)
  library(WGCNA)
  library(hdWGCNA)
  library(scPred)
  library(magrittr)
  library(plot1cell)
  library(EnhancedVolcano)
  library(viridis)
  library(msigdbr)
  library(GSVA)
  library(SCENIC)
  library(report)
  library(enrichR)
  library(igraph)
  library(JASPAR2020)
  library(JASPAR2024)
  library(motifmatchr)
  library(TFBSTools)
  library(EnsDb.Hsapiens.v86)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(GenomicRanges)
  library(xgboost)
})

rm(list = ls())
options(stringsAsFactors = FALSE)



# loading reference (2018 nature MTG)
MTG18_exon <- fread("reference/Human MTG 2018/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_exon-matrix.csv", header = T, data.table = F)
MTG18_exon <- column_to_rownames(MTG18_exon, var = "V1")

MTG18_gene_info <- read.csv(file = "reference/Human MTG 2018/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_genes-rows.csv", header = TRUE)
table(duplicated(MTG18_gene_info$gene))
table(duplicated(MTG18_gene_info$entrez_id))


MTG18_meta.data <- read.csv(file = "reference/Human MTG 2018/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_samples-columns.csv", header = TRUE)
table(MTG18_meta.data$class)
## GABAergic Glutamatergic      no class  Non-neuronal 
## 4164         10525           325           914

txdb <- makeTxDbFromGFF(file = "reference/Human MTG 2018/rsem_GRCh38.p2.gtf/rsem_GRCh38.p2.gtf", format = "gtf")
exons.list.per.gene <- exonsBy(txdb, by = "gene")

gene_len_from_exon <- sum(width(exons.list.per.gene))
gene_len_from_exon <- data.frame(ENTREZID = names(gene_len_from_exon), length = as.numeric(gene_len_from_exon))

MTG18_exons.counts <- subset(MTG18_exon, rownames(MTG18_exon) %in% gene_len_from_exon$ENTREZID)
dim(MTG18_exons.counts)
## [1] 50267 15928

gene_len_sub = gene_len_from_exon[match(rownames(MTG18_exons.counts), gene_len_from_exon$ENTREZID), ]
gene_len_sub = data.frame(length = gene_len_sub$length, row.names = gene_len_sub$ENTREZID)
head(gene_len_sub)

identical(rownames(MTG18_exons.counts), rownames(gene_len_sub))

counts_to_TPM <- function(counts, geneLength) {
  rpk <- counts / geneLength
  scaling_factor <- colSums(rpk, na.rm = TRUE)
  tpm <- sweep(rpk, 2, scaling_factor, "/") * 1e6
  tpm[, scaling_factor == 0] <- NA 
  return(tpm)
}

MTG18_TPM <- counts_to_TPM(MTG18_exons.counts, gene_len_sub$length)

colSums(MTG18_TPM)
MTG18_TPM[4:8,4:8]

rownames(MTG18_TPM) <- MTG18_gene_info[match(rownames(MTG18_TPM), MTG18_gene_info$entrez_id), 1]

identical(colnames(MTG18_TPM), MTG18_meta.data$sample_name)

## prepare MTG18GluN ref data for celltype mapping
MTG18GluN_meta.data <- subset(MTG18_meta.data, class == "Glutamatergic")
MTG18GluN_TPM <- MTG18_TPM[,MTG18GluN_meta.data$sample_name]
dim(MTG18GluN_TPM)
#[1] 50267 10525

identical(colnames(MTG18GluN_TPM), MTG18GluN_meta.data$sample_name)

MTG18GluN <- CreateSeuratObject(counts = log2(MTG18GluN_TPM + 1), project = "MTG18GluN", min.cells = 1)
MTG18GluN_TPM4mapping <- MTG18GluN@assays$RNA@counts
write.table(MTG18GluN_TPM4mapping, file = "MTG18GluN_TPM4mapping.txt", col.names = TRUE, sep = "\t", quote = FALSE)


## Identify ON and OFF markers
MTG18ref_meta.data <- subset(MTG18_meta.data, !class %in% "no class")
MTG18ref_meta.data <- subset(MTG18ref_meta.data, !(brain_subregion %in% c("L1", "L2", "L3") & class == "Glutamatergic"))

table(MTG18ref_meta.data$class)
## GABAergic     Glutamatergic   Non-neuronal 
## 4164          6861            914

table(MTG18ref_meta.data$cluster)

MTG18ref_meta.data$cell_type <- MTG18ref_meta.data$cluster
MTG18ref_meta.data$cell_type <- sub("\\s.*", "", MTG18ref_meta.data$cell_type)
table(MTG18ref_meta.data$cell_type)
## Astro  Endo   Exc    Inh     Micro  Oligo   OPC 
## 291    9      6861   4164    63     313     238

MTG18.refdata <- MTG18_TPM[,MTG18ref_meta.data$sample_name]

MTG18_ref <- CreateSeuratObject(counts = log2(MTG18.refdata + 1), project = "MTG18_ref", min.cells = 1)
MTG18_ref$cell_type <- MTG18ref_meta.data[match(colnames(MTG18.refdata), MTG18ref_meta.data$sample_name), 35]

MTG18_ref
## An object of class Seurat 
## 47991 features across 11939 samples within 1 assay 
## Active assay: RNA (47991 features, 0 variable features)

Idents(MTG18_ref) <- "cell_type"

MTG18_Allmarkers <- FindAllMarkers(MTG18_ref, only.pos = TRUE, min.pct = 0.9, logfc.threshold = 3.2)

ON_proportion <- rowSums(MTG18.refdata > 10) / ncol(MTG18.refdata)
ON_matrix <- MTG18.refdata[ON_proportion > 0.75, ]

cat("Ori_gene:", nrow(MTG18.refdata), "\n")
cat("filtered_gene:", nrow(ON_matrix), "\n")
cat("filtered_gene_example:", head(rownames(ON_matrix)), "\n")

MTG_deepPC_ONmarkers <- subset(MTG18_Allmarkers, !rownames(MTG18_Allmarkers) %in% rownames(ON_matrix))
MTG_deepPC_ONmarkers <- subset(MTG_deepPC_ONmarkers, cluster == "Exc")
table(duplicated(MTG_deepPC_ONmarkers$gene))

MTG_deepPC_ONmarkers_AverExp <- MTG18.refdata[MTG_deepPC_ONmarkers$gene,]
MTG18deepPC_ID <- subset(MTG18ref_meta.data, cell_type == "Exc")
MTG_deepPC_ONmarkers_AverExp <- MTG_deepPC_ONmarkers_AverExp[,MTG18deepPC_ID$sample_name]
MTG_deepPC_ONmarkers_AverExp <- as.data.frame(rowMeans(MTG_deepPC_ONmarkers_AverExp, na.rm = TRUE))
MTG_deepPC_ONmarkers_AverExp <-  subset(MTG_deepPC_ONmarkers_AverExp, MTG_deepPC_ONmarkers_AverExp > 100)

MTG_deepPC_ONmarkers <- subset(MTG_deepPC_ONmarkers, rownames(MTG_deepPC_ONmarkers) %in% rownames(MTG_deepPC_ONmarkers_AverExp))

Astro_OFF_markers <- FindMarkers(MTG18_ref, ident.1 = "Astro", ident.2 = "Exc", only.pos = TRUE, logfc.threshold = 3.2)
Endo_OFF_markers <- FindMarkers(MTG18_ref, ident.1 = "Endo", ident.2 = "Exc", only.pos = TRUE, logfc.threshold = 3.2)
Micro_OFF_markers <- FindMarkers(MTG18_ref, ident.1 = "Micro", ident.2 = "Exc", only.pos = TRUE, logfc.threshold = 3.2)
Oligo_OFF_markers <- FindMarkers(MTG18_ref, ident.1 = "Oligo", ident.2 = "Exc", only.pos = TRUE, logfc.threshold = 3.2)
OPC_OFF_markers <- FindMarkers(MTG18_ref, ident.1 = "OPC", ident.2 = "Exc", only.pos = TRUE, logfc.threshold = 3.2)
Inh_OFF_markers <- FindMarkers(MTG18_ref, ident.1 = "Inh", ident.2 = "Exc", only.pos = TRUE, logfc.threshold = 3.2)

Astro_OFF_markers <- subset(Astro_OFF_markers, pct.2 < 0.1)
Endo_OFF_markers <- subset(Endo_OFF_markers, pct.2 < 0.1)
Micro_OFF_markers <- subset(Micro_OFF_markers, pct.2 < 0.1)
Oligo_OFF_markers <- subset(Oligo_OFF_markers, pct.2 < 0.1)
OPC_OFF_markers <- subset(OPC_OFF_markers, pct.2 < 0.1)
Inh_OFF_markers <- subset(Inh_OFF_markers, pct.2 < 0.1)

Astro_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "Astro")
Astro_OFF_matrix <- MTG18.refdata[,Astro_OFF_matrix$sample_name]

Astro_OFF_proportion <- rowSums(Astro_OFF_matrix < 10) / ncol(Astro_OFF_matrix)
Astro_OFF_matrix <- Astro_OFF_matrix[Astro_OFF_proportion > 0.5, ]
Astro_OFF_markers <- subset(Astro_OFF_markers, !rownames(Astro_OFF_markers) %in% rownames(Astro_OFF_matrix))

Endo_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "Endo")
Endo_OFF_matrix <- MTG18.refdata[,Endo_OFF_matrix$sample_name]

Endo_OFF_proportion <- rowSums(Endo_OFF_matrix < 10) / ncol(Endo_OFF_matrix)
Endo_OFF_matrix <- Endo_OFF_matrix[Endo_OFF_proportion > 0.5, ]
Endo_OFF_markers <- subset(Endo_OFF_markers, !rownames(Endo_OFF_markers) %in% rownames(Endo_OFF_matrix))

Micro_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "Micro")
Micro_OFF_matrix <- MTG18.refdata[,Micro_OFF_matrix$sample_name]

Micro_OFF_proportion <- rowSums(Micro_OFF_matrix < 10) / ncol(Micro_OFF_matrix)
Micro_OFF_matrix <- Micro_OFF_matrix[Micro_OFF_proportion > 0.5, ]
Micro_OFF_markers <- subset(Micro_OFF_markers, !rownames(Micro_OFF_markers) %in% rownames(Micro_OFF_matrix))

Oligo_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "Oligo")
Oligo_OFF_matrix <- MTG18.refdata[,Oligo_OFF_matrix$sample_name]

Oligo_OFF_proportion <- rowSums(Oligo_OFF_matrix < 10) / ncol(Oligo_OFF_matrix)
Oligo_OFF_matrix <- Oligo_OFF_matrix[Oligo_OFF_proportion > 0.5, ]
Oligo_OFF_markers <- subset(Oligo_OFF_markers, !rownames(Oligo_OFF_markers) %in% rownames(Oligo_OFF_matrix))

OPC_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "OPC")
OPC_OFF_matrix <- MTG18.refdata[,OPC_OFF_matrix$sample_name]

OPC_OFF_proportion <- rowSums(OPC_OFF_matrix < 10) / ncol(OPC_OFF_matrix)
OPC_OFF_matrix <- OPC_OFF_matrix[OPC_OFF_proportion > 0.5, ]
OPC_OFF_markers <- subset(OPC_OFF_markers, !rownames(OPC_OFF_markers) %in% rownames(OPC_OFF_matrix))

Inh_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "Inh")
Inh_OFF_matrix <- MTG18.refdata[,Inh_OFF_matrix$sample_name]

Inh_OFF_proportion <- rowSums(Inh_OFF_matrix < 10) / ncol(Inh_OFF_matrix)
Inh_OFF_matrix <- Inh_OFF_matrix[Inh_OFF_proportion > 0.5, ]
Inh_OFF_markers <- subset(Inh_OFF_markers, !rownames(Inh_OFF_markers) %in% rownames(Inh_OFF_matrix))

Glu_OFF_matrix <- subset(MTG18ref_meta.data, cell_type == "Exc")
Glu_OFF_matrix <- MTG18.refdata[,Glu_OFF_matrix$sample_name]

Glu_OFF_proportion <- rowSums(Glu_OFF_matrix > 10) / ncol(Glu_OFF_matrix)
Glu_OFF_matrix <- Glu_OFF_matrix[Glu_OFF_proportion > 0.33, ]

Astro_OFF_markers <- subset(Astro_OFF_markers, !rownames(Astro_OFF_markers) %in% rownames(Glu_OFF_matrix))
Endo_OFF_markers <- subset(Endo_OFF_markers, !rownames(Endo_OFF_markers) %in% rownames(Glu_OFF_matrix))
Micro_OFF_markers <- subset(Micro_OFF_markers, !rownames(Micro_OFF_markers) %in% rownames(Glu_OFF_matrix))
Oligo_OFF_markers <- subset(Oligo_OFF_markers, !rownames(Oligo_OFF_markers) %in% rownames(Glu_OFF_matrix))
OPC_OFF_markers <- subset(OPC_OFF_markers, !rownames(OPC_OFF_markers) %in% rownames(Glu_OFF_matrix))
Inh_OFF_markers <- subset(Inh_OFF_markers, !rownames(Inh_OFF_markers) %in% rownames(Glu_OFF_matrix))

MTG_deepPC_ONmarkers$cluster <- "Exc on"

Astro_OFF_markers$cluster <- "Ast"
Astro_OFF_markers$gene <- rownames(Astro_OFF_markers)

Endo_OFF_markers$cluster <- "End"
Endo_OFF_markers$gene <- rownames(Endo_OFF_markers)

Micro_OFF_markers$cluster <- "Mic"
Micro_OFF_markers$gene <- rownames(Micro_OFF_markers)

Oligo_OFF_markers$cluster <- "Oli"
Oligo_OFF_markers$gene <- rownames(Oligo_OFF_markers)

OPC_OFF_markers$cluster <- "OPC"
OPC_OFF_markers$gene <- rownames(OPC_OFF_markers)

Inh_OFF_markers$cluster <- "Inh"
Inh_OFF_markers$gene <- rownames(Inh_OFF_markers)

MTG18deepPC_ONOFF_markers <- rbind(MTG_deepPC_ONmarkers, Astro_OFF_markers, Endo_OFF_markers, Micro_OFF_markers, Oligo_OFF_markers, OPC_OFF_markers, Inh_OFF_markers)

table(duplicated(MTG18deepPC_ONOFF_markers$gene))
## FALSE  TRUE 
## 182    22

duplicate_rows <- MTG18deepPC_ONOFF_markers[duplicated(MTG18deepPC_ONOFF_markers$gene), ]
dup_markers <- unique(duplicate_rows$gene)

MTG18_dup_markers_Exp <- as.data.frame(AverageExpression(MTG18_ref, features = dup_markers, assays = "RNA", slot = "counts")) 

Astro_OFF_markers <- subset(Astro_OFF_markers, !gene %in% c("HTRA1", "FGFR2", "PLEKHB1"))
Endo_OFF_markers <- subset(Endo_OFF_markers, !gene %in% c("MT2A", "SLC7A11"))
OPC_OFF_markers <- subset(OPC_OFF_markers, !gene %in% c("S100B", "COL9A2", "PLLP", "PON2", "DOCK10"))
Micro_OFF_markers <- subset(Micro_OFF_markers, !gene %in% c("TLR4", "SLC1A3", "HTRA1", "SPP1"))
Oligo_OFF_markers <- subset(Oligo_OFF_markers, !gene %in% c("HTRA1", "NDRG2", "S100B", "TSC22D4", "BCAS1", "OLIG1"))
Inh_OFF_markers <- subset(Inh_OFF_markers, !gene %in% c("DOCK10", "ZNF536"))

MTG18deepPC_ONOFF_markers_unique <- rbind(MTG_deepPC_ONmarkers, Astro_OFF_markers, Endo_OFF_markers, Micro_OFF_markers, Oligo_OFF_markers, OPC_OFF_markers, Inh_OFF_markers)



# loading data for cell type mapping
patchseq.data.1 <- fread("data/NCBI ID feature_counts 1.csv", header = T, data.table = F)
patchseq.data.2 <- fread("data/NCBI ID feature_counts 2.csv", header = T, data.table = F)

patchseq.data <- merge(patchseq.data.1[,1:116], patchseq.data.2[,1:16], by = "gene_id", all = TRUE)

table(duplicated(patchseq.data$gene_id))
## FALSE 
## 48354  

table(duplicated(patchseq.data.1$gene_name))
## FALSE  TRUE 
## 44091  4263

patchseq.count <- column_to_rownames(patchseq.data, var = "gene_id")

geneid_efflen <- subset(patchseq.data.1, select = c("gene_id", "gene_length"))
colnames(geneid_efflen) <- c("geneid", "efflen")  
efflen <- geneid_efflen[match(rownames(patchseq.count), geneid_efflen$geneid), "efflen"]

counts2TPM <- function(count = patchseq.count, efflength = efflen){
  RPK <- count/(efflength/1000)
  PMSC_rpk <- sum(RPK)/1e6
  RPK/PMSC_rpk              
}

patchseq.TPM <- as.data.frame(apply(patchseq.count, 2, counts2TPM))
colSums(patchseq.TPM) 

NCBI_gene_symbol <- patchseq.data.1[match(rownames(patchseq.count), patchseq.data.1$gene_id), "gene_name"]
table(duplicated(NCBI_gene_symbol))

patchseq.count <- aggregate(patchseq.count, by = list(NCBI_gene_symbol), FUN = sum)
patchseq.count$Group.1 <- gsub("\t", "", patchseq.count$Group.1)
patchseq.count <- column_to_rownames(patchseq.count,'Group.1')

patchseq.TPM <- aggregate(patchseq.TPM, by = list(NCBI_gene_symbol), FUN = sum)
patchseq.TPM$Group.1 <- gsub("\t", "", patchseq.TPM$Group.1)
patchseq.TPM <- column_to_rownames(patchseq.TPM,'Group.1')

head(patchseq.TPM)


## QC for patch-seq mapping data (calculate contamination score)
MTG18deepPC4score <- subset(MTG18_ref, ident = "Exc")
MTG18deepPC4score <- MTG18deepPC4score@assays$RNA@counts

PC_Ast.mtx <- MTG18deepPC4score[rownames(MTG18deepPC4score) %in% rownames(Astro_OFF_markers), ]
PC_Ast.sums <- colSums(PC_Ast.mtx, na.rm = TRUE)
d_PC_Ast <- median(PC_Ast.sums, na.rm = TRUE)

PC_End.mtx <- MTG18deepPC4score[rownames(MTG18deepPC4score) %in% rownames(Endo_OFF_markers), ]
PC_End.sums <- colSums(PC_End.mtx, na.rm = TRUE)
d_PC_End <- median(PC_End.sums, na.rm = TRUE)

PC_Mic.mtx <- MTG18deepPC4score[rownames(MTG18deepPC4score) %in% rownames(Micro_OFF_markers), ]
PC_Mic.sums <- colSums(PC_Mic.mtx, na.rm = TRUE)
d_PC_Mic <- median(PC_Mic.sums, na.rm = TRUE)

PC_Oli.mtx <- MTG18deepPC4score[rownames(MTG18deepPC4score) %in% rownames(Oligo_OFF_markers), ]
PC_Oli.sums <- colSums(PC_Oli.mtx, na.rm = TRUE)
d_PC_Oli <- median(PC_Oli.sums, na.rm = TRUE)

PC_OPC.mtx <- MTG18deepPC4score[rownames(MTG18deepPC4score) %in% rownames(OPC_OFF_markers), ]
PC_OPC.sums <- colSums(PC_OPC.mtx, na.rm = TRUE)
d_PC_OPC <- median(PC_OPC.sums, na.rm = TRUE)

PC_Inh.mtx <- MTG18deepPC4score[rownames(MTG18deepPC4score) %in% rownames(Inh_OFF_markers), ]
PC_Inh.sums <- colSums(PC_Inh.mtx, na.rm = TRUE)
d_PC_Inh <- median(PC_Inh.sums, na.rm = TRUE)

table(MTG18_ref@active.ident)
## Inh   Exc    Oligo   OPC   Astro  Micro  Endo 
## 4164  6861   313     238   291    63     9

Ast_Ast.mtx <- subset(MTG18_ref, ident = "Astro")
Ast_Ast.mtx <- Ast_Ast.mtx@assays$RNA@counts
Ast_Ast.mtx <- Ast_Ast.mtx[rownames(Ast_Ast.mtx) %in% rownames(Astro_OFF_markers), ]
Ast_Ast.sums <- colSums(Ast_Ast.mtx, na.rm = TRUE)
d_Ast_Ast <- median(Ast_Ast.sums, na.rm = TRUE)

End_End.mtx <- subset(MTG18_ref, ident = "Endo")
End_End.mtx <- End_End.mtx@assays$RNA@counts
End_End.mtx <- End_End.mtx[rownames(End_End.mtx) %in% rownames(Endo_OFF_markers), ]
End_End.sums <- colSums(End_End.mtx, na.rm = TRUE)
d_End_End <- median(End_End.sums, na.rm = TRUE)

Mic_Mic.mtx <- subset(MTG18_ref, ident = "Micro")
Mic_Mic.mtx <- Mic_Mic.mtx@assays$RNA@counts
Mic_Mic.mtx <- Mic_Mic.mtx[rownames(Mic_Mic.mtx) %in% rownames(Micro_OFF_markers), ]
Mic_Mic.sums <- colSums(Mic_Mic.mtx, na.rm = TRUE)
d_Mic_Mic <- median(Mic_Mic.sums, na.rm = TRUE)

Oli_Oli.mtx <- subset(MTG18_ref, ident = "Oligo")
Oli_Oli.mtx <- Oli_Oli.mtx@assays$RNA@counts
Oli_Oli.mtx <- Oli_Oli.mtx[rownames(Oli_Oli.mtx) %in% rownames(Oligo_OFF_markers), ]
Oli_Oli.sums <- colSums(Oli_Oli.mtx, na.rm = TRUE)
d_Oli_Oli <- median(Oli_Oli.sums, na.rm = TRUE)

OPC_OPC.mtx <- subset(MTG18_ref, ident = "OPC")
OPC_OPC.mtx <- OPC_OPC.mtx@assays$RNA@counts
OPC_OPC.mtx <- OPC_OPC.mtx[rownames(OPC_OPC.mtx) %in% rownames(OPC_OFF_markers), ]
OPC_OPC.sums <- colSums(OPC_OPC.mtx, na.rm = TRUE)
d_OPC_OPC <- median(OPC_OPC.sums, na.rm = TRUE)

Inh_Inh.mtx <- subset(MTG18_ref, ident = "Inh")
Inh_Inh.mtx <- Inh_Inh.mtx@assays$RNA@counts
Inh_Inh.mtx <- Inh_Inh.mtx[rownames(Inh_Inh.mtx) %in% rownames(Inh_OFF_markers), ]
Inh_Inh.sums <- colSums(Inh_Inh.mtx, na.rm = TRUE)
d_Inh_Inh <- median(Inh_Inh.sums, na.rm = TRUE)


## calculating contamination score for patchseq data
meta.data <- read.csv(file = "data/meta.data.csv", header = TRUE)
patchseq.data4QC <- patchseq.TPM
colnames(patchseq.data4QC) <- meta.data[match(colnames(patchseq.data4QC), meta.data$Sample.ID), 5]

patchseq4QC <- CreateSeuratObject(counts = log2(patchseq.data4QC + 1), project = "VENQC_project", min.cells = 1)

patchseq4QC@meta.data$sample_ID <- meta.data[match(colnames(patchseq.data4QC), meta.data$Cell.name), 1]
patchseq4QC@meta.data$cell_name <- meta.data[match(colnames(patchseq.data4QC), meta.data$Cell.name), 5]
patchseq4QC@meta.data$M_type <- meta.data[match(colnames(patchseq.data4QC), meta.data$Cell.name), 6]
patchseq4QC@meta.data$brain_region <- meta.data[match(colnames(patchseq.data4QC), meta.data$Cell.name), 10]

patchseq_TPM4QC <- patchseq4QC@assays$RNA@counts

patch_Ast.mtx <- patchseq_TPM4QC[rownames(patchseq_TPM4QC) %in% rownames(Astro_OFF_markers), ]
patch_Ast.sums <- colSums(patch_Ast.mtx, na.rm = TRUE)
CS.patch_Ast <- as.matrix((patch_Ast.sums - d_PC_Ast)/(d_Ast_Ast - d_PC_Ast))
patchseq4QC@meta.data$Ast_contam <- CS.patch_Ast[match(colnames(patchseq.data4QC), rownames(CS.patch_Ast)), 1]

patch_End.mtx <- patchseq_TPM4QC[rownames(patchseq_TPM4QC) %in% rownames(Endo_OFF_markers), ]
patch_End.sums <- colSums(patch_End.mtx, na.rm = TRUE)
CS.patch_End <- as.matrix((patch_End.sums - d_PC_End)/(d_End_End - d_PC_End))
patchseq4QC@meta.data$End_contam <- CS.patch_End[match(colnames(patchseq.data4QC), rownames(CS.patch_End)), 1]

patch_Mic.mtx <- patchseq_TPM4QC[rownames(patchseq_TPM4QC) %in% rownames(Micro_OFF_markers), ]
patch_Mic.sums <- colSums(patch_Mic.mtx, na.rm = TRUE)
CS.patch_Mic <- as.matrix((patch_Mic.sums - d_PC_Mic)/(d_Mic_Mic - d_PC_Mic))
patchseq4QC@meta.data$Mic_contam <- CS.patch_Mic[match(colnames(patchseq.data4QC), rownames(CS.patch_Mic)), 1]

patch_Oli.mtx <- patchseq_TPM4QC[rownames(patchseq_TPM4QC) %in% rownames(Oligo_OFF_markers), ]
patch_Oli.sums <- colSums(patch_Oli.mtx, na.rm = TRUE)
CS.patch_Oli <- as.matrix((patch_Oli.sums - d_PC_Oli)/(d_Oli_Oli - d_PC_Oli))
patchseq4QC@meta.data$Oli_contam <- CS.patch_Oli[match(colnames(patchseq.data4QC), rownames(CS.patch_Oli)), 1]

patch_OPC.mtx <- patchseq_TPM4QC[rownames(patchseq_TPM4QC) %in% rownames(OPC_OFF_markers), ]
patch_OPC.sums <- colSums(patch_OPC.mtx, na.rm = TRUE)
CS.patch_OPC <- as.matrix((patch_OPC.sums - d_PC_OPC)/(d_OPC_OPC - d_PC_OPC))
patchseq4QC@meta.data$OPC_contam <- CS.patch_OPC[match(colnames(patchseq.data4QC), rownames(CS.patch_OPC)), 1]

patch_Inh.mtx <- patchseq_TPM4QC[rownames(patchseq_TPM4QC) %in% rownames(Inh_OFF_markers), ]
patch_Inh.sums <- colSums(patch_Inh.mtx, na.rm = TRUE)
CS.patch_Inh <- as.matrix((patch_Inh.sums - d_PC_Inh)/(d_Inh_Inh - d_PC_Inh))
patchseq4QC@meta.data$Inh_contam <- CS.patch_Inh[match(colnames(patchseq.data4QC), rownames(CS.patch_Inh)), 1]

patchseq4QC@meta.data <- patchseq4QC@meta.data %>%
  mutate(contamination_score = rowSums(across(all_of(c("Ast_contam", "End_contam", "Mic_contam", "Oli_contam", "OPC_contam", "Inh_contam"))), na.rm = TRUE))

MTG18deepPC <- subset(MTG18_ref, ident = "Exc")
MTG18deepPC_ONOFF_Exp <- AverageExpression(MTG18deepPC, features = MTG18deepPC_ONOFF_markers_unique$gene, assays = "RNA", slot = "counts")

typeof(MTG18deepPC_ONOFF_Exp)
head(MTG18deepPC_ONOFF_Exp$RNA)
MTG18deepPC_ONOFF_Exp <- as.data.frame(MTG18deepPC_ONOFF_Exp)
colnames(MTG18deepPC_ONOFF_Exp) <- "MTG18deepPC"
MTG18deepPC_ONOFF_Exp$gene <- rownames(MTG18deepPC_ONOFF_Exp)

patchseq_ONOFF_Exp <- AverageExpression(patchseq4QC, features = MTG18deepPC_ONOFF_markers_unique$gene, assays = "RNA", slot = "counts")

patchseq_ONOFF_Exp <- as.data.frame(patchseq_ONOFF_Exp)
colnames(patchseq_ONOFF_Exp) <- "patchseq"
patchseq_ONOFF_Exp$gene <- rownames(patchseq_ONOFF_Exp)

ONOFFgene_cor_mtx <- merge(MTG18deepPC_ONOFF_Exp, patchseq_ONOFF_Exp, by = "gene", all = TRUE)
ONOFFgene_cor_mtx <- column_to_rownames(ONOFFgene_cor_mtx, var = "gene")

ONOFFgene_cor_result <- shapiro.test(ONOFFgene_cor_mtx$patchseq)
print(ONOFFgene_cor_result)

ONOFFgene_cor <- cor.test(ONOFFgene_cor_mtx$MTG18deepPC, ONOFFgene_cor_mtx$patchseq, method = "spearman")
cor_value <- round(ONOFFgene_cor$estimate, 2)
cor_p_value <- format.pval(ONOFFgene_cor$p.value, digits = 2)

patchseq4QC@meta.data <- patchseq4QC@meta.data %>% mutate(scaled_score = contamination_score * 0.3105976)

ONOFFgene_cor_mtx$cluster <- MTG18deepPC_ONOFF_markers_unique[match(rownames(ONOFFgene_cor_mtx), MTG18deepPC_ONOFF_markers_unique$gene), 6]

ggplot(ONOFFgene_cor_mtx, aes(x = MTG18deepPC, y = patchseq, color = cluster)) + geom_point(size = 3) +
  annotate("text", x = min(ONOFFgene_cor_mtx$MTG18deepPC), y = max(ONOFFgene_cor_mtx$patchseq),
           label = paste0("spearman r = ", cor_value, "\nP = ", cor_p_value), hjust = 0, vjust = 1, size = 4, color = "black") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black") +
  theme(axis.ticks.length = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.5), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(breaks = seq(0,10, by = 2), limits = c(0,10)) + scale_x_continuous(breaks = seq(0,10, by = 2), limits = c(0,10))

ggplot(ONOFFgene_cor_mtx, aes(x = MTG18deepPC, y = patchseq, color = cluster)) + geom_point(size = 3) + facet_wrap(~ cluster, scales = "free") +
  theme(axis.ticks.length = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.5), panel.background = element_rect(fill = "white"))


## cell type mapping data preparation (filter patchseq data)
OFF_markers4QC <- subset(MTG18deepPC_ONOFF_markers_unique, cluster != "Exc on")
OFF_markers4QC_ID <- subset(patchseq.data.1, gene_name %in% OFF_markers4QC$gene)
table(duplicated(OFF_markers4QC_ID$gene_name))
## FALSE  TRUE 
## 176    47

table(duplicated(OFF_markers4QC_ID$gene_id))

patchseq.count <- column_to_rownames(patchseq.data, var = "gene_id")
patchseq.count_QC <- subset(patchseq.count, !rownames(patchseq.count) %in% OFF_markers4QC_ID$gene_id)

efflen <- geneid_efflen[match(rownames(patchseq.count_QC), geneid_efflen$geneid), "efflen"]

counts2TPM <- function(count = patchseq.count_QC, efflength = efflen){
  RPK <- count/(efflength/1000)
  PMSC_rpk <- sum(RPK)/1e6
  RPK/PMSC_rpk              
}

patchseq.TPM_QC <- as.data.frame(apply(patchseq.count_QC, 2, counts2TPM))
colSums(patchseq.TPM_QC) 

NCBI_gene_symbol <- patchseq.data.1[match(rownames(patchseq.count_QC), patchseq.data.1$gene_id), "gene_name"]
table(duplicated(NCBI_gene_symbol))

patchseq.count_QC <- aggregate(patchseq.count_QC, by = list(NCBI_gene_symbol), FUN = sum)
patchseq.count_QC$Group.1 <- gsub("\t", "", patchseq.count_QC$Group.1)
patchseq.count_QC <- column_to_rownames(patchseq.count_QC,'Group.1')

patchseq.TPM_QC <- aggregate(patchseq.TPM_QC, by = list(NCBI_gene_symbol), FUN = sum)
patchseq.TPM_QC$Group.1 <- gsub("\t", "", patchseq.TPM_QC$Group.1)
patchseq.TPM_QC <- column_to_rownames(patchseq.TPM_QC,'Group.1')

head(patchseq.TPM_QC)

colnames(patchseq.TPM_QC) <- meta.data[match(colnames(patchseq.TPM_QC), meta.data$Sample.ID), 5]
patchseq.TPM_QC <- subset(patchseq.TPM_QC, select = -c(PC43, PC64, TRI33, TRI32, TRI20, PC52))

patchseq4mapping <- CreateSeuratObject(counts = log2(patchseq.TPM_QC + 1), project = "patchseq2mapping", min.cells = 1)

patchseq4mapping@meta.data$sample_ID <- meta.data[match(colnames(patchseq.TPM_QC), meta.data$Cell.name), 1]
patchseq4mapping@meta.data$cell_name <- meta.data[match(colnames(patchseq.TPM_QC), meta.data$Cell.name), 5]
patchseq4mapping@meta.data$M_type <- meta.data[match(colnames(patchseq.TPM_QC), meta.data$Cell.name), 6]
patchseq4mapping@meta.data$brain_region <- meta.data[match(colnames(patchseq.TPM_QC), meta.data$Cell.name), 10]

patchseq_TPM4mapping <- patchseq4mapping@assays$RNA@counts
write.table(patchseq_TPM4mapping, file = "VEN4mapping.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(patchseq4mapping@meta.data, file = "VEN4mapping_meta_data.txt", col.names = TRUE, sep = "\t", quote = FALSE)

patchseq4mapping <- FindVariableFeatures(patchseq4mapping, selection.method = "vst", nfeatures = 3000)
VariableFeaturePlot(patchseq4mapping, log = TRUE)
patchseq4mapping <- ScaleData(patchseq4mapping, features = rownames(patchseq4mapping), verbose = FALSE)

VEN4mapping_mat <- as.data.frame(patchseq4mapping@assays$RNA@scale.data)

ETIT_markers <- c("POU3F1", "FAM84B", "BCL11B", "NRP1", "DSCAML1", "NFIB", "NTNG1", "SYT2", "SEMA6A", "DAB1", "SEZ6",
                  "MYO16", "ASAP1", "SEMA3D", "GFRA1", "PTPRM", "CRTAC1", "IGSF21", "SERPINE2", "ABAT", "GRIK2", "SLC1A1",
                  "ETV5", "SEZ6L", "AKAP12", "PTPRF", "SORCS2", "ADRA1A", "SPARC", "SLC6A1", "CDHR3", "GRIK1", "CACNA1H",
                  "CACNA2D2", "SCN4B", "TRPC4", "SCN9A", "GRM4", "SLC5A8", "PCSK6", "LOC101929728", "PDE9A", "LOC105378657",
                  "ATP6V1C2", "RNF152", "MEIS2", "DGKD", "ALCAM", "LOC105378653", "FEZF2", "EYA4", "ADCY8", "RXFP1", "NPTX1",
                  "PDZRN3", "DACH1", "LY86-AS1", "SPATS2L", "FIGN", "LINC01202", "MGAT5B", "THEMIS")

ETIT_markers <- c("POU3F1", "FAM84B", "BCL11B", "NRP1", "DSCAML1", "NFIB", "NTNG1", "SYT2", "SEMA6A", "DAB1", "SEZ6",
                  "MYO16", "ASAP1", "SEMA3D", "GFRA1", "PTPRM", "CRTAC1", "IGSF21", "SERPINE2", "ABAT", "GRIK2", "SLC1A1",
                  "ETV5", "SEZ6L", "AKAP12", "PTPRF", "SORCS2", "ADRA1A", "SPARC", "SLC6A1", "CDHR3", "CACNA1H",
                  "CACNA2D2", "SCN4B", "TRPC4", "SCN9A", "GRM4", "SLC5A8", "PCSK6", "LOC101929728", "PDE9A", "LOC105378657",
                  "ATP6V1C2", "RNF152", "MEIS2", "DGKD", "ALCAM", "LOC105378653", "FEZF2", "EYA4", "ADCY8", "RXFP1", "NPTX1",
                  "PDZRN3", "DACH1", "LY86-AS1", "SPATS2L", "FIGN", "LINC01202", "MGAT5B", "THEMIS")

ETIT_markers_mat <- VEN4mapping_mat[ETIT_markers,]
pheatmap(ETIT_markers_mat, fontsize_row = 5, clustering_method = "ward.D", fontsize_col = 5, cluster_rows = F,
         color = c("blue", "lightgrey", "firebrick3"))



# loading and QC VEN patchseq data for clustering and other downstream analysis
## remove OFF markers
VEN_data_1 <- fread("data/feature_counts 1.csv", header = T, data.table = F)
VEN_data_2 <- fread("data/feature_counts 2.csv", header = T, data.table = F)
VEN_data_3 <- fread("data/feature_counts 3.csv", header = T, data.table = F)
VEN_data_4 <- fread("data/feature_counts 4.csv", header = T, data.table = F)
VEN_data_5 <- fread("data/feature_counts 5.csv", header = T, data.table = F)

VEN_counts1 <- VEN_data_1[,1:84]
VEN_counts2 <- VEN_data_2[,1:24]
VEN_counts3 <- VEN_data_3[,1:7]
VEN_counts4 <- VEN_data_4[,1:19]
VEN_counts5 <- VEN_data_5[,1:16]

VEN_counts <- merge(VEN_counts1, VEN_counts2, by = "gene_id", all = TRUE)
VEN_counts <- merge(VEN_counts, VEN_counts3, by = "gene_id", all = TRUE)
VEN_counts <- merge(VEN_counts, VEN_counts4, by = "gene_id", all = TRUE)
VEN_counts <- merge(VEN_counts, VEN_counts5, by = "gene_id", all = TRUE)
VEN_counts <- column_to_rownames(VEN_counts, var = "gene_id")

OFF_markers_miss <- subset(OFF_markers4QC, !rownames(OFF_markers4QC) %in% VEN_data_1$gene_name)
OFF_markers_miss$gene
## "LOC101930275", "PPAP2B", "LOC105376917", "LOC105369345", "FYB", "LOC101929249", "LOC105379054", "LPPR1"

OFF_markers_miss_NCBI_ID <- patchseq.data.1[match(OFF_markers_miss$gene, patchseq.data.1$gene_name), "gene_id"]
OFF_markers_miss_NCBI_ID
## "101930275", "8613", "105376917", "105369345", "2533", "101929249", "105379054", "54886_1" 
## "ENSG00000253944", "ENSG00000162407", "ENSG00000082074", "ENSG00000249835", "ENSG00000148123"

OFF_markers_overlap <- subset(VEN_data_1, gene_name %in% OFF_markers4QC$gene)
table(duplicated(OFF_markers_overlap$gene_name))
## FALSE 
## 168

table(duplicated(OFF_markers_overlap$gene_id))
table(duplicated(OFF_markers4QC$gene))
## FALSE 
## 176

VEN_counts_QC <- subset(VEN_counts, !rownames(VEN_counts) %in% c("ENSG00000253944", "ENSG00000162407", "ENSG00000082074", "ENSG00000249835", "ENSG00000148123"))
VEN_counts_QC <- subset(VEN_counts_QC, !rownames(VEN_counts_QC) %in% OFF_markers_overlap$gene_id)

geneid_efflen <- subset(VEN_data_1, select = c("gene_id", "gene_length"))
colnames(geneid_efflen) <- c("geneid", "efflen")  
efflen <- geneid_efflen[match(rownames(VEN_counts_QC), geneid_efflen$geneid), "efflen"]

counts2TPM <- function(count = VEN_counts_QC, efflength = efflen){
  RPK <- count/(efflength/1000)
  PMSC_rpk <- sum(RPK)/1e6
  RPK/PMSC_rpk              
}

VEN.TPM_QC <- as.data.frame(apply(VEN_counts_QC, 2, counts2TPM))
colSums(VEN.TPM_QC) 

gene_symbol <- VEN_data_1[match(rownames(VEN_counts_QC), VEN_data_1$gene_id), "gene_name"]
table(duplicated(gene_symbol))
## FALSE  TRUE 
## 56996  1566 

VEN_counts_QC <- aggregate(VEN_counts_QC, by = list(gene_symbol), FUN = sum)
VEN_counts_QC$Group.1 <- gsub("\t", "", VEN_counts_QC$Group.1)
VEN_counts_QC <- column_to_rownames(VEN_counts_QC,'Group.1')

VEN.TPM_QC <- aggregate(VEN.TPM_QC, by = list(gene_symbol), FUN = sum)
VEN.TPM_QC$Group.1 <- gsub("\t", "", VEN.TPM_QC$Group.1)
VEN.TPM_QC <- column_to_rownames(VEN.TPM_QC,'Group.1')

head(VEN.TPM_QC)

colnames(VEN.TPM_QC) <- meta.data[match(colnames(VEN.TPM_QC), meta.data$Sample.ID), 5]


## Sample filter
VEN.TPM4clust <- VEN.TPM_QC[,colnames(VEN.TPM_QC) %in% colnames(patchseq4mapping)]

all_genes <- subset(VEN_data_1, select = c("gene_id", "gene_name", "gene_chr", "gene_biotype"))
table(duplicated(all_genes$gene_name))

all_genes <- subset(all_genes, all_genes$gene_name %in% rownames(VEN.TPM4clust))
table(duplicated(all_genes$gene_name))
all_genes <- all_genes[!duplicated(all_genes[c("gene_name")]), ]
all_genes$gene_name <- gsub("\t", "", all_genes$gene_name)

allgenes_biotype <- as.data.frame(table(all_genes$gene_biotype))
colnames(allgenes_biotype) <- c("biotype", "count")  

ggplot() + geom_bar(data = allgenes_biotype, aes(x = biotype, y = count), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(), axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 10), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,20000, by = 4000), limits = c(0,20000))

all_genes4clustering <- subset(all_genes, !gene_chr == "X")
all_genes4clustering <- subset(all_genes4clustering, !gene_chr == "Y")
all_genes4clustering <- subset(all_genes4clustering, !gene_chr == "MT")

Human_Mito <- fread("data/Human.MitoCarta3.0.csv", header = T, data.table = F)
table(duplicated(Human_Mito$Symbol))

all_genes4clustering <- subset(all_genes4clustering, !gene_name %in% Human_Mito$Symbol)

VEN.data4clust <- VEN.TPM4clust[all_genes4clustering$gene_name,]
write.table(VEN.data4clust, file = "VEN.data4clustTPM.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VEN4FS <- CreateSeuratObject(counts = log2(VEN.TPM4clust + 1), project = "VEN4FS", min.cells = 1)

VEN4FS@meta.data$sample_ID <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 1]
VEN4FS@meta.data$batch <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 2]
VEN4FS@meta.data$donor <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 3]
VEN4FS@meta.data$cell_name <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 5]
VEN4FS@meta.data$M_type <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 6]
VEN4FS@meta.data$PM_type <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 8]
VEN4FS@meta.data$E_type <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 9]
VEN4FS@meta.data$brain_region <- meta.data[match(colnames(VEN.TPM4clust), meta.data$Cell.name), 10]
VEN4FS@meta.data$contamination_score <- patchseq4QC@meta.data[match(colnames(VEN.TPM4clust), colnames(patchseq.data4QC)), "scaled_score"]

VEN4FS <- ScaleData(VEN4FS, features = rownames(VEN4FS), verbose = FALSE)
write.table(VEN4FS@assays$RNA@scale.data, file = "VEN4FS_scaled.data.txt", col.names = TRUE, sep = "\t", quote = FALSE)



# downstream analysis for clustering
VEN4clust <- CreateSeuratObject(counts = log2(VEN.data4clust + 1), project = "VEN_project", min.cells = 2)
write.table(VEN4clust@assays$RNA@data, file = "VEN_datalog2TPM.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VEN4clust@meta.data$sample_ID <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 1]
VEN4clust@meta.data$batch <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 2]
VEN4clust@meta.data$donor <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 3]
VEN4clust@meta.data$cell_name <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 5]
VEN4clust@meta.data$M_type <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 6]
VEN4clust@meta.data$PM_type <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 8]
VEN4clust@meta.data$E_type <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 9]
VEN4clust@meta.data$brain_region <- meta.data[match(colnames(VEN.data4clust), meta.data$Cell.name), 10]
VEN4clust@meta.data$contamination_score <- patchseq4QC@meta.data[match(colnames(VEN.data4clust), colnames(patchseq.data4QC)), "scaled_score"]

VEN4clust[["percent.mt"]] <- PercentageFeatureSet(VEN4clust, pattern = "^MT-")
write.table(VEN4clust@meta.data, file = "VEN4clust_meta.data.txt", col.names = TRUE, sep = "\t", quote = FALSE)
Idents(VEN4clust) <- "M_type"

QCP1 <- VlnPlot(VEN4clust, features = "nFeature_RNA", pt.size = 1) + scale_y_continuous(breaks = seq(0,50000, by = 10000), limits = c(0,50000)) + NoLegend()
QCP2 <- VlnPlot(VEN4clust, features = "nCount_RNA", pt.size = 1) + scale_y_continuous(breaks = seq(0,200000, by = 40000), limits = c(0,200000)) + NoLegend()
wrap_plots(plots = list(QCP1, QCP2), ncol = 2)

VEN4clust
## An object of class Seurat 
## 50890 features across 124 samples within 1 assay 
## Active assay: RNA (50890 features, 0 variable features)

table(VEN4clust$PM_type)
##  ET_PC  ET_VEN  IT_PC   IT_VEN 
##  31     13      71      9 

table(VEN4clust$E_type)
## 1  2   3 
## 38 12  6   

table(VEN4clust$brain_region)
## Frontal  Parietal   Temporal 
## 41       18         65 

VEN4clust <- FindVariableFeatures(VEN4clust, selection.method = "vst", nfeatures = 2500)
VariableFeaturePlot(VEN4clust, log = TRUE)

HVG_2500 <- subset(all_genes, all_genes$gene_name %in% VEN4clust@assays$RNA@var.features)
HVG_2500_biotype <- as.data.frame(table(HVG_2500$gene_biotype))
colnames(HVG_2500_biotype) <- c("biotype", "count")

ggplot() + geom_bar(data = HVG_2500_biotype, aes(x = biotype, y = count), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(), axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 10), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,500, by = 100), limits = c(0,500))

VEN4clust <- ScaleData(VEN4clust, features = rownames(VEN4clust), verbose = FALSE)

VEN4clust <- RunPCA(VEN4clust, npcs = 50, verbose = FALSE)
VEN4clust
## An object of class Seurat 
## 50890 features across 124 samples within 1 assay 
## Active assay: RNA (50890 features, 2500 variable features)
## 1 dimensional reduction calculated: pca

VizDimLoadings(VEN4clust, dims = 1:9, reduction = "pca")
DimPlot(object = VEN4clust, reduction = "pca", pt.size = 2)
FeatureScatter(VEN4clust, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 2)

VEN4clust <- JackStraw(VEN4clust, num.replicate = 100)
VEN4clust <- ScoreJackStraw(VEN4clust, dims = 1:20)
JackStrawPlot(VEN4clust, dims = 1:20)
ElbowPlot(VEN4clust)

VEN4clust$PM_type <- factor(x = VEN4clust$PM_type, levels = c("ET_VEN", "ET_PC", "IT_VEN", "IT_PC"))

VEN4clust <- FindNeighbors(VEN4clust, reduction = "pca", dims = c(1,2,3,4,5,6,9))
VEN4clust <- FindClusters(VEN4clust, resolution = 0.7)
VEN4clust <- RunUMAP(VEN4clust, reduction = "pca", dims = c(1,2,3,4,5,6,9))

DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "seurat_clusters")
DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "M_type")
DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "PM_type")
DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "E_type")
DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "brain_region")
DimPlot(VEN4clust, reduction = "umap", label = TRUE, pt.size = 3, group.by = "sample_ID") + NoLegend()

VlnPlot(VEN4clust, features = c("SLC17A6", "SLC17A7", "POU3F1", "COL22A1", "PVALB", "SLC32A1", "GAD1", "GAD2", "TH",
                                "SLC6A3", "SLC18A2", "AQP4", "CX3CR1", "MBP"), group.by = "M_type", same.y.lims = TRUE, add.noise = FALSE, slot = "counts")

write.table(VEN4clust@assays$RNA@scale.data, file = "VEN4clust_scaled.data.txt", col.names = TRUE, sep = "\t", quote = FALSE)


## DEG among clusters
All_cluster_markers <- FindAllMarkers(VEN4clust, logfc.threshold = 1, only.pos = TRUE, min.pct = 1)
All_cluster_markers <- subset(All_cluster_markers, p_val_adj < 0.05)
table(All_cluster_markers$cluster)

cluster1_markers <- subset(All_cluster_markers, cluster == 0)
cluster2_markers <- subset(All_cluster_markers, cluster == 1)
cluster1_markers_biotype <- all_genes[match(cluster1_markers$gene, all_genes$gene_name), 4]
cluster2_markers_biotype <- all_genes[match(cluster2_markers$gene, all_genes$gene_name), 4]

table(cluster1_markers_biotype)
## antisense  lincRNA    protein_coding    sense_intronic    sense_overlapping 
## 3          4          176               1                 1 

table(cluster2_markers_biotype)
## antisense              lincRNA              processed_pseudogene 
## 3                      9                    1 
## processed_transcript   protein_coding       sense_intronic 
## 2                      261                  3 
## sense_overlapping      snoRNA               transcribed_unprocessed_pseudogene 
## 1                      1                    1

write.table(All_cluster_markers, file = "All_cluster_markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

All_cluster_markers_cosg <- cosg(VEN4clust, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                                 expressed_pct = 1, n_genes_user = 300)

clust_marker4plot <- All_cluster_markers_cosg$names[1:5,]
clust_marker4plot <- c(clust_marker4plot[,1], clust_marker4plot[,2])

DoHeatmap(VEN4clust, features = All_cluster_markers$gene, draw.lines = FALSE, hjust = 0) + scale_fill_gradientn(colors = c("#2f58a7", "lightgrey", "#FF6347"))

VlnPlot(VEN4clust, features = clust_marker4plot, same.y.lims = TRUE, stack = TRUE, flip = TRUE)
RidgePlot(VEN4clust, features = clust_marker4plot, fill.by = "ident", ncol = 5)

RidgePlot(VEN4clust, features = c("SULF2", "FAM84B", "ADRA1A", "VAT1L", "POU3F1", "FEZF2", "SCN4B", "BCL11B"), stack = TRUE, fill.by = "ident")
RidgePlot(VEN4clust, features = c("NPTX1", "RXFP1", "PDZRN3", "LY86-AS1", "SPATS2L", "MGAT5B","SCN3B", "THEMIS"), stack = TRUE, fill.by = "ident")

All_cluster_markers.df <- bitr(All_cluster_markers$gene, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
## 5.57% of input gene IDs are fail to map...

All_cluster_markers_ID <- unique(All_cluster_markers.df$ENTREZID)

All_cluster_markers_ego <- enrichGO(gene = All_cluster_markers_ID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
d <- godata('org.Hs.eg.db', ont = "BP")

goCls <- pairwise_termsim(All_cluster_markers_ego, method = "Resnik", semData = d)
treeplot(goCls)

All_cluster_markers_ego <- as.data.frame(All_cluster_markers_ego)
All_cluster_markers_ego_top5 <- All_cluster_markers_ego %>% group_by(ONTOLOGY) %>% arrange(p.adjust) %>% slice_head(n = 5)

All_cluster_markers_ego_top5$ONTOLOGY <- factor(x = All_cluster_markers_ego_top5$ONTOLOGY, levels = c("MF", "BP", "CC"))
All_cluster_markers_ego_top5$Description <- factor(All_cluster_markers_ego_top5$Description, levels = rev(All_cluster_markers_ego_top5$Description))

mycol3 <- c('#FF7F50', '#6B8E23', '#6BA5CE')

p <- ggplot(data = All_cluster_markers_ego_top5, aes(x = Count, y = Description, fill = ONTOLOGY)) +
  geom_bar(width = 0.5, stat = "identity") + theme_classic() + 
  scale_x_continuous(expand = c(0,0.5), breaks = seq(0,70, by = 14), limits = c(0,70)) +
  scale_fill_manual(values = alpha(mycol3, 0.8))

p <- p + theme(axis.text.y = element_blank()) +
  geom_text(data = All_cluster_markers_ego_top5, aes(x = 0.1, y = Description, label = Description), size = 5, hjust = 0)

p <- p + geom_text(data = All_cluster_markers_ego_top5, aes(x = 0.1, y = Description, label = geneID, color = -log10(pvalue)), size = 3,
                   fontface = 'italic', hjust = 0, vjust = 2.7) + scale_colour_viridis(option = "G", direction = -1, end = 0.8)

p <- p + labs(title = 'Enriched top 5 GO terms') +
  theme(plot.title = element_text(size = 14, face = 'bold'), axis.title = element_text(size = 14),
        axis.line = element_line(colour = "black"), axis.text = element_text(size = 12, colour = "black"), 
        axis.ticks.y = element_blank())

p

All_cluster_markers_kegg <- enrichKEGG(gene = All_cluster_markers_ID, keyType = "kegg", organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                       qvalueCutoff = 0.2)
clusterProfiler::dotplot(All_cluster_markers_kegg, title = "All_cluster_markers_KEGG")

ETIT_markers <- c("POU3F1", "FAM84B", "BCL11B", "NRP1", "DSCAML1", "NFIB", "NTNG1", "SYT2", "SEMA6A", "DAB1", "SEZ6",
                  "MYO16", "ASAP1", "SEMA3D", "GFRA1", "PTPRM", "CRTAC1", "IGSF21", "SERPINE2", "ABAT", "GRIK2", "SLC1A1",
                  "ETV5", "SEZ6L", "AKAP12", "PTPRF", "SORCS2", "ADRA1A", "SPARC", "SLC6A1", "CDHR3", "GRIK1", "CACNA1H",
                  "CACNA2D2", "SCN4B", "TRPC4", "SCN9A", "GRM4", "SLC5A8", "PCSK6", "LOC101929728", "PDE9A", "LOC105378657",
                  "ATP6V1C2", "RNF152", "MEIS2", "DGKD", "ALCAM", "AL139158.2", "FEZF2", "EYA4", "ADCY8", "RXFP1", "NPTX1",
                  "PDZRN3", "DACH1", "LY86-AS1", "SPATS2L", "FIGN", "LINC01202", "MGAT5B", "THEMIS")

DoHeatmap(VEN4clust, features = ETIT_markers, draw.lines = FALSE, hjust = 0) + scale_fill_gradientn(colors = c("#2f58a7", "lightgrey", "#FF6347"))


## channels and receptors
Idents(VEN4clust) <- "E_type"
VEN4Etype <- subset(VEN4clust, ident = c(1, 2, 3))

VEN4Etype
## An object of class Seurat 
## 50890 features across 56 samples within 1 assay 
## Active assay: RNA (50890 features, 2500 variable features)
## 2 dimensional reductions calculated: pca, umap

DimPlot(VEN4Etype, reduction = "umap", pt.size = 3, group.by = "E_type", shape.by = "seurat_clusters")

VlnPlot(VEN4Etype, features = c("CACNA1A", "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E", "CACNA1F", "CACNA1G", "CACNA1H", "CACNA1I", "CACNA1S",
                                "CACNA2D1", "CACNA2D2", "CACNA2D3", "CACNA2D4", "CACNB1", "CACNB2", "CACNB3", "CACNB4", "CACNG1", "CACNG2",
                                "CACNG3", "CACNG4", "CACNG5", "CACNG6", "CACNG7", "CACNG8"), fill.by = "ident", adjust = 1, stack = TRUE,
        same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("CLCC1", "CLCF1", "CLCN1", "CLCN2", "CLCN3", "CLCN4", "CLCN5", "CLCN6", "CLCN7", "CLCNKA", "CLCNKB"),
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("SCN1A", "SCN1B", "SCN2A", "SCN2B", "SCN3A", "SCN3B", "SCN4A", "SCN4B", "SCN5A", "SCN7A", "SCN8A",
                                "SCN9A", "SCN10A", "SCN11A", "SCNM1", "SCNN1A", "SCNN1B", "SCNN1G"), fill.by = "ident", stack = TRUE,
        same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("HCN1", "HCN2", "HCN3", "HCN4", "KCNA1", "KCNA2", "KCNA3", "KCNA4", "KCNA5", "KCNA7", "KCNA10", "KCNAB1",
                                "KCNAB2", "KCNAB3", "KCNB1", "KCNB2", "KCNC1", "KCNC2", "KCNC3", "KCNC4", "KCND1", "KCND2", "KCND3"),
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("KCNE1", "KCNE2", "KCNE3", "KCNE4", "KCNF1", "KCNG1", "KCNG3", "KCNG4", "KCNH1", "KCNH2", "KCNH3", "KCNH5",
                                "KCNH6", "KCNH7", "KCNH8", "KCNIP1", "KCNIP2", "KCNIP3", "KCNIP4", "KCNJ1", "KCNJ2", "KCNJ3", "KCNJ4"),
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("KCNJ5", "KCNJ6", "KCNJ8", "KCNJ9", "KCNJ10", "KCNJ11", "KCNJ12", "KCNJ13", "KCNJ14", "KCNJ15", "KCNJ16",
                                "KCNK1", "KCNK2", "KCNK3", "KCNK4", "KCNK5", "KCNK6", "KCNK7", "KCNK9", "KCNK10", "KCNK12", "KCNK13",
                                "KCNK15"), fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("KCNK18", "KCNMA1", "KCNMB1", "KCNMB2", "KCNMB4", "KCNN1", "KCNN2", "KCNN3", "KCNN4", "KCNQ1", "KCNQ2",
                                "KCNQ3", "KCNQ5", "KCNRG", "KCNS1", "KCNS2", "KCNS3", "KCNT1", "KCNT2", "KCNU1", "KCNV1", "KCNV2"),
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("GABBR1", "GABBR2", "GABRA1", "GABRA2", "GABRA3", "GABRA4", "GABRA5", "GABRA6", "GABRB1", "GABRB2", "GABRB3",
                                "GABRD", "GABRE", "GABRG1", "GABRG2", "GABRG3", "GABRP",  "GABRQ",  "GABRR1", "GABRR2", "GABRR3"),
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8", "GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRID1",
                                "GRID2", "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5", "GRIN1", "GRIN2A", "GRIN2B", "GRIN2D", "GRIN3A",
                                "GRIN3B", "GRINA"), fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c("DRD1", "DRD2", "DRD3", "DRD4", "DRD5", "HTR1A", "HTR1B", "HTR1D", "HTR1F", "HTR2A", "HTR2B", "HTR2C",
                                "HTR3A", "HTR3B", "HTR4", "HTR5A", "HTR6", "HTR7"), 
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()

VlnPlot(VEN4Etype, features = c('CHRM1', 'CHRM2', 'CHRM3', 'CHRM4', 'CHRM5', 'CHRNA1', 'CHRNA10', 'CHRNA2', 'CHRNA3', 'CHRNA4', 'CHRNA5',
                                'CHRNA6', 'CHRNA7', 'CHRNA9', 'CHRNB1', 'CHRNB2', 'CHRNB3', 'CHRNB4', 'CHRND', 'CHRNE', 'CHRNG'), 
        fill.by = "ident", stack = TRUE, same.y.lims = TRUE, flip = TRUE)  + NoLegend()


## VEN markers
Idents(VEN4clust) <- "M_type"

VEN4clust <- RenameIdents(VEN4clust, `cVEN` = "VEN", `PC` = "non_VEN", `TRI` = "non_VEN")
VEN4clust$VEN_type <- VEN4clust@active.ident
Idents(VEN4clust) <- "VEN_type"

VEN_markers.cad_1 <- FindMarkers(VEN4clust, ident.1 = "VEN", ident.2 = "non_VEN", only.pos = TRUE, min.pct = 1, logfc.threshold = 0.5)
VEN_markers.cad_2 <- cosg(VEN4clust, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                          expressed_pct = 1, n_genes_user = 1939)
VEN_markers.cad_1 <- subset(VEN_markers.cad_1, rownames(VEN_markers.cad_1) %in% VEN_markers.cad_2$names$VEN)

VEN_Negmarkers.cad_1 <- FindMarkers(VEN4clust, ident.1 = "VEN", ident.2 = "non_VEN", only.pos = FALSE, min.pct = 0.8, logfc.threshold = 0.5)
VEN_Negmarkers.cad_1 <- subset(VEN_Negmarkers.cad_1, avg_log2FC < 0)
VEN_Negmarkers.cad_1 <- subset(VEN_Negmarkers.cad_1, pct.2 == 1)

VEN_Negmarkers.cad_2 <- cosg(VEN4clust, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                             expressed_pct = 1, n_genes_user = 1279)
VEN_Negmarkers.cad_1 <- subset(VEN_Negmarkers.cad_1, rownames(VEN_Negmarkers.cad_1) %in% VEN_Negmarkers.cad_2$names$non_VEN)

VEN_markers.cad <- rbind(VEN_markers.cad_1, VEN_Negmarkers.cad_1)

EnhancedVolcano(VEN_markers.cad, lab = rownames(VEN_markers.cad), x = 'avg_log2FC', y = 'p_val', pCutoff = 0.05, shape = c(4, 7, 2, 1),
                colAlpha = 1, selectLab = c("COL5A2", "EPHA6", "SLC35F2", "LYPD1", "SLC5A8", "PCDH9", "SLC7A14", "NPTXR", "DNAJA4", "HSPA1A"),
                labSize = 3, labFace = "italic", boxedLabels = TRUE, title = "VEN markers", subtitle = "vs. non-VEN",
                cutoffLineCol = "black", gridlines.minor = F, gridlines.major = F, FCcutoff = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-4,3, by = 1), limits = c(-4,3)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,9, by = 3), limits = c(0,9)) +
  theme(axis.text.x = element_text(size = 16, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        axis.line = element_line(color = "black", size = 1))

VEN_markers.cad_flit <- subset(VEN_markers.cad, p_val < 0.05)
VEN_markers.cad_flit <- subset(VEN_markers.cad_flit, abs(avg_log2FC) > 1)

write.table(VEN_markers.cad_flit, file = "VENvsNonVEN_markers.cad_flit.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VEN_markers.cad_ETIT <- subset(VEN_markers.cad_flit, rownames(VEN_markers.cad_flit) %in% ETIT_markers)

VEN_known_markers <- c("POU3F1", "ITGA4", "BMP3", "VAT1L", "SULF2", "LYPD1", "CHST8", "ADRA1A", "GABRQ", "SLC18A2", "DRD3",
                       "HTR2B", "NMB", "DISC1", "ATF3", "IL4R", "AVPR1A")

VEN_markers.cad_known <- subset(VEN_markers.cad_flit, rownames(VEN_markers.cad_flit) %in% VEN_known_markers)

Idents(VEN4clust) <- "M_type"
VEN4clust$M_type <- factor(x = VEN4clust$M_type, levels = c("cVEN", "PC", "TRI"))

RidgePlot(VEN4clust, features = rownames(VEN_markers.cad_ETIT), stack = TRUE, fill.by = "ident")
RidgePlot(VEN4clust, features = c("LYPD1", "SULF2", "ATF3", "ADRA1A", "CHST8"), stack = TRUE, fill.by = "ident")

VEN_markers_biotype <- subset(all_genes, all_genes$gene_name %in% rownames(VEN_markers.cad_flit))
VEN_markers_biotype <- as.data.frame(table(VEN_markers_biotype$gene_biotype))
colnames(VEN_markers_biotype) <- c("biotype", "count")

ggplot() + geom_bar(data = VEN_markers_biotype, aes(x = biotype, y = count), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(), axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 10), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,400, by = 80), limits = c(0,400))

write.table(VEN_markers_biotype, file = "VEN_markers_biotype.txt", col.names = TRUE, sep = "\t", quote = FALSE)


## VEN vs IT-PC
Idents(VEN4clust) <- "seurat_clusters"
VEN4clust$clust_PM <- paste(Idents(VEN4clust), VEN4clust$PM_type, sep = " ")
Idents(VEN4clust) <- "clust_PM"
DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "clust_PM")

VENIT <- subset(VEN4clust, ident = c("0 IT_PC", "0 ET_VEN", "0 IT_VEN", "1 ET_VEN", "1 IT_VEN"))
DimPlot(VENIT, reduction = "umap", label = FALSE, pt.size = 3, group.by = "clust_PM")
Idents(VENIT) <- "VEN_type"
DimPlot(VENIT, reduction = "umap", label = FALSE, pt.size = 3, group.by = "VEN_type")

VENvsITPC.cad_1 <- FindMarkers(VENIT, ident.1 = "VEN", ident.2 = "non_VEN", only.pos = TRUE, min.pct = 1, logfc.threshold = 0.5)
VENvsITPC.cad_2 <- cosg(VENIT, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                         expressed_pct = 1, n_genes_user = 1597)
VENvsITPC.cad_1 <- subset(VENvsITPC.cad_1, rownames(VENvsITPC.cad_1) %in% VENvsITPC.cad_2$names$VEN)

ITPCvsVEN.cad_1 <- FindMarkers(VENIT, ident.1 = "VEN", ident.2 = "non_VEN", only.pos = FALSE, min.pct = 0.8, logfc.threshold = 0.5)
ITPCvsVEN.cad_1 <- subset(ITPCvsVEN.cad_1, avg_log2FC < 0)
ITPCvsVEN.cad_1 <- subset(ITPCvsVEN.cad_1, pct.2 == 1)

ITPCvsVEN.cad_2 <- cosg(VENIT, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T, expressed_pct = 1,
                         n_genes_user = 1644)
ITPCvsVEN.cad_1 <- subset(ITPCvsVEN.cad_1, rownames(ITPCvsVEN.cad_1) %in% ITPCvsVEN.cad_2$names$non_VEN)

VENvsITPC_markers.cad <- rbind(VENvsITPC.cad_1, ITPCvsVEN.cad_1)

EnhancedVolcano(VENvsITPC_markers.cad, lab = rownames(VENvsITPC_markers.cad), x = 'avg_log2FC', y = 'p_val', pCutoff = 0.05, 
                shape = c(4, 7, 2, 1), colAlpha = 1, labSize = 3, labFace = "italic", boxedLabels = TRUE,
                selectLab = c("CCK", "NPTXR", "NPTX1", "NCALD", "CCL4L2", "CCL3L1", "KCTD16", "POU3F1", "ADRA1A", "SULF2", "PCP4L1",
                              "VAT1L", "SYT2", "COL5A2"), title = "VEN markers", subtitle = "vs. IT_PC", cutoffLineCol = "black",
                gridlines.minor = F, gridlines.major = F, FCcutoff = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-8,6, by = 2), limits = c(-8,6)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,12, by = 3), limits = c(0,12)) +
  theme(axis.text.x = element_text(size = 16, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        axis.line = element_line(color = "black", size = 1))

VENvsITPC_markers.cad_filt <- subset(VENvsITPC_markers.cad, p_val < 0.05)
VENvsITPC_markers.cad_filt <- subset(VENvsITPC_markers.cad_filt, abs(avg_log2FC) > 1)

write.table(VENvsITPC_markers.cad_filt, file = "VENvsITPC_markers.cad_filt.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VENvsITPC_markers_biotype <- subset(all_genes, all_genes$gene_name %in% rownames(VENvsITPC_markers.cad_filt))
VENvsITPC_markers_biotype <- as.data.frame(table(VENvsITPC_markers_biotype$gene_biotype))
colnames(VENvsITPC_markers_biotype) <- c("biotype", "count")  

ggplot() + geom_bar(data = VENvsITPC_markers_biotype, aes(x = biotype, y = count), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(), axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 10), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,600, by = 120), limits = c(0,600))

write.table(VENvsITPC_markers_biotype, file = "VENvsITPC_markers_biotype.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VENETIT <- subset(VEN4clust, ident = c("0 IT_PC", "0 ET_VEN", "0 IT_VEN", "1 ET_VEN", "1 IT_VEN", "1 ET_PC"))
Idents(VENETIT) <- "PM_type"
DimPlot(VENETIT, reduction = "umap", label = FALSE, pt.size = 3, group.by = "PM_type")

VENETIT <- RenameIdents(VENETIT, `ET_VEN` = "VEN", `IT_VEN` = "VEN")
VENETIT@meta.data$VEN_M_type <- VENETIT@active.ident
DimPlot(VENETIT, reduction = "umap", label = FALSE, pt.size = 3)

VENvsIT_markers.cad_ETIT <- subset(VENvsITPC_markers.cad_filt, rownames(VENvsITPC_markers.cad_filt) %in% ETIT_markers)
VENvsIT_markers.cad_known <- subset(VENvsITPC_markers.cad_filt, rownames(VENvsITPC_markers.cad_filt) %in% VEN_known_markers)

RidgePlot(VENETIT, features = c("POU3F1", "SYT2", "SPARC", "FEZF2", "NPTX1", "RNF152"), stack = TRUE, fill.by = "ident")
RidgePlot(VENETIT, features = c("SULF2", "ADRA1A", "LYPD1", "VAT1L", "CHST8", "BMP3"), stack = TRUE, fill.by = "ident")


## VEN vs ET-PC
DimPlot(VEN4clust, reduction = "umap", label = FALSE, pt.size = 3, group.by = "clust_PM")
VENET <- subset(VEN4clust, ident = c("0 ET_VEN", "0 IT_VEN",  "1 ET_PC", "1 ET_VEN", "1 IT_VEN"))
DimPlot(VENET, reduction = "umap", label = FALSE, pt.size = 3, group.by = "clust_PM")
Idents(VENET) <- "M_type"
DimPlot(VENET, reduction = "umap", label = FALSE, pt.size = 3, group.by = "M_type")

VENvsET.cad_1 <- FindMarkers(VENET, ident.1 = "cVEN", ident.2 = "PC", only.pos = TRUE, min.pct = 1, logfc.threshold = 0.5)
VENvsET.cad_2 <- cosg(VENET, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                         expressed_pct = 1, n_genes_user = 5017)
VENvsET.cad_1 <- subset(VENvsET.cad_1, rownames(VENvsET.cad_1) %in% VENvsET.cad_2$names$cVEN)

ETvsVEN.cad_1 <- FindMarkers(VENET, ident.1 = "cVEN", ident.2 = "PC", only.pos = FALSE, min.pct = 0.7, logfc.threshold = 0.5)
ETvsVEN.cad_1 <- subset(ETvsVEN.cad_1, avg_log2FC < 0)
ETvsVEN.cad_1 <- subset(ETvsVEN.cad_1, pct.2 == 1)

ETvsVEN.cad_2 <- cosg(VENET, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                         expressed_pct = 1, n_genes_user = 1747)
ETvsVEN.cad_1 <- subset(ETvsVEN.cad_1, rownames(ETvsVEN.cad_1) %in% ETvsVEN.cad_2$names$PC)

VENvsET_markers.cad <- rbind(VENvsET.cad_1, ETvsVEN.cad_1)

EnhancedVolcano(VENvsET_markers.cad, lab = rownames(VENvsET_markers.cad), x = 'avg_log2FC', y = 'p_val', pCutoff = 0.05, 
                shape = c(4, 7, 2, 1), colAlpha = 1, labSize = 3, labFace = "italic", boxedLabels = TRUE,
                selectLab = c("SYNJ1", "DIRAS1", "PEG10", "SYT2", "SCN4B", "LINC02002", "AL512590.1", "CCK", "CARD16", "CLLU1"),
                title = "VEN markers", subtitle = "vs. ET_PC", cutoffLineCol = "black",
                gridlines.minor = F, gridlines.major = F, FCcutoff = 1) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(-3,3, by = 1), limits = c(-3,3)) + 
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,8, by = 2), limits = c(0,8)) +
  theme(axis.text.x = element_text(size = 16, colour = "black"), axis.text.y = element_text(size = 16, colour = "black"),
        axis.line = element_line(color = "black", size = 1))

VENvsET_markers.cad_filt <- subset(VENvsET_markers.cad, p_val < 0.05)
VENvsET_markers.cad_filt <- subset(VENvsET_markers.cad_filt, abs(avg_log2FC) > 1)
write.table(VENvsET_markers.cad_filt, file = "VENvsET_markers.cad_filt.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VENvsET_markers_biotype <- subset(all_genes, all_genes$gene_name %in% rownames(VENvsET_markers.cad_filt))
VENvsET_markers_biotype <- as.data.frame(table(VENvsET_markers_biotype$gene_biotype))
colnames(VENvsET_markers_biotype) <- c("biotype", "count")  

ggplot() + geom_bar(data = VENvsET_markers_biotype, aes(x = biotype, y = count), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(), axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 10), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,1000, by = 200), limits = c(0,1000))

write.table(VENvsET_markers_biotype, file = "VENvsET_markers_biotype.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VENvsET_markers.cad_ETIT <- subset(VENvsET_markers.cad_filt, rownames(VENvsET_markers.cad_filt) %in% ETIT_markers)
VENvsET_markers.cad_known <- subset(VENvsET_markers.cad_filt, rownames(VENvsET_markers.cad_filt) %in% VEN_known_markers)

RidgePlot(VENETIT, features = c("AVPR1A", "CCK", "LINC02002", "AL512590.1", "CARD16", "RXFP2"), stack = TRUE, fill.by = "ident")
RidgePlot(VENETIT, features = c("SCN4B", "ETV5", "AKAP12", "SYNJ1", "DIRAS1", "LONRF2"), stack = TRUE, fill.by = "ident")

VENvsET.cad_1_chr <- subset(VENvsET.cad_1, p_val < 0.05)
VENvsET.cad_1_chr <- subset(VENvsET.cad_1_chr, avg_log2FC > 1)

VENvsET.cad_1_chr <- rownames(VENvsET.cad_1_chr)

VENvsET.cad_1_chr.df <- bitr(VENvsET.cad_1_chr, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
## 21.97% of input gene IDs are fail to map...

VENvsET.cad_1_ID <- unique(VENvsET.cad_1_chr.df$ENTREZID)

VENvsET.cad_1_ego <- enrichGO(gene = VENvsET.cad_1_ID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)

VENvsET.cad_1_kegg <- enrichKEGG(gene = VENvsET.cad_1_ID, keyType = "kegg", organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2)

ETvsVEN.cad_1_chr <- subset(ETvsVEN.cad_1, p_val < 0.05)
ETvsVEN.cad_1_chr <- subset(ETvsVEN.cad_1_chr, avg_log2FC < -1)

ETvsVEN.cad_1_chr <- rownames(ETvsVEN.cad_1_chr)

ETvsVEN.cad_1_chr.df <- bitr(ETvsVEN.cad_1_chr, fromType = "SYMBOL", toType = c("ENTREZID", "ENSEMBL"), OrgDb = org.Hs.eg.db)
## 3.01% of input gene IDs are fail to map...

ETvsVEN.cad_1_ID <- unique(ETvsVEN.cad_1_chr.df$ENTREZID)

ETvsVEN.cad_1_ego <- enrichGO(gene = ETvsVEN.cad_1_ID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
ETvsVEN.cad_1_ego <- clusterProfiler::simplify(ETvsVEN.cad_1_ego, cutoff = 0.7, by = "p.adjust", select_fun = min)

goCls <- pairwise_termsim(ETvsVEN.cad_1_ego, method = "Resnik", semData = d)
treeplot(goCls)

ETvsVEN.cad_1_kegg <- enrichKEGG(gene = ETvsVEN.cad_1_ID, keyType = "kegg", organism = 'hsa', pAdjustMethod = "BH", pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.2)
clusterProfiler::dotplot(ETvsVEN.cad_1_kegg, title = "ETvsVEN.cad_1_kegg")


## VEN subtypes
Idents(VEN4clust) <- "M_type"
VEN <- subset(VEN4clust, ident = "cVEN")
Idents(VEN) <- "seurat_clusters"

VEN <- RenameIdents(VEN, `0` = "subtype_1", `1` = "subtype_2")
VEN$subtype <- VEN@active.ident
Idents(VEN) <- "subtype"
VEN$subtype <- factor(x = VEN$subtype, levels = c("subtype_1", "subtype_2"))
DimPlot(VEN, reduction = "umap", label = TRUE, pt.size = 4) + NoLegend()

Idents(VEN) <- "sample_ID"
DimPlot(VEN, reduction = "umap", label = TRUE, pt.size = 4, shape.by = "subtype") + NoLegend()

Idents(VEN) <- "subtype"
VEN_submarkers <- FindAllMarkers(VEN, only.pos = TRUE, min.pct = 1, logfc.threshold = 0.5)
VEN_submarkers_1 <- subset(VEN_submarkers, cluster == "subtype_1")
VEN_submarkers_2 <- subset(VEN_submarkers, cluster == "subtype_2")

VEN_sub1markers <- FindMarkers(VEN, ident.1 = "subtype_1", ident.2 = "subtype_2", only.pos = FALSE, logfc.threshold = 0.5)
VEN_sub1markers <- subset(VEN_sub1markers, avg_log2FC < 0)
VEN_sub1markers <- subset(VEN_sub1markers, rownames(VEN_sub1markers) %in% VEN_submarkers_2$gene)
VEN_sub1markers$cluster <- "subtype_2"
VEN_sub1markers$gene <- rownames(VEN_sub1markers)

VEN_submarkers_cosg <- cosg(VEN, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                             expressed_pct = 1, n_genes_user = 268)

VEN_submarkers_1 <- subset(VEN_submarkers_1, rownames(VEN_submarkers_1) %in% VEN_submarkers_cosg$names$subtype_1)

VEN_submarkers_cosg <- cosg(VEN, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                            expressed_pct = 1, n_genes_user = 1187)

VEN_sub1markers <- subset(VEN_sub1markers, rownames(VEN_sub1markers) %in% VEN_submarkers_cosg$names$subtype_2)

VEN_submarkers_1_biotype <- subset(all_genes, all_genes$gene_name %in% rownames(VEN_submarkers_1))
VEN_submarkers_1_biotype <- as.data.frame(table(VEN_submarkers_1_biotype$gene_biotype))
colnames(VEN_submarkers_1_biotype) <- c("biotype", "count") 

VEN_submarkers_2_biotype <- subset(all_genes, all_genes$gene_name %in% rownames(VEN_sub1markers))
VEN_submarkers_2_biotype <- as.data.frame(table(VEN_submarkers_2_biotype$gene_biotype))
colnames(VEN_submarkers_2_biotype) <- c("biotype", "count") 

write.table(VEN_submarkers_1_biotype, file = "VEN_submarkers_1_biotype.txt", col.names = TRUE, sep = "\t", quote = FALSE)
write.table(VEN_submarkers_2_biotype, file = "VEN_submarkers_2_biotype.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VEN_submarkers_filt <- rbind(VEN_submarkers_1, VEN_sub1markers)
write.table(VEN_submarkers_filt, file = "VEN_submarkers_filt.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VEN_submarkers_mean <- AverageExpression(VEN, features = VEN_submarkers_filt$gene, assays = "RNA", slot = "counts")
VEN_submarkers_mean <- as.data.frame(VEN_submarkers_mean)
VEN_submarkers_mean$group <- VEN_submarkers_filt[match(rownames(VEN_submarkers_mean), VEN_submarkers_filt$gene), 6]
VEN_submarkers_mean$L2FC  <- VEN_submarkers_filt[match(rownames(VEN_submarkers_mean), VEN_submarkers_filt$gene), 2]

ggplot(VEN_submarkers_mean, aes(x = RNA.subtype_1, y = RNA.subtype_2, color = L2FC)) + geom_point(size = 3) + 
  scale_color_gradient2(low = "blue", mid = "lightgrey", high = "#d7191c", midpoint = 0, name = "L2FC") +
  geom_abline(intercept = -1, slope = 1, col = 'black', linetype = 'dashed', size = 0.8) +
  geom_abline(intercept = 1, slope = 1, col = 'black', linetype = 'dashed', size = 0.8) +
  geom_abline(intercept = 0, slope = 1, col = 'black', linetype = 'dashed', size = 0.8) +
  scale_x_continuous(expand = c(0, 0), breaks = seq(0,12, by = 3), limits = c(0,12)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,12, by = 3), limits = c(0,12)) +
  theme(axis.ticks.length = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.5), panel.background = element_rect(fill = "white"))

RidgePlot(VEN, features = c("NEUROD6", "CACNA1H", "MIAT", "SNURF", "KCNQ5-AS1", "GRIK2"), stack = TRUE, fill.by = "ident")
RidgePlot(VEN, features = c("OR7E21P", "HLX-AS1", "MIR4791", "RN7SL826P", "MIR4276", "PEX12P1"), stack = TRUE, fill.by = "ident")

DoHeatmap(VEN, features = ETIT_markers, draw.lines = FALSE, hjust = 0) + scale_fill_gradientn(colors = c("#2f58a7", "lightgrey", "#FF6347"))

Idents(VENETIT) <- "clust_PM"
DimPlot(VENETIT, reduction = "umap", label = FALSE, pt.size = 3)

VENETIT <- RenameIdents(VENETIT, `0 IT_PC` = "ITPC", `1 ET_PC` = "ETPC", `1 ET_VEN` = "S2VEN", `1 IT_VEN` = "S2VEN",
                               `0 IT_VEN` = "S1VEN", `0 ET_VEN` = "S1VEN")

VENETIT@meta.data$CM_type <- VENETIT@active.ident
VENETIT$CM_type <- factor(x = VENETIT$CM_type , levels = c("ITPC", "ETPC", "S1VEN", "S2VEN"))
Idents(VENETIT) <- "CM_type"

DimPlot(VENETIT, reduction = "umap", pt.size = 3)

RidgePlot(VENETIT, features = c("POU3F1", "SULF2", "ADRA1A", "VAT1L", "LYPD1", "AVPR1A", "CHST8", "ITGA4"), stack = TRUE, fill.by = "ident")
RidgePlot(VENETIT, features = c("BMP3", "SLC18A2", "DRD3", "HTR2B", "NMB", "DISC1", "ATF3", "IL4R"), stack = TRUE, fill.by = "ident")

AverageExp <- AverageExpression(VENETIT, features = VENETIT@assays$RNA@var.features, assays = "RNA", slot = "counts")
typeof(AverageExp)
head(AverageExp$RNA)
coorda <- corr.test(AverageExp$RNA, AverageExp$RNA, method = "spearman")
pheatmap(coorda$r, color = c("#2f58a7", "lightgrey", "#FF6347"))


## extract data for M-type feature selecting
Idents(VEN4FS) <- "cell_name"

M_level2_62cells <- c('PC17', 'TRI31', 'PC25', 'PC26', 'PC27', 'PC05', 'TRI09', 'TRI23', 'cVEN01', 'cVEN02', 'cVEN03', 'cVEN04',
                      'cVEN05', 'cVEN06', 'cVEN07', 'cVEN08', 'cVEN14', 'PC37', 'PC38', 'PC39', 'PC40', 'PC41', 'PC42', 'cVEN09',
                      'cVEN10', 'cVEN11', 'PC33', 'PC34', 'cVEN12', 'cVEN13', 'PC36', 'PC44', 'PC45', 'PC46', 'PC47', 'PC48',
                      'PC49', 'PC50', 'PC51', 'PC53', 'PC54', 'PC55', 'PC56', 'PC57', 'PC58', 'PC59', 'PC60', 'PC61', 'PC62',
                      'PC63', 'cVEN15', 'cVEN16', 'cVEN17', 'cVEN18', 'cVEN19', 'cVEN20', 'PC65', 'PC66', 'cVEN21', 'PC67',
                      'PC68', 'cVEN22')

VEN4Mtype_L2 <- subset(VEN4FS, idents = M_level2_62cells)
Idents(VEN4Mtype_L2) <- "M_type"

VEN4Mtype_L2 <- RenameIdents(VEN4Mtype_L2, `cVEN` = "VEN", `PC` = "non_VEN", `TRI` = "non_VEN")
VEN4Mtype_L2$VEN_type <- VEN4Mtype_L2@active.ident
Idents(VEN4Mtype_L2) <- "VEN_type"

Mlevel2_markers.cad_1 <- FindMarkers(VEN4Mtype_L2, ident.1 = "VEN", ident.2 = "non_VEN", only.pos = TRUE, min.pct = 1, logfc.threshold = 0.5)
Mlevel2_markers.cad_1_chr <- rownames(Mlevel2_markers.cad_1)

Mlevel2_markers.cad_2 <- cosg(VEN4Mtype_L2, groups = "all", assay = "RNA", slot = "data", mu = 1, remove_lowly_expressed = T,
                              expressed_pct = 1, n_genes_user = 4200)
Mlevel2_markers.cad_2_chr <- Mlevel2_markers.cad_2$names[,1]

Mlevel2_markers.cad <- generics ::intersect(Mlevel2_markers.cad_1_chr, Mlevel2_markers.cad_2_chr)

Mlevel2_markers <- subset(all_genes, all_genes$gene_name %in% Mlevel2_markers.cad)
write.table(Mlevel2_markers, file = "62cells_Mlevel2_markers.txt", col.names = TRUE, sep = "\t", quote = FALSE)

Mlevel2_markers_biotype <- as.data.frame(table(Mlevel2_markers$gene_biotype))
colnames(Mlevel2_markers_biotype) <- c("biotype", "count")  

ggplot() + geom_bar(data = Mlevel2_markers_biotype, aes(x = biotype, y = count), position = position_dodge2(padding = 0.3), stat = "identity") + 
  theme(axis.title = element_blank(), axis.text.x = element_text(size = 6, angle = 45, vjust = 1, hjust = 1),
        axis.ticks.x = element_blank(), axis.ticks.length.y = unit(.25, "cm"), axis.line = element_line(color = "black", size = 0.4), 
        axis.text.y = element_text(size = 10), panel.background = element_rect(fill = "white")) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0,2500, by = 500), limits = c(0,2500))

VEN4Mtype_L2_mat <- as.data.frame(VEN4Mtype_L2@assays$RNA@scale.data)
VEN4Mtype_L2_mat <- VEN4Mtype_L2_mat[Mlevel2_markers.cad,]
pheatmap(VEN4Mtype_L2_mat, fontsize_row = 1, clustering_method = "ward.D2", fontsize_col = 5) 

VEN4clust_mat <- as.data.frame(VEN4clust@assays$RNA@scale.data)
M_typeFS_38genes <- c("GRIA4", "SYTL2", "NPTX1", "RXFP1", "MMD", "SPOCK1", "AC025287.3", "DNAJA4", "FOS", "SLITRK1", "ATP8A1",
                      "MAP2K1", "NUDT4", "SLC7A14", "PCDH9", "CAMK1D", "OPTN", "RNF41", "INPP5A", "FAM126B", "PLEKHB2", "MEIS2",
                      "POU3F1", "ADRA1A", "SULF2", "COL21A1", "WNT5A", "C6orf163", "SETBP1", "FAM111A-DT", "AL512590.1",
                      "HCG17", "PLA2G12B", "SLC27A5", "CARD16", "LINC01845", "AC092634.5", "SLC7A2")

M_typeFS_38genes_mat <- VEN4clust_mat[M_typeFS_38genes,]
pheatmap(M_typeFS_38genes_mat, fontsize_row = 5, clustering_method = "ward.D", fontsize_col = 5, color = c("#2f58a7", "lightgrey", "#FF6347"))



# hdWGCNA and TF
theme_set(theme_cowplot())

VENnetwork <- SetupForWGCNA(VENETIT, gene_select = "variable", wgcna_name = "VENnetwork")
length(VENnetwork@misc$VENnetwork$wgcna_genes)

VENnetwork <- MetacellsByGroups(seurat_obj = VENnetwork, group.by = "VEN_M_type", reduction = "pca", ident.group = "VEN_M_type",
                                min_cells = 22, k = 15)

VENnetwork <- NormalizeMetacells(VENnetwork)

VENnetwork <- SetDatExpr(VENnetwork, group_name = "VEN", group.by = "VEN_M_type", assay = "RNA", slot = "data")
VENnetwork <- TestSoftPowers(VENnetwork, networkType = "signed")
plot_list <- PlotSoftPowers(VENnetwork)
wrap_plots(plot_list, ncol = 2)

power_table <- GetPowerTable(VENnetwork)
head(power_table)


## construct co-expression network
VENnetwork <- ConstructNetwork(VENnetwork, soft_power = 12, setDatExpr = FALSE, tom_name = "VEN")
PlotDendrogram(VENnetwork, main = "VEN hdWGCNA Dendrogram")

VENnetwork@misc$VENnetwork$wgcna_modules %>% head
write.table(VENnetwork@misc$VENnetwork$wgcna_modules, file = "VEN_modules.txt", col.names = TRUE, sep = "\t", quote = FALSE)

table(VENnetwork@misc$VENnetwork$wgcna_modules$module)
## green    turquoise    purple     brown    grey         midnightblue    black    red       salmon     cyan 
## 113      369          69         292      494          57              92       99        58         58 
## grey60   blue         pink       yellow   greenyellow  lightcyan       tan      magenta 
## 51       292          89         115      67           54              61       70

Tom <- GetTOM(VENnetwork)


## Compute harmonized module eigengenes
VENnetwork <- ScaleData(VENnetwork, features =  VariableFeatures(VENnetwork))
VENnetwork <- ModuleEigengenes(VENnetwork, reduction.use = "pca")


## harmonized module eigengenes
hMEs <- GetMEs(VENnetwork)
head(hMEs)
MEs <- GetMEs(VENnetwork, harmonized = FALSE)
head(MEs)


## Compute module connectivity
VENnetwork <- ModuleConnectivity(VENnetwork, group.by = "VEN_M_type", group_name = "VEN")
VENnetwork <- ResetModuleNames(VENnetwork, new_name = "VEN-M")

PlotKMEs(VENnetwork, ncol = 5)

modules <- GetModules(VENnetwork)

hub_df <- GetHubGenes(VENnetwork, n_hubs = 10)
head(hub_df)
write.table(hub_df, file = "hub_genes.txt", col.names = TRUE, sep = "\t", quote = FALSE)

VENnetwork <- ModuleExprScore(VENnetwork, n_genes = 25, method = "Seurat")


## Visualization
plot_list <- ModuleFeaturePlot(VENnetwork, features = "hMEs", order = TRUE, point_size = 1)
wrap_plots(plot_list, ncol = 5) 

plot_list <- ModuleFeaturePlot(VENnetwork, features = "scores", order = "shuffle", ucell = TRUE, point_size = 1)
wrap_plots(plot_list, ncol = 5)

MEs <- GetMEs(VENnetwork, harmonized = TRUE)
mods <- colnames(MEs); mods <- mods[mods != "grey"]

VENnetwork@meta.data <- cbind(VENnetwork@meta.data, MEs)
P <- DotPlot(VENnetwork, features = mods, group.by = "VEN_M_type")
P <- P + coord_flip() + RotatedAxis() + scale_color_gradient2(high = "red", mid = "grey95", low = "blue")
P

theme_set(theme_cowplot())
ModuleNetworkPlot(VENnetwork)

HubGeneNetworkPlot(VENnetwork, n_hubs = 3, n_other = 5, edge_prop = 0.75, mods = 'all')

hubFS_genes <- intersect(hub_df$gene_name, M_typeFS_38genes)

VENnetwork <- RunModuleUMAP(VENnetwork, n_hubs = 10, n_neighbors = 15, min_dist = 0.1)
ModuleUMAPPlot(VENnetwork, edge.alpha = 0.5, sample_edges = TRUE, edge_prop = 0.1, label_hubs = 1,
               keep_grey_edges = FALSE, vertex.label.cex = 0.1, label_genes = c("ADRA1A", "SULF2"))


## Enrichment analysis
set.seed(12345)
dbs <- c('GO_Biological_Process_2023','GO_Cellular_Component_2023','GO_Molecular_Function_2023')

VENnetwork <- RunEnrichr(VENnetwork, dbs = dbs, max_genes = Inf)
enrich_df <- GetEnrichrTable(VENnetwork)

EnrichrDotPlot(VENnetwork, mods = "all", database = "GO_Biological_Process_2023", n_terms = 2, term_size = 8, p_adj = TRUE) +
  scale_color_stepsn(colors = rev(viridis::magma(256)))

EnrichrBarPlot(VENnetwork, outdir = "enrichr_plots", n_terms = 10, plot_size = c(5,7), logscale = TRUE)


## TF
theme_set(theme_cowplot())
set.seed(12345)

pfm_core <- TFBSTools::getMatrixSet(x = JASPAR2020,
                                    opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

VENnetwork <- MotifScan(VENnetwork, species_genome = 'hg38', pfm = pfm_core, EnsDb = EnsDb.Hsapiens.v86)
motif_df <- GetMotifs(VENnetwork)

tf_genes <- unique(motif_df$gene_name)
modules <- GetModules(VENnetwork)
nongrey_genes <- subset(modules, module != 'grey') %>% .$gene_name
genes_use <- c(tf_genes, nongrey_genes)

VENnetwork <- SetWGCNAGenes(VENnetwork, genes_use)
VENnetwork <- SetDatExpr(VENnetwork, group.by = 'VEN_M_type', group_name = "VEN", assay = "RNA")

model_params <- list(objective = 'reg:squarederror', max_depth = 1, eta = 0.1, nthread = 16, alpha = 0.5)

VENnetwork <- ConstructTFNetwork(VENnetwork, model_params)

results <- GetTFNetwork(VENnetwork)
head(results)

VENnetwork <- AssignTFRegulons(VENnetwork, strategy = "A", reg_thresh = 0.01, n_tfs = 10)
RegulonBarPlot(VENnetwork, selected_tf = 'TBP')

VENnetwork <- RegulonScores(VENnetwork, target_type = 'positive', ncores = 8)
VENnetwork <- RegulonScores(VENnetwork, target_type = 'negative', cor_thresh = -0.05, ncores = 8)

pos_regulon_scores <- GetRegulonScores(VENnetwork, target_type = 'positive')
neg_regulon_scores <- GetRegulonScores(VENnetwork, target_type = 'negative')

tf_regulons <- GetTFRegulons(VENnetwork)
hub_df <- GetHubGenes(VENnetwork, n_hubs = 25) %>% subset(gene_name %in% tf_regulons$tf)

Idents(VENnetwork) <- "VEN_M_type"
marker_tfs <- FindAllMarkers(VENnetwork, features = unique(tf_regulons$tf))

top_tfs <- marker_tfs %>% subset(cluster == 'VEN') %>% slice_max(n = 25, order_by = avg_log2FC)
intersect(top_tfs$gene, hub_df$gene_name)

cur_tf <- 'FOXF2'

TFNetworkPlot(VENnetwork, selected_tfs = cur_tf)

group1 <- VENnetwork@meta.data %>% subset(VEN_M_type == "VEN") %>% rownames
group2 <- VENnetwork@meta.data %>% subset(VEN_M_type == "ET_PC") %>% rownames

dregs <- FindDifferentialRegulons(VENnetwork, barcodes1 = group1, barcodes2 = group2)

head(dregs)




# WD ^ ^ written by Shuxuan Lyu 2025/07/17




