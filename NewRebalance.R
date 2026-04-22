# this script reproduces Fig 6c,d,e,f,g,h,i,j from the EpiTrace paper + a rebalancing experiment 

# in the rebalancing experiment, EC/Peric. has only 31 cells but the cap is 350, 
# leaving an ~11x imbalance against the smallest major type. 
# Possible change is to lower the cap substantially (closer to 31–50) or explicitly exclude this cell type and note it.

# Sys.setenv(XML_CONFIG = "/share/pkg.8/libxml2/2.15.1/install/bin/xml2-config")
# pak::pkg_install('MagpiePKU/EpiTrace')


library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(ArchR)
library(dplyr)
library(tidyr)
library(readr)
library(GenomicRanges)
library(reshape2)
library(openxlsx)
library(ggplot2)
library(Matrix)
library(EpiTrace)
library(Seurat)
library(SeuratObject)
library(ggtree)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(ggbeeswarm)
library(ggthemes)
library(ggsci)
library(ggpubr)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(zoo)
library(RColorBrewer)
library(openxlsx)
library(WGCNA)
library(parallel)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(SingleCellExperiment) # added. required for line: assays(atac_1)[['counts']] -> mtx
library(easyLift) # added. required for this line: clock_combine_hg38 <- plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])) %>% easyLift::easyLiftOver('hg19_hg38')

options(bitmapType='cairo') # solve the PNG/X11 issue
data("clock_gr_list")
set.seed(1234)

# build Epitrace object 

# readRDS('./Packed_Figure_Reproducibility/input_data_for_EpiTrace/Figure6/matrix//Multiome_ATAC_SCE.RDS') -> atac_1
readRDS('/projectnb/ds596/students/jishan/input_data_for_EpiTrace/Figure6/Multiome_ATAC_SCE.RDS') -> atac_1

atac_1@rowRanges -> GR
assays(atac_1)[['counts']] -> mtx
atac_1@colData$Cell.ID -> cellname_vec

splicedm <- readRDS('/projectnb/ds596/students/jishan/input_data_for_EpiTrace/Figure6/Brain_scMultiome_spliced_mtx_for_CreateAssay.rds')

# Epitrace 
inputpeak <- GR
clock_combine_hg38 <- plyranges::reduce_ranges(c(clock_gr_list[[1]],clock_gr_list[[2]])) %>% easyLift::easyLiftOver('hg19_hg38')
initiated_peaks <- Init_Peakset(inputpeak)
initiated_peaks_df <- as.data.frame(initiated_peaks,row.names = NULL)
paste0(initiated_peaks_df$seqnames,'_',initiated_peaks_df$start,'_',initiated_peaks_df$end) -> initiated_peaks_df$peakName
initiated_mm <- Init_Matrix(cellname = cellname_vec,peakname = initiated_peaks_df$peakName,matrix = mtx)

epitrace_obj_age_conv_estimated <- EpiTraceAge_Convergence(initiated_peaks,initiated_mm,ref_genome = 'hg38',non_standard_clock = F,qualnum=1,Z_cutoff = 3,parallel = T,ncore_lim = 46,iterative_time=10, normalization_method = 'randomized')
epitrace_obj_age_conv_estimated[['all']] <- Seurat::CreateAssayObject(counts = initiated_mm[,epitrace_obj_age_conv_estimated$cell],min.cells = 0,min.features = 0,check.matrix = F)
DefaultAssay(epitrace_obj_age_conv_estimated) <- 'all'


# Add meta 
readRDS('/projectnb/ds596/students/jishan/input_data_for_EpiTrace/Figure6/meta/Figure6_Brain_scMultiome_meta.rds') -> cell_meta 

# cbind(epitrace_obj_age_conv_estimated@meta.data, temp22 %>% dplyr::select(Cluster.Name, celltype, cytotrace_rna)) -> epitrace_obj_age_conv_estimated@meta.data


# ----------------added the following bc seurat_clusters not found in seurat obj ------------------
# bc of error no surat obj at : Idents(epitrace_obj_age_estimated_multiome) <- epitrace_obj_age_estimated_multiome$seurat_clusters


# added the following code bc suerat obj was missing 
# Then pull seurat_clusters directly from atac_1 and join separately
atac_1@colData[, c("Cell.ID", "seurat_clusters")] %>% 
  as.data.frame() -> clusters_df

epitrace_obj_age_conv_estimated@meta.data <- left_join(
  epitrace_obj_age_conv_estimated@meta.data,
  clusters_df,
  by = c("cell" = "Cell.ID")
)

# Clean join replacing the broken cbind block
epitrace_obj_age_conv_estimated@meta.data -> temp

# Drop any pre-existing duplicate columns
temp <- temp[, !colnames(temp) %in% c("celltype", "Cluster.Name", "cytotrace_rna")]

# Join properly by cell ID
temp22 <- left_join(
  temp,
  cell_meta[, c("cell", "celltype", "Cluster.Name", "cytotrace_rna")],
  by = "cell"
)

# Verify join worked before attaching
print(table(temp22$celltype))  # should show real cell types, not "unlabeled"

# Re-attach
rownames(temp22) <- temp22$cell
epitrace_obj_age_conv_estimated@meta.data <- temp22

# replaced bc meta.data rownames are "1", "2", "3"... but Seurat expects them to be the actual cell barcodes. Fix it:
# Set meta.data rownames to the actual cell barcodes
rownames(epitrace_obj_age_conv_estimated@meta.data) <- as.character(epitrace_obj_age_conv_estimated$cell)

# Now try again
epitrace_obj_age_conv_estimated[['rna_spliced']] <- Seurat::CreateAssayObject(
  counts = splicedm[, as.character(epitrace_obj_age_conv_estimated$cell)],
  min.cells = 0,
  min.features = 0,
  check.matrix = F
)


epitrace_obj_age_conv_estimated -> epitrace_obj_age_estimated_multiome


# added lsi 
epitrace_obj_age_estimated_multiome <- RunTFIDF(epitrace_obj_age_estimated_multiome)
epitrace_obj_age_estimated_multiome <- FindTopFeatures(epitrace_obj_age_estimated_multiome, min.cutoff = 'q0')
epitrace_obj_age_estimated_multiome <- RunSVD(epitrace_obj_age_estimated_multiome)
names(epitrace_obj_age_estimated_multiome@reductions)
# =====================================================

epitrace_obj_age_estimated_multiome <- RunUMAP(
  epitrace_obj_age_estimated_multiome,
  reduction = 'lsi',
  dims = 2:30
)

# Save UMAP  (can ignore this UMAP; not part of the paper)
pdf('/projectnb/ds596/students/jishan/Plots/UMAP_clusters.pdf', height = 6, width = 7)
DimPlot(epitrace_obj_age_estimated_multiome, reduction = 'umap', 
        group.by = 'seurat_clusters', label = TRUE, pt.size = 0.1) +
  ggtitle('UMAP - Seurat Clusters')
dev.off()

# Figure 6c and 6d
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

Idents(epitrace_obj_age_estimated_multiome) <- epitrace_obj_age_estimated_multiome$seurat_clusters
DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'all'

# This is just a Seurat v5 API change, "slot" was renamed to "layer"
Signac::CreateChromatinAssay(counts = GetAssayData(epitrace_obj_age_estimated_multiome, layer = 'data'),min.cells = 0,min.features = 0,genome = BSgenome.Hsapiens.UCSC.hg38@seqinfo) -> epitrace_obj_age_estimated_multiome[['all_chrom_assay']]




DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'all_chrom_assay'

epitrace_obj_age_estimated_multiome <- AddMotifs(
  object = epitrace_obj_age_estimated_multiome,
  genome = BSgenome.Hsapiens.UCSC.hg38,
  pfm = pfm
)

# ex of changes made: changed "3" into "c3", "0" to "c0" bc those are the name in levels(Idents(epitrace_obj_age_estimated_multiome)):
# [1] "c10" "c2"  "c0"  "c4"  "c11" "c1"  "c9"  "c5"  "c6"  "c3"  "c7"  "c8"  "c13" "c12"

da_peaks <- FindMarkers(
  object = epitrace_obj_age_estimated_multiome,
  ident.1 = 'c3',
  ident.2 = 'c0',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
open.peaks <- AccessiblePeaks(epitrace_obj_age_estimated_multiome, idents = c("c3", "c0"))

# match the overall GC content in the peak set
meta.feature <- GetAssayData(epitrace_obj_age_estimated_multiome, assay = "all_chrom_assay", layer = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)

enriched.motifs <- FindMotifs(
  object = epitrace_obj_age_estimated_multiome,
  features = top.da.peak,
  background=peaks.matched
)

epitrace_obj_age_estimated_multiome <- RunChromVAR(
  object = epitrace_obj_age_estimated_multiome,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'chromvar'

differential.activity <- FindMarkers(
  object = epitrace_obj_age_estimated_multiome,
  ident.1 = c('c0','c11'),
  ident.2 = c('c3','c8'),
  only.pos = F,
  mean.fxn = rowMeans,
  fc.name = "avg_diff",features = epitrace_obj_age_estimated_multiome@assays$chromvar@data %>% rownames,
  logfc.threshold = 0,min.pct = 0,min.diff.pct = 0
)

differential.activity <- FindMarkers(
  object = epitrace_obj_age_estimated_multiome,
  ident.1 = c('c0'),
  ident.2 = c('c3'),
  only.pos = F,
  test.use='LR',
  fc.name = "avg_diff",features = epitrace_obj_age_estimated_multiome@assays$chromvar@data %>% rownames,
  logfc.threshold = 0,min.pct = 0,min.diff.pct = 0
)

differential.activity$motif <- rownames(differential.activity)
differential.activity.with.enrichment <- left_join(differential.activity,enriched.motifs)
differential.activity.with.enrichment$gene_symbol <- gsub('\\(.+','',differential.activity.with.enrichment$motif.name) %>% toupper()
differential.activity.with.enrichment$motif_origin <- ifelse(differential.activity.with.enrichment$gene_symbol == gsub('\\(.+','',differential.activity.with.enrichment$motif.name), 'human','others')
differential.activity.with.enrichment <- differential.activity.with.enrichment[!grepl(':',differential.activity.with.enrichment$gene_symbol),]
c(differential.activity.with.enrichment$gene_symbol) %>%unique -> da_tf_names

DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'rna_spliced'
epitrace_obj_age_estimated_multiome <- NormalizeData(epitrace_obj_age_estimated_multiome)
epitrace_obj_age_estimated_multiome <- ScaleData(epitrace_obj_age_estimated_multiome)

subset(epitrace_obj_age_estimated_multiome,seurat_clusters %in% c("c3","c0","c11","c9","c5","c8")) -> tt1
tt1$class <- ifelse(tt1$seurat_clusters %in% c('c0','c11','c5','c9'),'NR2F1+','LMO3+')
Idents(tt1) <- tt1$class

tt1@assays$rna_spliced@data@Dimnames[[1]] -> features_name
da_tf_names[da_tf_names %in% features_name] -> sel_features
DefaultAssay(tt1) <- 'rna_spliced'
FindAllMarkers(tt1,features = sel_features, logfc.threshold = 0,min.pct=0,min.cells.feature = 0,min.cells.group = 0,min.diff.pct = 0,only.pos = F) -> mrk_3
mrk_3$gene_symbol <- mrk_3$gene
da_tf_with_activity <- left_join(mrk_3 %>% dplyr::select(avg_log2FC,cluster,gene_symbol,p_val_adj),differential.activity.with.enrichment %>% dplyr::select(avg_diff,p_val_adj,motif,gene_symbol,motif_origin),by='gene_symbol')

da_tf_with_activity_2 <- lapply(unique(da_tf_with_activity$gene_symbol),function(x){
  da_tf_with_activity[da_tf_with_activity$gene_symbol %in% x & da_tf_with_activity$cluster %in% 'NR2F1+',] %>% arrange(p_val_adj.x,p_val_adj.y) %>% head(1)
}) %>% bind_rows()

# Plot selected TF pairs
sel_gene <- c('TCF4','NR2F1')
ggplot(da_tf_with_activity_2[da_tf_with_activity_2$motif_origin %in% 'human',],aes(x=-log(p_val_adj.x,10),y=-log(p_val_adj.y,10))) + geom_point() + ggrepel::geom_label_repel(data=da_tf_with_activity_2[-log(da_tf_with_activity_2$p_val_adj.x,10)>3 & -log(da_tf_with_activity_2$p_val_adj.y,10)>100,],mapping=aes(label=gene_symbol)) + theme_classic() + xlab('-logP Gene Exp') + ylab('-logP TF activity')+ theme(text=element_text(size=10)) + ggtitle('Statistics for TF activity and Gene Exp') -> p1

DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'chromvar'
FeaturePlot(epitrace_obj_age_estimated_multiome,features=da_tf_with_activity_2[da_tf_with_activity_2$gene_symbol %in% sel_gene,]$motif,min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1,combine = F) -> pl1

lapply(c(1:length(pl1)), function(x){
  pl1[[x]]$data -> xdata
  colnames(xdata)[4] -> motif_id
  colnames(xdata)[4] <- "Motif_Activity"
  xdata$Motif <- motif_id
  xdata$SYMBOL <- da_tf_with_activity_2$gene_symbol[da_tf_with_activity_2$motif %in% motif_id] %>% unique
  xdata$cell <- rownames(xdata)
  rownames(xdata) <- NULL
  xdata
}) %>% bind_rows() -> TFBS_activity_on_umap_plotdf

TFBS_activity_on_umap_plotdf[TFBS_activity_on_umap_plotdf$SYMBOL %in% "NR2F1", ] -> NR2F1_tfbs_plot
NR2F1_tfbs_plot <- arrange(NR2F1_tfbs_plot,Motif_Activity)

DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'rna_spliced'
FeaturePlot(epitrace_obj_age_estimated_multiome,features=sel_gene,min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1,combine = F) -> pl1_exp

DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'chromvar'
pl2 <- lapply(pl1,function(x){
  x$data -> xdata
  colnames(xdata)[4] <- 'value'
  x$labels$title -> featureid
  feature_new_id <- paste0(da_tf_with_activity_2$gene_symbol[da_tf_with_activity_2$motif %in% featureid]%>%unique,' TFBS_ChrAcc ',featureid)
  xdata$cell <- rownames(xdata)
  xdata <- left_join(xdata,tt1@meta.data %>% as.data.frame,by='cell')
  xdata <- arrange(xdata,value)
  ggplot(xdata,aes(x=umap_1,y=umap_2,fill=value))+geom_point(size=1.5,pch=21,color='darkgray')+xlim(c(-10.2,4.5)) + ylim(c(-9,4)) + scale_fill_gradientn(colors=c('cornflowerblue','beige','red')) + theme_void() + ggtitle(feature_new_id) + theme(text=element_text(size=10))
})

pl3 <- list(p1,pl2[[1]])

pdf('/projectnb/ds596/students/jishan/Plots//Figure6c_6d_NR2F1_TCF4_TF_activity.pdf',height=4,width=9)
wrap_plots(pl3,ncol = 2) + plot_layout(guides = 'collect')
dev.off()
#


# Figure 6e and 6f
DefaultAssay(epitrace_obj_age_estimated_multiome) <- 'rna_spliced'
epitrace_obj_age_estimated_multiome <- NormalizeData(epitrace_obj_age_estimated_multiome)
epitrace_obj_age_estimated_multiome <- ScaleData(epitrace_obj_age_estimated_multiome)

FeaturePlot(epitrace_obj_age_estimated_multiome,features=c('NR2F1','LMO3'),reduction = 'umap',order = T,pt.size = 0.1,cols = c('beige','red'),combine = F) -> pl3
pl4 <- lapply(pl3,function(x){
  x$data -> xdata
  colnames(xdata)[4] <- 'value'
  x$labels$title -> featureid
  xdata <- arrange(xdata,value)
  ggplot(xdata,aes(x=umap_1,y=umap_2,fill=value))+geom_point(size=1.5,pch=21,color='darkgray')+xlim(c(-10.2,4.5)) + ylim(c(-9,4)) + scale_fill_gradientn(colors=c('cornflowerblue','beige','red')) + theme_void() + ggtitle(featureid) + theme(text=element_text(size=20))
})

pdf('/projectnb/ds596/students/jishan/Plots/Figure6e_6f_Brain_NR2F1_LMO3_RNA_Expression.pdf',height=4,width=8)
wrap_plots(pl4,nrow=1,byrow = T) + plot_layout(guides = 'collect')
dev.off()



# Drop broken cytotrace_rna
epitrace_obj_age_estimated_multiome@meta.data$cytotrace_rna <- NULL

# Join from cell_meta which has the real values
meta_temp <- epitrace_obj_age_estimated_multiome@meta.data
meta_temp <- left_join(meta_temp, cell_meta[, c("cell", "cytotrace_rna")], by = "cell")

# Verify
sum(!is.na(meta_temp$cytotrace_rna))  # should be ~8981

# Re-attach
rownames(meta_temp) <- meta_temp$cell
epitrace_obj_age_estimated_multiome@meta.data <- meta_temp

# Verify on the object
sum(!is.na(epitrace_obj_age_estimated_multiome$cytotrace_rna))
summary(epitrace_obj_age_estimated_multiome$cytotrace_rna)

# Figure 6g and 6h
subset(epitrace_obj_age_estimated_multiome,seurat_clusters %in% c("c3","c0","c11","c9","c5","c8")) -> tt1
tt1$class <- ifelse(tt1$seurat_clusters %in% c("c0","c11","c5","c9"),'NR2F1+','LMO3+')
tt1@meta.data[tt1$seurat_clusters %in% c("c0","c3"),] -> nIPC_cytotrace_Epitrace_df

my_comparision <- list(c('LMO3+','NR2F1+'))
ggplot(nIPC_cytotrace_Epitrace_df,aes(x=class,y=EpiTraceAge_iterative)) + geom_violin(aes(fill=class)) + ggbeeswarm::geom_beeswarm(size=3,pch=21) + geom_boxplot(outlier.alpha = 0,width=0.3,fill='gray') +  ggpubr::stat_compare_means(comparisons =my_comparision,label = 'p.signif' ) + theme_classic()  + theme(text=element_text(size=20),axis.title.x = element_blank()) + ylab('EpiTrace Age') + stat_compare_means(label.y = 1.2,size=5) + ggsci::scale_fill_jco() + ggtitle('Age') -> p1

my_comparision <- list(c('LMO3+','NR2F1+'))
ggplot(nIPC_cytotrace_Epitrace_df,aes(x=class,y=cytotrace_rna)) + geom_violin(aes(fill=class)) + ggbeeswarm::geom_beeswarm(size=3,pch=21) + geom_boxplot(outlier.alpha = 0,width=0.3,fill='gray') +  ggpubr::stat_compare_means(comparisons =my_comparision,label = 'p.signif' ) + theme_classic()  + theme(text=element_text(size=20),axis.title.x = element_blank()) + ylab('CytoTRACE by RNA') + stat_compare_means(label.y = 1.2,size=5) + ggsci::scale_fill_jco()  + ggtitle('Stemness') -> p2
(p1|p2) + plot_layout(guides='collect')

pdf('/projectnb/ds596/students/jishan/Plots/Figure6g_6h_compare_Epitrace_and_Cytotrace_two_nIPC_population.pdf',height=4.5,width=9)
(p1|p2) + plot_layout(guides='collect')
dev.off()

# fix for 6i:

# Check celltype values match exactly
print(unique(epitrace_obj_age_estimated_multiome$celltype))

#: convert to factor first
epitrace_obj_age_estimated_multiome$celltype <- factor(epitrace_obj_age_estimated_multiome$celltype)
print(levels(epitrace_obj_age_estimated_multiome$celltype))

# Figure 6i
brain_meta <- epitrace_obj_age_estimated_multiome@meta.data

colorlist <- c('GluN5' = "cadetblue4",
               'IN1'="#A2CD5A",
               'nIPC/GluN1' = 'cornflowerblue',
               'IN2' = "#BCEE68",
               'SP' = "darkgreen",
               'GluN2' = 'cadetblue1',
               'IN3'= "#CAFF70",
               'RG' = "red",
               'GluN4' = 'cadetblue3',
               'GluN3'= 'cadetblue2',
               'Cyc. Prog.' = "#FFA500",
               'mGPC/OPC'= "#68228B",
               'EC/Peric.' = "#8B0000")


ggplot(brain_meta,aes(y=celltype,x=EpiTraceAge_iterative)) + geom_violin(scale='width', aes(fill = celltype), alpha = 0.8) + geom_boxplot(width=0.15, outlier.alpha = 0, fill = 'black') + theme_classic() + theme(text=element_text(size=20), legend.position = "none") + scale_fill_manual(values = colorlist) + xlab("\n EpiTrace") + ylab("Cell type") -> pp1
ggplot(brain_meta,aes(y=celltype,x=cytotrace_rna)) + geom_violin(scale='width', aes(fill = celltype), alpha = 0.8) + geom_boxplot(width=0.15, outlier.alpha = 0, fill = 'black') + theme_classic() + theme(text=element_text(size=20), axis.title.y = element_blank(), axis.text.y = element_blank(), legend.position = "none") + scale_fill_manual(values = colorlist) + scale_x_reverse() + xlab("\n CytoTrace") -> pp2

pdf('/projectnb/ds596/students/jishan/Plots/Figure6i_Brain_scMultiome_Epitrace_Cytotrace.pdf',height=9,width=8)
pp1|pp2
dev.off()



# Figure 6j
colors_for_plot <- colorRampPalette(RColorBrewer::brewer.pal(9,'Set1'))(length(levels(epitrace_obj_age_estimated_multiome$celltype )))
names(colors_for_plot) <- levels(epitrace_obj_age_estimated_multiome$celltype)
colors2 <- c('red','orange','cornflowerblue','cadetblue1','cadetblue2','cadetblue3','cadetblue4','darkgreen','darkorchid4','darkolivegreen1','darkolivegreen2','darkolivegreen3','darkred')
names(colors2) <- names(colors_for_plot)

Idents(epitrace_obj_age_estimated_multiome) <- epitrace_obj_age_estimated_multiome$celltype
phylores_multiome_exN <- RunEpiTracePhylogeny(subset(epitrace_obj_age_estimated_multiome,subset=celltype %in% c('RG','nIPC/GluN1','GluN2','GluN3','GluN4','GluN5','SP')))

phylores_multiome_exN[['iterative']][[2]] -> data.tree_clock
data.tree_clock <- ape::root(data.tree_clock,outgroup='RG')
tree_plot_clock <- ggtree::ggtree(data.tree_clock,
                                  layout = "rectangular", ladderize = FALSE) +
  geom_tiplab(aes(), color='black', size = 5, offset = 40) + geom_tippoint(aes(fill=label,color=label),size=7) +
  scale_color_manual(values = colors2)
xmax <- (tree_plot_clock$data$branch.length %>% max(na.rm = T)) *
  1.4
tree_plot_clock <- tree_plot_clock + xlim(c(NA, xmax))

pdf('/projectnb/ds596/students/jishan/Plots/Figure6j_Brain_Epitrace_tree.pdf',width=6,height=3)
print(tree_plot_clock)
dev.off()




# ======== add following to check for imbalanced data ============
gini_coeff <- function(counts) {
  counts <- sort(counts[counts > 0])
  n <- length(counts)
  2 * sum(seq_len(n) * counts) / (n * sum(counts)) - (n + 1) / n
}

ct_counts <- table(epitrace_obj_age_estimated_multiome$celltype)
ct_df     <- as.data.frame(ct_counts, stringsAsFactors = FALSE)
colnames(ct_df) <- c("celltype", "n")
ct_df$pct <- round(100 * ct_df$n / sum(ct_df$n), 2)
ct_df      <- ct_df[order(-ct_df$n), ]

gini  <- gini_coeff(ct_df$n)
ratio <- max(ct_df$n) / min(ct_df$n)

print(ct_df, row.names = FALSE)
message(sprintf("Gini  : %.3f  (flag if > 0.4)", gini))
message(sprintf("Ratio : %.1fx (flag if > 5x)",  ratio))
# ======  added abovve code to check for imbalance data =========

rownames(epitrace_obj_age_estimated_multiome@meta.data) <- as.character(epitrace_obj_age_estimated_multiome$cell)

epitrace_balanced <- resample_cells(
  epitrace_obj_age_estimated_multiome,
  alpha = 0.7,
  mode  = "down"
)

# ── RE-RUN EpiTrace on balanced dataset ────────────────────────────────────

cells_keep <- epitrace_balanced$cell
mtx_balanced <- initiated_mm[, cells_keep]

epitrace_obj_age_balanced <- EpiTraceAge_Convergence(
  initiated_peaks,
  mtx_balanced,
  ref_genome           = 'hg38',
  non_standard_clock   = FALSE,
  qualnum              = 1,
  Z_cutoff             = 3,
  parallel             = TRUE,
  ncore_lim            = 46,
  iterative_time       = 10,
  normalization_method = 'randomized'
)

# Transfer metadata from the balanced Seurat object
meta_bal <- epitrace_balanced@meta.data[
  match(epitrace_obj_age_balanced$cell, epitrace_balanced$cell), ]
rownames(meta_bal) <- meta_bal$cell
epitrace_obj_age_balanced@meta.data <- meta_bal

message(sprintf("Balanced EpiTrace done: %d cells", ncol(epitrace_obj_age_balanced)))


# ── COMPARISON PLOTS: original vs balanced ──────────────────────────────────

library(patchwork)

# shared colour list (same as your original script)
colorlist <- c('GluN5'      = "cadetblue4",
               'IN1'        = "#A2CD5A",
               'nIPC/GluN1' = 'cornflowerblue',
               'IN2'        = "#BCEE68",
               'SP'         = "darkgreen",
               'GluN2'      = 'cadetblue1',
               'IN3'        = "#CAFF70",
               'RG'         = "red",
               'GluN4'      = 'cadetblue3',
               'GluN3'      = 'cadetblue2',
               'Cyc. Prog.' = "#FFA500",
               'mGPC/OPC'   = "#68228B",
               'EC/Peric.'  = "#8B0000")

# ── Scatter: original vs balanced age per shared cell ────────────────────

shared_cells <- intersect(epitrace_obj_age_estimated_multiome$cell,
                          epitrace_obj_age_balanced$cell)

age_compare <- data.frame(
  cell         = shared_cells,
  age_original = epitrace_obj_age_estimated_multiome@meta.data[
    match(shared_cells, epitrace_obj_age_estimated_multiome$cell),
    "EpiTraceAge_iterative"],
  age_balanced = epitrace_obj_age_balanced@meta.data[
    match(shared_cells, epitrace_obj_age_balanced$cell),
    "EpiTraceAge_iterative"],
  celltype     = epitrace_obj_age_estimated_multiome@meta.data[
    match(shared_cells, epitrace_obj_age_estimated_multiome$cell),
    "celltype"]
)

cor_val <- cor(age_compare$age_original, age_compare$age_balanced,
               use = "complete.obs")
message(sprintf("Pearson r (original vs balanced): %.4f", cor_val))

p_scatter <- ggplot(age_compare,
                    aes(x = age_original, y = age_balanced, colour = celltype)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "gray40") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5, size = 3.5,
           label = sprintf("r = %.4f", cor_val)) +
  scale_colour_manual(values = colorlist) +
  labs(title = "EpiTrace age: original vs balanced",
       x = "EpiTraceAge (original)", y = "EpiTraceAge (balanced)",
       colour = NULL) +
  theme_classic() +
  theme(text = element_text(size = 11))

# ── Figure 6i equivalent: violin per cell type, both runs side by side ───

# original
brain_meta_orig <- epitrace_obj_age_estimated_multiome@meta.data
brain_meta_orig$run <- "Original"

# balanced — need cytotrace_rna too, pull from original meta by cell ID
brain_meta_bal <- epitrace_obj_age_balanced@meta.data
brain_meta_bal$cytotrace_rna <- epitrace_obj_age_estimated_multiome@meta.data[
  match(brain_meta_bal$cell, epitrace_obj_age_estimated_multiome$cell),
  "cytotrace_rna"]
brain_meta_bal$run <- "Balanced"

cols_needed <- c("cell", "celltype", "EpiTraceAge_iterative", "cytotrace_rna", "run")
brain_combined <- rbind(brain_meta_orig[, cols_needed],
                        brain_meta_bal[, cols_needed])
brain_combined$run <- factor(brain_combined$run, levels = c("Original", "Balanced"))

pp1_orig <- ggplot(brain_meta_orig,
                   aes(y = celltype, x = EpiTraceAge_iterative)) +
  geom_violin(scale = 'width', aes(fill = celltype), alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.alpha = 0, fill = 'black') +
  scale_fill_manual(values = colorlist) +
  scale_y_discrete(limits = names(colorlist)) +
  theme_classic() +
  theme(text = element_text(size = 14), legend.position = "none") +
  xlab("EpiTrace (original)") + ylab("Cell type")

pp1_bal <- ggplot(brain_meta_bal,
                  aes(y = celltype, x = EpiTraceAge_iterative)) +
  geom_violin(scale = 'width', aes(fill = celltype), alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.alpha = 0, fill = 'black') +
  scale_fill_manual(values = colorlist) +
  scale_y_discrete(limits = names(colorlist)) +
  theme_classic() +
  theme(text = element_text(size = 14),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y  = element_blank()) +
  xlab("EpiTrace (balanced)")

pp2_orig <- ggplot(brain_meta_orig,
                   aes(y = celltype, x = cytotrace_rna)) +
  geom_violin(scale = 'width', aes(fill = celltype), alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.alpha = 0, fill = 'black') +
  scale_fill_manual(values = colorlist) +
  scale_y_discrete(limits = names(colorlist)) +
  scale_x_reverse() +
  theme_classic() +
  theme(text = element_text(size = 14),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y  = element_blank()) +
  xlab("CytoTRACE (original)")

pp2_bal <- ggplot(brain_meta_bal,
                  aes(y = celltype, x = cytotrace_rna)) +
  geom_violin(scale = 'width', aes(fill = celltype), alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.alpha = 0, fill = 'black') +
  scale_fill_manual(values = colorlist) +
  scale_y_discrete(limits = names(colorlist)) +
  scale_x_reverse() +
  theme_classic() +
  theme(text = element_text(size = 14),
        legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y  = element_blank()) +
  xlab("CytoTRACE (balanced)")

pdf('/projectnb/ds596/students/jishan/Plots/Comparison_Figure6i_Original_vs_Balanced.pdf',
    height = 9, width = 16)
pp1_orig | pp1_bal | pp2_orig | pp2_bal
dev.off()

# ── Figure 6g/6h equivalent: nIPC NR2F1+ vs LMO3+ ───────────────────────

subset(epitrace_obj_age_balanced,
       seurat_clusters %in% c("c3","c0","c11","c9","c5","c8")) -> tt1_bal
tt1_bal$class <- ifelse(tt1_bal$seurat_clusters %in% c("c0","c11","c5","c9"),
                        'NR2F1+', 'LMO3+')

# pull cytotrace from original meta
tt1_bal$cytotrace_rna <- epitrace_obj_age_estimated_multiome@meta.data[
  match(tt1_bal$cell, epitrace_obj_age_estimated_multiome$cell),
  "cytotrace_rna"]

tt1_bal@meta.data[tt1_bal$seurat_clusters %in% c("c0","c3"), ] -> nIPC_bal_df

my_comparision <- list(c('LMO3+', 'NR2F1+'))

make_gh_plot <- function(df, y_col, y_label, title_label) {
  ggplot(df, aes(x = class, y = .data[[y_col]])) +
    geom_violin(aes(fill = class)) +
    ggbeeswarm::geom_beeswarm(size = 3, pch = 21) +
    geom_boxplot(outlier.alpha = 0, width = 0.3, fill = 'gray') +
    ggpubr::stat_compare_means(comparisons = my_comparision, label = 'p.signif') +
    stat_compare_means(label.y = 1.2, size = 5) +
    ggsci::scale_fill_jco() +
    theme_classic() +
    theme(text = element_text(size = 20), axis.title.x = element_blank()) +
    ylab(y_label) +
    ggtitle(title_label)
}

# original (replot cleanly for direct visual comparison)
nIPC_orig_df <- epitrace_obj_age_estimated_multiome@meta.data[
  epitrace_obj_age_estimated_multiome$seurat_clusters %in% c("c0","c3"), ]
nIPC_orig_df$class <- ifelse(nIPC_orig_df$seurat_clusters == "c0", 'NR2F1+', 'LMO3+')

p_age_orig <- make_gh_plot(nIPC_orig_df, "EpiTraceAge_iterative",
                           "EpiTrace Age", "Age — Original")
p_cyt_orig <- make_gh_plot(nIPC_orig_df, "cytotrace_rna",
                           "CytoTRACE by RNA", "Stemness — Original")
p_age_bal  <- make_gh_plot(nIPC_bal_df,  "EpiTraceAge_iterative",
                           "EpiTrace Age", "Age — Balanced")
p_cyt_bal  <- make_gh_plot(nIPC_bal_df,  "cytotrace_rna",
                           "CytoTRACE by RNA", "Stemness — Balanced")

pdf('/projectnb/ds596/students/jishan/Plots/Comparison_Figure6gh_Original_vs_Balanced.pdf',
    height = 5, width = 18)
(p_age_orig | p_cyt_orig | p_age_bal | p_cyt_bal) + plot_layout(guides = 'collect')
dev.off()

# ── Figure 6j equivalent: phylogeny tree on balanced object ───────────────

# need celltype set as Idents and assay available
epitrace_obj_age_balanced[['all']] <- Seurat::CreateAssayObject(
  counts      = initiated_mm[, epitrace_obj_age_balanced$cell],
  min.cells   = 0,
  min.features = 0,
  check.matrix = FALSE
)
DefaultAssay(epitrace_obj_age_balanced) <- 'all'
Idents(epitrace_obj_age_balanced) <- epitrace_obj_age_balanced$celltype

colors2 <- c('red','orange','cornflowerblue','cadetblue1','cadetblue2',
             'cadetblue3','cadetblue4','darkgreen','darkorchid4',
             'darkolivegreen1','darkolivegreen2','darkolivegreen3','darkred')
names(colors2) <- names(colorlist)

phylores_bal <- RunEpiTracePhylogeny(
  subset(epitrace_obj_age_balanced,
         subset = celltype %in% c('RG','nIPC/GluN1','GluN2','GluN3',
                                  'GluN4','GluN5','SP'))
)

data.tree_bal <- phylores_bal[['iterative']][[2]]
data.tree_bal <- ape::root(data.tree_bal, outgroup = 'RG')

tree_plot_bal <- ggtree(data.tree_bal, layout = "rectangular", ladderize = FALSE) +
  geom_tiplab(color = 'black', size = 5, offset = 40) +
  geom_tippoint(aes(fill = label, color = label), size = 7) +
  scale_color_manual(values = colors2) +
  xlim(c(NA, max(ggtree(data.tree_bal)$data$branch.length, na.rm = TRUE) * 1.4))

pdf('/projectnb/ds596/students/jishan/Plots/Comparison_Figure6j_Balanced_Tree.pdf',
    width = 6, height = 3)
print(tree_plot_bal)
dev.off()




# ── Cell type counts: before vs after balancing ──────────────────────────

ct_before <- as.data.frame(table(epitrace_obj_age_estimated_multiome$celltype),
                           stringsAsFactors = FALSE)
colnames(ct_before) <- c("celltype", "n_before")

ct_after <- as.data.frame(table(epitrace_balanced$celltype),
                          stringsAsFactors = FALSE)
colnames(ct_after) <- c("celltype", "n_after")

ct_comparison <- merge(ct_before, ct_after, by = "celltype", all = TRUE)
ct_comparison$n_before[is.na(ct_comparison$n_before)] <- 0
ct_comparison$n_after[is.na(ct_comparison$n_after)]   <- 0
ct_comparison$cells_removed <- ct_comparison$n_before - ct_comparison$n_after
ct_comparison$pct_retained  <- round(100 * ct_comparison$n_after / ct_comparison$n_before, 1)

ct_comparison <- ct_comparison[order(-ct_comparison$n_before), ]

print(ct_comparison, row.names = FALSE)

# celltype n_before n_after cells_removed pct_retained
# nIPC/GluN1     2348     350          1998         14.9
# GluN2     1546     350          1196         22.6
# IN1      959     350           609         36.5
# GluN3      798     350           448         43.9
# IN2      780     350           430         44.9
# RG      646     350           296         54.2
# GluN4      459     350           109         76.3
# mGPC/OPC      359     350             9         97.5
# Cyc. Prog.      341     341             0        100.0
# IN3      301     301             0        100.0
# GluN5      223     223             0        100.0
# SP      190     190             0        100.0
# EC/Peric.       31      31             0        100.0

message(sprintf("Total before: %d", sum(ct_comparison$n_before)))
message(sprintf("Total after : %d", sum(ct_comparison$n_after)))

# Total before: 8981
# Total after : 3886

