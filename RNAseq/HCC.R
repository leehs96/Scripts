library(tidyverse)
library(EnhancedVolcano)
library(DESeq2)
library(ggrepel)
library(tximport)
library(RColorBrewer)
library(pheatmap)
library(pathfindR)
library(knitr)
library(reshape)
library(sva)
library(umapr)

#setwd('/users/hslee/HCC/HCC')

###################################################################################################

#..........Data pre-processing (batch effect correction, remove lowly expressed gene).............#

###################################################################################################

# bring RAW DATA

Raw_Data <- read.csv('merged_NOT_batch_corrected_HCC_data.csv', header = T, row.names = 1)

# 1. Remove lowly expressed gene (row sum = 0)

rm_low_EXP_gene_Data <- Raw_Data[rowSums(Raw_Data) != 0,]

# 2. Batch effect correction

# non-biological factors that have to be corrected ( different sources )
# biological factors that have to be conserved  (sample type (Cancer, normal), inflammation...etc) 

### a. bring condition matrix

condition <- read.csv('HCC_condition.csv',row.names = 1, header = T)

### b. matching order between data and condition 

before_batch <- rm_low_EXP_gene_Data[,rownames(condition)]

### c. assign batch information (0=PlosGenetics, 1=Nature communication, 2=Cancer Cell)
###    assign sample type as a group (normal, normal_LT. FA.....wiHCC)

batch <- condition$batch
group <- condition$type
cov_mat <- as.matrix(condition %>% dplyr::select(type, type3))

###d. batch correction using ComBat_seq

after_batch <- ComBat_seq(as.matrix(before_batch), batch=batch, covar_mod = cov_mat)
write.csv(after_batch,'merged_batch_corrected_HCC_data.csv')


###################################################################################################

#........................Differential expressed gene analysis using DESeq2........................#

###################################################################################################


# assign batch corrected data matrix
data <- read.csv('merged_batch_corrected_HCC_data.csv', header = T, row.names = 1)

# assign ensGene to GeneID matching information
GENE <- read.csv('ensGene_GeneID.csv', header = T, row.names = 1)
dev.off()

###################################################################################################
###########.............................!!!!IMPORTANT!!!!............................##############
#                                                                                                 #
## assign condition                                                                               #
Condition <- 'miFTC_vs_miHCC'          
#                    ###### this variable is used as plot file_name & DESeq2 parameter            #
#                    ###### Condition must be 'control_vs_case' form  ex: FA_vs_HA                #
###################################################################################################

#sample data wrangling 
adjusted <- data
sample <- read.csv('HCC_condition.csv',header = T, row.names = 1)

sample$condition <- factor(sample$condition) # condition


sample <- sample[sample$condition %in% c(gsub('.*vs_','',Condition),gsub('_vs_.*','',Condition)),]


adjusted <- adjusted[, colnames(adjusted) %in% base::rownames(sample)]
adjusted <- as.matrix(adjusted)
all(rownames(sample) == colnames(adjusted))
adjusted <- adjusted[, rownames(sample)]
all(rownames(sample) == colnames(adjusted))

sample$condition <- factor(sample$condition)

dds <- DESeqDataSetFromMatrix(countData = adjusted,
                              colData = sample,
                              design = ~ condition)


#remove lowly exp ( >=10 ). gene
keep <- rowSums(counts(dds)) >= 10

dds <- dds[keep,]

#create sub-directory of each condition if not available
dir.create(paste0('./plot/',Condition))


#export normalized exp values
ddsN <- estimateSizeFactors(dds)
ddsN <- estimateDispersions(ddsN)

Ndds <- counts(ddsN, normalized = T)
#write.csv(Ndds,paste0('./Norm/',Condition,".csv"))

#DEG
ddsDE <- DESeq(dds)

#results change p-adj sig value 0.05
res <- results(ddsDE, contrast = c("condition",gsub('.*vs_','',Condition), gsub('_vs_.*','',Condition)) , alpha = 0.05)


#order by p-adj
resOrdered <- as.data.frame(res[order(res$padj),])
write.csv(merge(resOrdered %>% dplyr::filter(padj <= 0.05),GENE,by = 0, all.x=T) %>% relocate(GENEid) %>% dplyr::rename(ensGene = Row.names) %>% arrange(padj) ,paste0('./DEG/',Condition,".csv")  , row.names = F)

#PCA
vsd <- vst(dds, blind = F)
pcaData <- plotPCA(vsd, intgroup = c("condition"), returnData = T)
percentVar <- round(100 * attr(pcaData, "percentVar"))

tiff(file=paste0('./plot/',Condition,'/',Condition,"_PCA.tiff"),width=10, height=10, units = 'in', res = 150)

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() + stat_ellipse(aes(fill = condition), geom = 'polygon', alpha = 0.05) + 
  geom_text(aes(label=rownames(pcaData)),hjust=.5, vjust=-1.5, size=2, color='black', fontface = 'bold')
dev.off()

#UMAP

# UMAP_DF <- as.matrix(t(before_batch))
# UMAP <- umap(log10(UMAP_DF+1), metric = 'manhattan')
# 
# UMAP <- cbind(UMAP,sample %>% dplyr::select(condition))
# 
# UMAP %>%
#   mutate(sample = gsub('_.*','',rownames(UMAP))) %>%
#   ggplot(aes(UMAP1, UMAP2, color = condition)) +
#   geom_point(size = 4) +
#   ggsci::scale_color_jco() +
#   geom_text_repel(label=rownames(UMAP_DF),hjust=.5, size=2.5, color='black', fontface = 'bold') +
#   ggtitle('UMAP : before batch correction ') +
#   theme(plot.title = element_text(color="black", size=17, face="bold"))# + facet_wrap(~ Distance, scales = 'free')
# 




#sample to sample dist.
sampleDist <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDist)
rownames(sampleDistMatrix) <- vsd$condition
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")))(255)

#heat-map sample sample distance
tiff(file=paste0('./plot/',Condition,'/',Condition,"_diff.tiff"),width=10, height=10, units = 'in', res = 150)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDist,
         clustering_distance_cols = sampleDist,
         cluster_cols = T,
         cluster_rows = T,
         border_color = NA,
         display_numbers = F,
         cutree_rows = 2,
         cutree_cols = 2,
         color = colors,
         show_rownames = T
)
dev.off()

#Volcano
resOrdered_GENEid <- merge(as.data.frame(resOrdered),GENE,by=0)
resOrdered_GENEid <- resOrdered_GENEid[order(resOrdered_GENEid$padj),]
resOrdered_vol <- resOrdered_GENEid %>% mutate(GENEid = ifelse(!resOrdered_GENEid$GENEid %in% resOrdered_GENEid$GENEid[1:10],' ',resOrdered_GENEid$GENEid))


FC <- 2
Q <- 0.05

keyvals <- rep('grey75', nrow(resOrdered_vol)) 
names(keyvals) <- rep('Not Sig', nrow(resOrdered_vol)) 

keyvals[which(abs(resOrdered_vol$log2FoldChange) > FC & resOrdered_vol$padj > Q)] <- 'grey50'
names(keyvals)[which(abs(resOrdered_vol$log2FoldChange) > FC & resOrdered_vol$padj > Q)] <- 'Log2foldchange'

keyvals[which(abs(resOrdered_vol$log2FoldChange) < FC & resOrdered_vol$padj < Q)] <- 'grey25'
names(keyvals)[which(abs(resOrdered_vol$log2FoldChange)  < FC & resOrdered_vol$padj < Q)] <- '-Log10Q'

keyvals[which(resOrdered_vol$log2FoldChange < -FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean > 100)] <- 'blue2'
names(keyvals)[which(resOrdered_vol$log2FoldChange  < -FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean > 100)] <- 'Down:basemean > 100'

keyvals[which(resOrdered_vol$log2FoldChange < -FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean <= 100)] <- 'royalblue'
names(keyvals)[which(resOrdered_vol$log2FoldChange  < -FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean <= 100)] <- 'Down:basemean < 100'

keyvals[which(resOrdered_vol$log2FoldChange > FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean > 100)] <- 'red3'
names(keyvals)[which(resOrdered_vol$log2FoldChange > FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean > 100)] <- 'UP:basemean > 100'

keyvals[which(resOrdered_vol$log2FoldChange > FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean <= 100)] <- 'pink'
names(keyvals)[which(resOrdered_vol$log2FoldChange > FC & resOrdered_vol$padj < Q & resOrdered_vol$baseMean <= 100)] <- 'UP:basemean < 100'

unique(keyvals)
unique(names(keyvals))

tiff(file=paste0('./plot/',Condition,'/',Condition,"_Volcano.tiff"),width=13, height=10, units = 'in', res = 150)
EnhancedVolcano(resOrdered_vol ,
                title = gsub('_',' ',Condition),
                lab = resOrdered_vol$GENEid,
                x = 'log2FoldChange',
                y = 'padj',
                labhjust = 0.5,
                selectLab = resOrdered_vol$GENEid[1:10],
                xlim = c(-(ceiling(max(abs(resOrdered_vol$log2FoldChange)))+1),ceiling(max(abs(resOrdered_vol$log2FoldChange)))+1),
                boxedLabels = F,
                xlab = bquote(~Log[2]~'(fold change)'),
                ylab = bquote(~-Log[10]~'adjusted'~italic(P)),
                pCutoff = Q,
                FCcutoff = FC,
                pointSize = 1.5,
                labSize = 3,
                labFace = 'bold',
                colCustom = keyvals,
                colAlpha = 0.75,
                legendPosition = 'bottom',
                legendLabSize = 13,
                legendIconSize = 5.0,
                drawConnectors = F,
                widthConnectors = 0.8,
                colConnectors = 'grey50',
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'partial',
                borderWidth = 1.5,
                borderColour = 'black',
                legendLabels = F)
dev.off()

#fold change plot

resOrdered_fold <- as.data.frame(resOrdered_GENEid)[as.data.frame(resOrdered_GENEid)$padj <= 0.05,] %>%
  arrange(log2FoldChange) %>% 
  na.omit()

TOP10_fold <- rbind(head(resOrdered_fold, n=10),tail(resOrdered_fold,n=10)) %>% 
  mutate(Trend = ifelse(log2FoldChange <0,'DOWN','UP'))

jpeg(file=paste0('./plot/',Condition,'/',Condition,"_FC.jpeg"),width=650, height=724)

TOP10_fold %>%
  ggplot(aes(x=log2FoldChange, y=factor(GENEid, levels = GENEid))) +
  geom_segment(aes(x = 0, y = factor(GENEid, levels = GENEid), xend = log2FoldChange, yend = GENEid, color = Trend)) +
  theme_bw() + 
  geom_point(aes(color = Trend, size = abs(padj)),alpha = 1,size=5) +
  xlab("Fold change") +
  ylab('Enriched gene : FDR < 0.05') +
  scale_colour_manual(name="Trend",labels = c("UP", "DOWN"), values = c("UP"="#d1495b", "DOWN"="#00798c")) +
  theme_minimal() +
  ggtitle(paste0(gsub('_',' ',Condition),'\nTOP10 UP & DOWN genes(Fold change)')) + 
  theme(axis.text=element_text(size=10.5),axis.title=element_text(size=16,face="bold"),text = element_text(size = 16, face = "bold"), legend.position = 'none')
dev.off()

###################################################################################################

#.................................pathway analysis using PathfindR................................#

###################################################################################################

term_gene_graph_manual <- function (result_df, num_terms = 11, layout = "stress", use_description = FALSE, node_size = "num_genes") {
  if (!is.numeric(num_terms) & !is.null(num_terms)) {
    stop("`num_terms` must either be numeric or NULL!")
  }
  if (!is.logical(use_description)) {
    stop("`use_description` must either be TRUE or FALSE!")
  }
  ID_column <- ifelse(use_description, "Term_Description",
                      "ID")
  val_node_size <- c("num_genes", "p_val")
  if (!node_size %in% val_node_size) {
    stop("`node_size` should be one of ", paste(dQuote(val_node_size),
                                                collapse = ", "))
  }
  if (!is.data.frame(result_df))
    stop("`result_df` should be a data frame")
  necessary_cols <- c(ID_column, "lowest_p", "Up_regulated",
                      "Down_regulated")
  if (!all(necessary_cols %in% colnames(result_df))) {
    stop(paste(c("All of", paste(necessary_cols, collapse = ", "),
                 "must be present in `results_df`!"), collapse = " "))
  }
  if (!is.null(num_terms)) {
    if (nrow(result_df) < num_terms) {
      num_terms <- NULL
    }
  }
  result_df <- result_df[order(result_df$lowest_p, decreasing = FALSE),
  ]
  if (!is.null(num_terms)) {
    result_df <- result_df[1:num_terms, ]
  }
  graph_df <- data.frame()
  for (i in base::seq_len(nrow(result_df))) {
    up_genes <- unlist(strsplit(result_df$Up_regulated[i],
                                ", "))
    down_genes <- unlist(strsplit(result_df$Down_regulated[i],
                                  ", "))
    genes <- c(up_genes, down_genes)
    for (gene in genes) {
      graph_df <- rbind(graph_df, data.frame(Term = result_df[i,
                                                              ID_column], Gene = gene))
    }
  }
  up_genes <- lapply(result_df$Up_regulated, function(x) unlist(strsplit(x,
                                                                         ", ")))
  up_genes <- unlist(up_genes)
  g <- igraph::graph_from_data_frame(graph_df, directed = FALSE)
  cond_term <- names(igraph::V(g)) %in% result_df[, ID_column]
  cond_up_gene <- names(igraph::V(g)) %in% up_genes
  igraph::V(g)$type <- ifelse(cond_term, "term", ifelse(cond_up_gene,
                                                        "up", "down"))
  if (node_size == "num_genes") {
    sizes <- igraph::degree(g)
    sizes <- ifelse(igraph::V(g)$type == "term", sizes, 2)
    size_label <- "# genes"
  }
  else {
    idx <- match(names(igraph::V(g)), result_df[, ID_column])
    sizes <- -log10(result_df$lowest_p[idx])
    sizes[is.na(sizes)] <- 2
    size_label <- "-log10(p)"
  }
  igraph::V(g)$size <- sizes
  igraph::V(g)$label.cex <- 0.15
  igraph::V(g)$frame.color <- "gray"
  igraph::V(g)$color <- ifelse(igraph::V(g)$type == "term",
                               "#E5D7BF", ifelse(igraph::V(g)$type == "up", "red",
                                                 "green"))
  p <- ggraph::ggraph(g, layout = layout)
  p <- p + ggraph::geom_edge_link(alpha = 0.8, colour = "darkgrey")
  p <- p + ggraph::geom_node_point(ggplot2::aes_(color = ~I(color),
                                                 size = ~size))
  p <- p + ggplot2::scale_size(range = c(5, 10), breaks = round(seq(round(min(igraph::V(g)$size)),
                                                                    round(max(igraph::V(g)$size)), length.out = 4)), name = size_label)
  p <- p + ggplot2::theme_void()
  p <- p + ggraph::geom_node_text(ggplot2::aes_(label = ~name),
                                  nudge_y = 0.2)
  p <- p + ggplot2::scale_colour_manual(values = unique(igraph::V(g)$color),
                                        name = NULL, labels = c("enriched term", "up-regulated gene",
                                                                "down-regulated gene"))
  if (is.null(num_terms)) {
    p <- p + ggplot2::ggtitle("Term-Gene Graph")
  }
  else {
    p <- p + ggplot2::ggtitle("Term-Gene Graph", subtitle = paste0(gsub('_',' ',Condition))
  }
  p <- p + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5),
                          plot.subtitle = ggplot2::element_text(hjust = 0.5))
  return(p)
}

# assign input data
input <- resOrdered_GENEid %>%
  mutate(Gene.symbol = GENEid) %>%
  dplyr::select(Gene.symbol,log2FoldChange, padj) %>%
  dplyr::filter(padj <= 0.05)

kable(head(input))



# using all DEG (FDR < 0.05)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_KEGG_Top10_pathway.tiff"),width=10, height=15, units = 'in', res = 150)
pathfindR_pathway_df <- run_pathfindR(
  input,
  output_dir = paste0('./pathway/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.01)

dev.off()


clustered_pathway <- cluster_enriched_terms(pathfindR_pathway_df, plot_clusters_graph = F, plot_dend = F)
clustered_pathway_top5 <- subset(clustered_pathway, Cluster %in% 1:5)


tiff(file=paste0('./plot/',Condition,'/',Condition,"_Top5_cluster.tiff"),width=10, height=15, units = 'in', res = 150)
enrichment_chart(clustered_pathway_top5, plot_by_cluster = TRUE)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_term_gene_heatmap.tiff"),width=20, height=20, units = 'in', res = 150)
term_gene_heatmap(result_df = pathfindR_pathway_df, genes_df = input)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_term_gene_graph.tiff"),width=15, height=15, units = 'in', res = 150)
term_gene_graph_manual(result_df = pathfindR_pathway_df, use_description = TRUE , node_size = 'p_val')
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_term_UpSet_plot.tiff"),width=20, height=20, units = 'in', res = 150)
UpSet_plot(result_df = pathfindR_pathway_df, genes_df = input, method = 'heatmap')
dev.off()


tiff(file=paste0('./plot/',Condition,'/',Condition,"_GO_BP_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_BP_df <- run_pathfindR(
  gene_sets = 'GO-BP',
  input,
  output_dir = paste0('./GO/BP/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_GO_CC_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_CC_df <- run_pathfindR(
  gene_sets = 'GO-CC',
  input,
  output_dir = paste0('./GO/CC/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_GO_MF_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_MF_df <- run_pathfindR(
  gene_sets = 'GO-MF',
  input,
  output_dir = paste0('./GO/MF/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_CM_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_CM_df <- run_pathfindR(
  gene_sets = 'cell_markers',
  input,
  output_dir = paste0('./cell_markers/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

dev.off()

clustered_BP <- cluster_enriched_terms(pathfindR_BP_df, plot_clusters_graph = F, plot_dend = F)
clustered_BP_top5 <- subset(clustered_BP, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_Top5_cluster_BP.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_BP_top5, plot_by_cluster = TRUE)
dev.off()

clustered_CC <- cluster_enriched_terms(pathfindR_CC_df, plot_clusters_graph = F, plot_dend = F)
clustered_CC_top5 <- subset(clustered_CC, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_Top5_cluster_CC.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_CC_top5, plot_by_cluster = TRUE)
dev.off()

clustered_MF <- cluster_enriched_terms(pathfindR_CC_df, plot_clusters_graph = F, plot_dend = F)
clustered_MF_top5 <- subset(clustered_MF, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_Top5_cluster_MF.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_MF_top5, plot_by_cluster = TRUE)
dev.off()

clustered_CM <- cluster_enriched_terms(pathfindR_CC_df, plot_clusters_graph = F, plot_dend = F)
clustered_CM_top5 <- subset(clustered_CM, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_Top5_cluster_CM.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_CM_top5, plot_by_cluster = TRUE)
dev.off()

dev.off()

########################################################################################################################################
# using Up regulated DEG (FDR < 0.05, Log2FC > 1)


dev.off()

input2 <- read.csv(paste0('DEG/',Condition,'.csv'), header = T) %>%
  mutate(Gene.symbol = ifelse(is.na(GENEid),ensGene,GENEid)) %>%
  dplyr::select(Gene.symbol,log2FoldChange, padj) %>%
  dplyr::filter(padj <= 0.05 & log2FoldChange > 1)


kable(head(input2))


tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_KEGG_Top10_pathway.tiff"),width=10, height=15, units = 'in', res = 150)
pathfindR_pathway_df_FC1 <- run_pathfindR(
  input2,
  output_dir = paste0('./pathway_FC1/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.01)

dev.off()


clustered_pathway_FC1 <- cluster_enriched_terms(pathfindR_pathway_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_pathway_top5_FC1 <- subset(clustered_pathway_FC1, Cluster %in% 1:5)


tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_Top5_cluster.tiff"),width=10, height=15, units = 'in', res = 150)
enrichment_chart(clustered_pathway_top5_FC1, plot_by_cluster = TRUE)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_term_gene_heatmap.tiff"),width=20, height=20, units = 'in', res = 150)
term_gene_heatmap(result_df = pathfindR_pathway_df_FC1, genes_df = input2)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_term_gene_graph.tiff"),width=15, height=15, units = 'in', res = 150)
term_gene_graph_manual(result_df = pathfindR_pathway_df_FC1, use_description = TRUE , node_size = 'p_val')
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_term_UpSet_plot.tiff"),width=20, height=20, units = 'in', res = 150)
UpSet_plot(result_df = pathfindR_pathway_df_FC1, genes_df = input2, method = 'heatmap')
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_GO_BP_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_BP_df_FC1 <- run_pathfindR(
  gene_sets = 'GO-BP',
  input2,
  output_dir = paste0('./GO/BP_FC1/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_GO_CC_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_CC_df_FC1 <- run_pathfindR(
  gene_sets = 'GO-CC',
  input2,
  output_dir = paste0('./GO/CC_FC1/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_GO_MF_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_MF_df_FC1 <- run_pathfindR(
  gene_sets = 'GO-MF',
  input2,
  output_dir = paste0('./GO/MF_FC1/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_CM_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_CM_df_FC1 <- run_pathfindR(
  gene_sets = 'cell_markers',
  input2,
  output_dir = paste0('./cell_markers_FC1/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

dev.off()

clustered_BP_FC1 <- cluster_enriched_terms(pathfindR_BP_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_BP_top5_FC1 <- subset(clustered_BP_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_Top5_cluster_BP.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_BP_top5_FC1, plot_by_cluster = TRUE)
dev.off()

clustered_CC_FC1 <- cluster_enriched_terms(pathfindR_CC_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_CC_top5_FC1 <- subset(clustered_CC_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_Top5_cluster_CC.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_CC_top5_FC1, plot_by_cluster = TRUE)
dev.off()

clustered_MF_FC1 <- cluster_enriched_terms(pathfindR_MF_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_MF_top5_FC1 <- subset(clustered_MF_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_Top5_cluster_MF.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_MF_top5_FC1, plot_by_cluster = TRUE)
dev.off()

clustered_CM_FC1 <- cluster_enriched_terms(pathfindR_CM_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_CM_top5_FC1 <- subset(clustered_CM_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_Top5_cluster_CM.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_CM_top5_FC1, plot_by_cluster = TRUE)
dev.off()

dev.off()
########################################################################################################################################
# using Up regulated DEG (FDR < 0.05, Log2FC < -1)


dev.off()

input3 <- read.csv(paste0('DEG/',Condition,'.csv'), header = T) %>%
  mutate(Gene.symbol = ifelse(is.na(GENEid),ensGene,GENEid)) %>%
  dplyr::select(Gene.symbol,log2FoldChange, padj) %>%
  dplyr::filter(padj <= 0.05 & log2FoldChange < -1)


kable(head(input2))


tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_KEGG_Top10_pathway.tiff"),width=10, height=15, units = 'in', res = 150)
pathfindR_pathway_df_FC1 <- run_pathfindR(
  input3,
  output_dir = paste0('./pathway_FC1_DOWN/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.01)

dev.off()


clustered_pathway_FC1 <- cluster_enriched_terms(pathfindR_pathway_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_pathway_top5_FC1 <- subset(clustered_pathway_FC1, Cluster %in% 1:5)


tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_Top5_cluster.tiff"),width=10, height=15, units = 'in', res = 150)
enrichment_chart(clustered_pathway_top5_FC1, plot_by_cluster = TRUE)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_term_gene_heatmap.tiff"),width=20, height=20, units = 'in', res = 150)
term_gene_heatmap(result_df = pathfindR_pathway_df_FC1, genes_df = input2)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_term_gene_graph.tiff"),width=15, height=15, units = 'in', res = 150)
term_gene_graph_manual(result_df = pathfindR_pathway_df_FC1[with(pathfindR_pathway_df_FC1, order(lowest_p,-Fold_Enrichment)),][][1:10,], use_description = TRUE , node_size = 'p_val')
dev.off()


tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_term_UpSet_plot.tiff"),width=20, height=20, units = 'in', res = 150)
UpSet_plot(result_df = pathfindR_pathway_df_FC1, genes_df = input2, method = 'heatmap')
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_GO_BP_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_BP_df_FC1 <- run_pathfindR(
  gene_sets = 'GO-BP',
  input3,
  output_dir = paste0('./GO/BP_FC1_DOWN/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_GO_CC_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_CC_df_FC1 <- run_pathfindR(
  gene_sets = 'GO-CC',
  input3,
  output_dir = paste0('./GO/CC_FC1_DOWN/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_GO_MF_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_MF_df_FC1 <- run_pathfindR(
  gene_sets = 'GO-MF',
  input3,
  output_dir = paste0('./GO/MF_FC1_DOWN/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_CM_Top10_pathway.tiff"),width=15, height=7, units = 'in', res = 150)
pathfindR_CM_df_FC1 <- run_pathfindR(
  gene_sets = 'cell_markers',
  input3,
  output_dir = paste0('./cell_markers_FC1_DOWN/',Condition),
  adj_method = 'fdr',
  enrichment_threshold = 0.05)
dev.off()

dev.off()

clustered_BP_FC1 <- cluster_enriched_terms(pathfindR_BP_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_BP_top5_FC1 <- subset(clustered_BP_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_Top5_cluster_BP.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_BP_top5_FC1, plot_by_cluster = TRUE)
dev.off()

clustered_CC_FC1 <- cluster_enriched_terms(pathfindR_CC_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_CC_top5_FC1 <- subset(clustered_CC_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_Top5_cluster_CC.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_CC_top5_FC1, plot_by_cluster = TRUE)
dev.off()

clustered_MF_FC1 <- cluster_enriched_terms(pathfindR_MF_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_MF_top5_FC1 <- subset(clustered_MF_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_Top5_cluster_MF.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_MF_top5_FC1, plot_by_cluster = TRUE)
dev.off()

clustered_CM_FC1 <- cluster_enriched_terms(pathfindR_CM_df_FC1, plot_clusters_graph = F, plot_dend = F)
clustered_CM_top5_FC1 <- subset(clustered_CM_FC1, Cluster %in% 1:5)
tiff(file=paste0('./plot/',Condition,'/',Condition,"_FC1_DOWN_Top5_cluster_CM.tiff"),width=15, height=7, units = 'in', res = 150)
enrichment_chart(clustered_CM_top5_FC1, plot_by_cluster = TRUE)
dev.off()

dev.off()

