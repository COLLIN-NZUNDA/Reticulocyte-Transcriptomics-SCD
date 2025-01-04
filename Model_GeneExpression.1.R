#!/usr/bin/env Rscript

#library("cqn") #optional for cqn normalisation
#BiocManager::install("tximportData")
library("BiocManager")
library("DESeq2")
#####  GENES COUNT DATA LOADINGS PROCESSS
library("tximport")
library("tximeta")
library("GO.db")
library("tximportData")
library("magrittr")

### GENES EXPREWSSION, GENES COUNTS PROCESSING
library("goseq")
library("piano")
library("apeglm")
library("edgeR")
library("DGEAR")
library("tidyverse")
library("canvasXpress")
library("clinfun")
library("GGally")
library("factoextra")

## RESULTS VISUALISATION
library("dplyr")
library("ggplot2")
library(tidyr)
library(reshape)
library("biomaRt")
library(pheatmap)

#### FUNCTION GENE EXPRESSION PROCESSING
library("WGCNA")

library(ggrepel)
library(DEGreport)
library(RColorBrewer)


# https://github.com/CebolaLab/RNA-seq
#Read in the files with the sample information
setwd("/home/echimusa/Documents/RNA_SEQ_SIANA/salmon/")
samples = read.table('/home/echimusa/Documents/RNA_SEQ_SIANA/SAMPLE.txt',header = TRUE)
names(samples)

gtf <- rtracklayer::import("/home/echimusa/Documents/RNA_SEQ_SIANA/gencode.v46.basic.annotation.gtf.gz")

head(gtf)

gtf_df<- as.data.frame(gtf)
#attach(gtf_df)
names(gtf_df)

tx2gene <- unique(data.frame(gtf[gtf$type=="transcript" | gtf$type =="mRNA"])[,c("transcript_id", "gene_id")])

head(tx2gene)

gene_name_map<- unique(data.frame(gtf[gtf$type =="gene"] )[,c("gene_id","gene_name", "gene_type")] ) 
head(gene_name_map)

#####Column 2 of samples, samples[,2], contains the paths to the quant.sf files

counts.imported = tximport(files = as.character(samples[,2]), type = 'salmon', tx2gene = tx2gene)
samples$Conditions <- as.factor(samples$Conditions)

#Import to DEseq2
counts.DEseq = DESeqDataSetFromTximport(counts.imported, colData = samples, design = ~Conditions)

dds <- DESeq(counts.DEseq)
resultsNames(dds) #lists the coefficients
dds
#plotDispEsts(dds)
#Add the normalisation offset from cqn
##normalizationFactors(dds) <- cqnNormFactors


#### Plotting Counts dispersion versus mean of mnormalised counts #######################################################################
slotNames(dds)
png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure1.png", res =600, width = 10080, height = 6580)
plotDispEsts(dds) #,ylim = c(1e-6, 1e1))
dev.off()
######################### End ! #########################################################################################################

######################## Write first result without normalisation steps #################################################################
res <- results(dds)
head(results(dds, tidy=TRUE)) 
res <- res[order(res$padj),]
res <- as.data.frame(res)
res <- as.data.frame(cbind(transcript_id=gene_name_map[match (gene_name_map[,1], rownames(res)),1], Gene_name=gene_name_map[match (gene_name_map[,1], rownames(res)),2], Gene_type=gene_name_map[ match(gene_name_map[,1], rownames(res)),3], res))
head(res, n=3)
RES.sig <- res[res$padj < 0.05 & !is.na(res$padj),] #subset the significant genes  7.9e-07
RES.sig <- RES.sig[order(RES.sig$log2FoldChange, decreasing = TRUE),]
RES.top <- LFC[ (RES.sig$baseMean > 50) & (abs(RES.sig$log2FoldChange) > 2),]
head(RES.sig,n=25)
write.table(res, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/Table1.txt', sep="\t", row.names = F)
#######################    End!     #####################################################################################################
#######################   Normalisation #################################################################################################
nm <- assays(dds)[["avgTxLength"]]

sf1 <- estimateSizeFactorsForMatrix(counts(dds))
sf2  <- estimateSizeFactorsForMatrix(normalizationFactors(dds)) 
sf3 <- estimateSizeFactorsForMatrix(counts(dds) / nm)

normalized_counts <- counts(dds, normalized=TRUE)
row.names(normalized_counts)

GeneList <- as.data.frame(cbind(gene_id=unlist(lapply(strsplit(gene_name_map[match (gene_name_map[,1], row.names(normalized_counts)),1], '\\.'), '[[',1)), gene_name=gene_name_map[match (gene_name_map[,1], rownames(normalized_counts)),2], transcript_id=gene_name_map[ match(gene_name_map[,1], rownames(normalized_counts)),1])) #, as.data.frame( normalized_counts)))

write.table(GeneList, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/geneList.txt',row.names = FALSE)

####################################  WRITING NORMALISED COUNTS FOR C0_EXPRESSION NETWORKS AND MODULES ##################################
t <- as.vector(samples[,1]) ; t1 <- c("sample_id") ; s <-c(t1, t) 
df_norm <-as.data.frame(t(normalized_counts))
df_norm <- as.data.frame(cbind(sample_id=samples[,1], df_norm))
write.csv(df_norm, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/expression_data.csv',row.names = FALSE)
############################################################ END #########################################################################

###plotMA(normalized_counts, main = '???', cex = 8.5, ylim=c(-2,100000), xlim=c(0,100000))

counts <- as.data.frame(cbind(gene_id=unlist(lapply(strsplit(gene_name_map[match (gene_name_map[,1], row.names(normalized_counts)),1], '\\.'), '[[',1)), gene_name=gene_name_map[match (gene_name_map[,1], rownames(normalized_counts)),2], Gene_type=gene_name_map[ match(gene_name_map[,1], rownames(normalized_counts)),3], as.data.frame(normalized_counts)))

write.table(counts, file="/home/echimusa/Documents/RNA_SEQ_SIANA/normalized_counts.txt", sep="\t", row.names = F)

########################################################### COUNTS AGAINST PHENOTYPES AND OTHERS SAMPLES CHARACTERISTICS #####################

LFC <- lfcShrink(dds, coef = "Conditions_WholeBlood_vs_CordBlood", type = "apeglm")
LFC <- LFC[order(LFC$padj),]
LFC <- LFC[order(LFC$log2FoldChange, decreasing = TRUE),]
LFC1 <- as.data.frame(cbind(transcript_id=gene_name_map[match (gene_name_map[,1], rownames(LFC)),1], Gene_name=gene_name_map[match (gene_name_map[,1], rownames(LFC)),2], Gene_type=gene_name_map[ match(gene_name_map[,1], rownames(LFC)),3], as.data.frame(LFC)))
write.table(LFC1, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/Table.2.LFC_Gene_Epression.txt', sep="\t", row.names = F)

############################################### TOP SIGNIFICANT GENES ########################################################################
LFC.sig <- LFC[LFC$padj < 0.05 & !is.na(LFC$padj),] #subset the significant genes  7.9e-07
LFC.sig <- LFC.sig[order(LFC.sig$log2FoldChange, decreasing = TRUE),]
LFC.sig <- as.data.frame(LFC.sig)
LFC1 <- as.data.frame(cbind( gene_name_map %>% filter(gene_name_map$gene_id %in% row.names(LFC.sig)), LFC.sig))

####################### COUNTS OF TO 631 GENES AND GENE LIST ###################
GeneList_top631 <- as.data.frame(cbind(gene_id= unlist(lapply(strsplit(LFC1$gene_id, '\\.'), '[[',1)), gene_name=LFC1$gene_name, transcript_id=LFC1$gene_id))
LFC1$gene_id <- unlist(lapply(strsplit(LFC1$gene_id, '\\.'), '[[',1))
write.table(LFC1, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/Table.3.LFC_ToP631.Gene_Epression.txt', sep="\t", row.names = F)
write.table(GeneList_top631, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/geneList_top631.txt',row.names = FALSE)


####################################  WRITING NORMALISED COUNTS FOR C0_EXPRESSION NETWORKS AND MODULES #####################################
t <- as.vector(samples[,1]) ; t1 <- c("sample_id") ; s <-c(t1, t) 
normalized_counts <-  as.data.frame( normalized_counts)
top631_sigOE_norm <- normalized_counts  %>% filter(rownames(normalized_counts) %in% GeneList_top631$transcript_id) 
df_norm_631 <-as.data.frame(t(top631_sigOE_norm ))
df_norm_631 <- as.data.frame(cbind(sample_id=samples[,1], df_norm_631))
write.csv(df_norm_631, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/expression_data_631.csv',row.names = FALSE)

############################################################    END     ######################################################################

df.top <- LFC.sig[ (LFC.sig$baseMean > 50) & (abs(LFC.sig$log2FoldChange) > 2),]
df.top <- as.data.frame(df.top[order(df.top$log2FoldChange, decreasing = TRUE),])
LFC2 <- as.data.frame(cbind( gene_name_map %>% filter(gene_name_map$gene_id %in% row.names(df.top)), df.top))
LFC2 <- as.data.frame(LFC2[order(LFC2$log2FoldChange, decreasing = TRUE),])
head(LFC2, n=25)
GeneList_top47 <- as.data.frame(cbind(gene_id=unlist(lapply(strsplit(LFC2$gene_id, '\\.'), '[[',1)), gene_name=LFC2$gene_name, transcript_id=LFC2$gene_id))
LFC2$gene_id <- unlist(lapply(strsplit(LFC2$gene_id, '\\.'), '[[',1))
write.table(LFC2, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/Table.4.LFC_ToP47.Gene_Epression.txt', sep="\t", row.names = F)
write.table(GeneList_top47, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/geneList_top47.txt',row.names = FALSE)

top47_sigOE_norm <- normalized_counts  %>% filter(rownames(normalized_counts) %in% GeneList_top47$transcript_id) 
df_norm_47 <-as.data.frame(t(top47_sigOE_norm ))
df_norm_47 <- as.data.frame(cbind(sample_id=samples[,1], df_norm_47))
write.csv(df_norm_47, file='/home/echimusa/Documents/RNA_SEQ_SIANA/Tables/expression_data_47.csv',row.names = FALSE)

#vsd <- varianceStabilizingTransformation(dds)
#plotMA(LFC)

################### HOW MANY GENES ARE DIFFERENTIALLY EXPRESSED ######################
#View the top 10 genes with the most significant (adjusted) p-values
head(LFC.sig, n = 10)
attach(as.data.frame(LFC))
#The total number of DEGs with an adjusted p-value<0.05
summary(LFC, alpha=0.05)
#The total number of DEGs with an adjusted p-value<0.05 AND absolute fold-change > 2
sum(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) >2)
#Decreased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange <(-2)) #fold-change greater than 2
#Increased expression:
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >0) #any fold-change
sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) #fold-change greater than 2

##########################I TOP GENES DIFFERENTIALLY EXPRESSED  FOR R WGCNA ################################

top631_sigOE_norm2 <- top631_sigOE_norm

############################################################## PREPERE HEATMAP INPUTS ############################
## normalized counts for top 20 significant genes

top47_sigOE_norm2 <- top47_sigOE_norm
top25_sigOE_genes <- head(rownames(LFC2), n=25)
names(LFC2) 

top25_sigOE_norm <- top47_sigOE_norm  %>% filter(rownames(top47_sigOE_norm) %in% top25_sigOE_genes)   #%>%  add_row(samples$Conditions,.before = 2)
names(top25_sigOE_norm) <- samples[,1]
top25_sigOE_norm$Gene <- LFC2[match(rownames(top25_sigOE_norm), rownames(LFC2)),2]
rownames(top25_sigOE_norm) <- unlist(lapply(strsplit(rownames(top25_sigOE_norm) , '\\.'), '[[',1))
top25_sigOE_norm2 <- top25_sigOE_norm 

top25_sigOE_norm <- data.frame(melt(top25_sigOE_norm))

colnames(top25_sigOE_norm) <- c("gene", "SAMPLEID", "normalized_count")
Sam <- samples[,-2]
top25_sigOE_norm <- inner_join(Sam, top25_sigOE_norm , by="SAMPLEID")


melted_top47_sigOE <-  as.data.frame(top47_sigOE_norm)
names(melted_top47_sigOE ) <- samples[,1]
melted_top47_sigOE$Gene <- GeneList_top47[match(rownames(melted_top47_sigOE), GeneList_top47$transcript_id),2]
rownames(melted_top47_sigOE) <- unlist(lapply(strsplit(rownames(melted_top47_sigOE) , '\\.'), '[[',1))
melted_top47_sigOE <- data.frame(melt(melted_top47_sigOE))
colnames(melted_top47_sigOE) <- c("gene", "SAMPLEID", "normalized_count")
melted_top47_sigOE <- inner_join(Sam, melted_top47_sigOE , by="SAMPLEID")

head(LFC2, n=25)

########################################################## PREPARE WGCNA ############################################


################################################## END ##############################################################


###################  PVALUUES DITRIBUTIONS CHECK  ##############################
png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure2.png", width = 780, height = 480)
par(mfrow=c(2,1))
hist(LFC$pvalue, breaks = 50, col = 'green' , main = 'Distribution of Pvalue Before adjusted', xlab = 'p-value')
hist(LFC$padj, breaks = 50, col = 'blue', main = 'Distribution of Pvalue After adjusted', xlab = 'Adjusted p-value')

dev.off()
 

### ###########################  To generate the PCA plot with any batch effects removed:
png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure3.png", width = 780, height = 480)
vst.r <- vst(dds,blind = TRUE)
vst_mat <- assay(vst.r)
pca <- prcomp(t(vst_mat))
#Batch effect (donor) removed
assay(vst.r) <- limma::removeBatchEffect(assay(vst.r), vst.r$batch)

z=plotPCA(vst.r, "Conditions")
nudge <- position_nudge(y = 1,x=4)
z + geom_text(size=4.5, aes(label = name), position = nudge) + theme()

#An example with no labels
#z=plotPCA(vst.r, "Conditions")
#nudge <- position_nudge(y = 1,x=4)
#z + geom_text(size=2.5, aes(label = NA), position = nudge) + theme()
dev.off()
###################### ALL VOLCANO PLOTS ########################################

#Allow for more space around the borders of the plot
png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure4.png", width = 780, height = 480)
par(mar = c(5, 4, 4, 4))
#Set your log-fold-change and p-value thresholds
lfc = 2
pval = 7.9e-07  #0.05
tab = data.frame(logFC = LFC$log2FoldChange, negLogPval = -log10(LFC$padj))#make a data frame with the log2 fold-changes and adjusted p-values
plot(tab, pch = 16, cex = 1.6, xlab = expression(log[2]~fold~change),
     ylab = expression(-log[10]~pvalue), main = 'Volcano') #replace main = with your title

#Genes with a fold-change greater than 2 and p-value<0.05:
signGenes = (abs(tab$logFC) > lfc & tab$negLogPval > -log10(pval))

#Colour these red
points(tab[signGenes, ], pch = 16, cex = 1.9, col = "red")
#Show the cut-off lines
abline(h = -log10(pval), col = "green3", lty = 4)
abline(v = c(-lfc, lfc), col = "blue", lty = 4)

mtext(paste("FDR =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.7, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.8, line = 0.7)

dev.off()


###################  Expression for the top genes ###################################
#Select your chosen gene 
png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure5.png", width = 780, height = 458)
tmp = plotCounts(dds, gene = grep('ENSG00000291072', names(dds), value = TRUE), intgroup = "Conditions", pch = 18, main = 'ENSG00000291072 expression in lncRNA', returnData = TRUE)

theme <- theme(panel.background = element_blank(), panel.border = element_rect(fill = NA),
               plot.title = element_text(hjust = 0.5))

p <- ggplot(tmp, aes(x = Conditions, y = count, fill=Conditions)) + geom_boxplot() + geom_boxplot() + 
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.6) + ggtitle('ENSG00000291072 expression in lncRNA') + theme

dev.off() 
#print(p)


############################  PlOTTING TOP 20 SIGNIFICANT GENES ########################################
png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure6.png", width = 880, height = 380)
ggplot(top25_sigOE_norm) +
  geom_point(aes(x = gene, y = normalized_count, color = Conditions,size = 35*normalized_count)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("Normalized Counts") +
  ggtitle("Top 25 Significant DE Genes") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title=element_text(hjust=0.5))
dev.off()


##################################  SAMPLE HEATMAP PLOTTING #####################
#norm_OEsig <- normalized_counts[rownames(sigOE),]
################# Optional #######################
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
converted3 <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id',
                    values = rownames(top20_sigOE_norm), mart = ensembl)

top20_sigOE_norm$Gene <- converted3[match (converted3[,2], rownames(top20_sigOE_norm)),1]
nrow(top20_sigOE_norm)
names(top20_sigOtop20_sigOE_normE_norm)
top20_sigOE_norm <-  as.data.frame(top20_sigOE_norm)

rownames(top20_sigOE_norm) <- converted3[match (converted3[,2], rownames(top20_sigOE_norm)),1]

names(top25_sigOE_norm)


### Annotate our heatmap (optional)
annotation <- data.frame(sampletype=as.character(samples[,3]),row.names=colnames(top25_sigOE_norm))
col_annotation = rownames(top25_sigOE_norm)
row_annotation <- colnames(top25_sigOE_norm)


################################################  HeatMAP Plot ########################
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
### Set a color palette
heat.colors <- brewer.pal(6, "YlOrRd") 

top25_sigOE_norm2 <- as.matrix(top25_sigOE_norm2[,-25])

mat.scaled <- t(apply(top25_sigOE_norm2 , 1, scale)) #center and scale each column (Z-score) then transpose
colnames(mat.scaled)<-colnames(top25_sigOE_norm2)


num_keep <- 25
#1 to num_keep len-num_keep to len
rows_keep <- c(seq(1:num_keep), seq((nrow(mat.scaled)-num_keep), nrow(mat.scaled)) )

l2_val <- as.matrix(df.top[rows_keep,]$log2FoldChange) #getting log2 value for each gene we are keeping
colnames(l2_val)<-"logFC"

mean <- as.matrix(df.top[rows_keep,]$baseMean) #getting mean value for each gene we are keeping
colnames(mean)<-"AveExpr"

col_logFC <- colorRamp2(c(min(l2_val),0, max(l2_val)), c("blue", "white", "red")) 
#maps between 0% quantile, and 75% quantile of mean values --- 0, 25, 50, 75, 100
col_AveExpr <- colorRamp2(c(quantile(mean)[1], quantile(mean)[4]), c("white", "red"))


png("/home/echimusa/Documents/RNA_SEQ_SIANA/Figures/Figure7.png", res =600, width = 10080, height = 6580)

ha <- HeatmapAnnotation(summary = anno_summary(gp = gpar(fill = 6), 
                                               height = unit(6, "cm")))

h1 <- Heatmap(mat.scaled[rows_keep,], cluster_rows = F, 
              column_labels = colnames(mat.scaled), name="Z-score",
              cluster_columns = T)
h2 <- Heatmap(l2_val, row_labels = df.top$Gene[rows_keep], 
              cluster_rows = F, name="logFC", top_annotation = ha, col = col_logFC,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(l2_val[i, j],2), x, y)
              })
h3 <- Heatmap(mean, row_labels = df.top$Gene[rows_keep], 
              cluster_rows = F, name = "AveExpr", col=col_AveExpr,
              cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
                grid.text(round(mean[i, j],2), x, y)
              })

h<-h1+h2+h3
h
dev.off()



######################################### ENRICHMENT ANALYSIS ON TOP SIGNFICANT GENES ##################
nrow(as.data.frame(LFC.sig))
df.top
#The largest fold-changes with a significant p-value
library(tibble)
library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
#library(tidyverse)
library(fgsea)
# you may have to install some of these libraries; use 
# BiocManager::install(c("org.Hs.eg.db","clusterProfiler","enrichplot","fgsea"))
res_tableOE_tb <- df.top%>%
  data.frame() %>%
  rownames_to_column(var="gene") %>%
  dplyr::filter(!is.na(log2FoldChange))  %>% as_tibble()

###names(res_tableOE_tb)
##res_tableOE_tb$gene<- unlist(lapply(strsplit(res_tableOE_tb$gene, '\\.'), '[[',1))

LFC.sig$Gene


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

converted <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = res_tableOE_tb$gene, mart = ensembl)

#Add gene names to the LFC.sig data-frame
names(converted)
allOE_genes <-converted  %>%
  filter(converted[,2] %in% res_tableOE_tb$gene) %>% select(hgnc_symbol)
allOE_genes <- as.character(allOE_genes[allOE_genes !=""])
#res_tableOE_tb$gene <- converted[match (converted[,2], res_tableOE_tb$gene),1]
#sigOE <- dplyr::filter(res_tableOE_tb, padj < 0.00001)
sigOE <- res_tableOE_tb$Gene #  sigOE$gene


#sigOE <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id',values = sigOE, mart = ensembl)
#sigOE <- as.character(sigOE$hgnc_symbol[sigOE$hgnc_symbol!=""])
## Run GO enrichment analysis 
sigOE <- as.character(sigOE[sigOE !=""])
Uni <- as.character(LFC.sig$Gene[LFC.sig$Gene !=""])
ego <- enrichGO(gene = sigOE, 
                universe =Uni,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                minGSSize = 10,
                maxGSSize = 500,
                ont = "MF", 
                pAdjustMethod = "none", 
                pvalueCutoff  = 0.05,
                qvalueCutoff = 0.5)

#ont           = "CC",
#pAdjustMethod = "BH",
#pvalueCutoff  = 0.01,
#qvalueCutoff  = 0.05,
## Output results from GO analysis to a table
cluster_summary <- data.frame(ego)
## Dotplot 
dotplot(ego, showCategory = 50)
goplot(ego)
heatplot(ego)
upsetplot(ego)


## Enrichmap clusters the 50 most significant (by padj) GO terms to visualize relationships between terms

pwt <- pairwise_termsim(
  ego,
  method = "JC",
  semData = NULL,
  showCategory = 50
)

emapplot(pwt, showCategory = 50)


## To color genes by log2 fold changes, we need to extract the log2 fold changes from our results table creating a named vector

sigOE <- dplyr::filter(res_tableOE_tb, padj < 0.00001)

OE_foldchanges <- df.top$log2FoldChange

#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

#converted <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id', values = df.top$Ge, mart = ensembl)
#names(OE_foldchanges) <- converted$hgnc_symbol

## Cnetplot details the genes associated with one or more terms - by default gives the top 5 significant terms (by padj)
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 5, 
         foldChange = OE_foldchanges, vertex.label.font=6)
         
######################################### GSEA ##################
## Extract the foldchanges
foldchanges <- res_tableOE_tb$log2FoldChange


ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

converted <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = res_tableOE_tb$gene, mart = ensembl)
## Name each fold change with the corresponding Entrez ID
names(foldchanges) <- converted  %>%
  filter(converted[,2] %in% res_tableOE_tb$gene) %>% select(hgnc_symbol)
#subset the significant genes

foldchanges <- sort(foldchanges, decreasing = TRUE)

#

# GSEA using gene sets associated with BP Gene Ontology terms

gseaGO <- clusterProfiler::gseGO(
  geneList = as.dataframe(df.top),
  ont = "BP",
  keyType = "SYMBOL",
  eps = 0,
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  verbose = TRUE,
  OrgDb = "org.Hs.eg.db",
  by = "fgsea"
)

gseaGO_results <- gseaGO@result

goplot(gseaGO)

gseaplot2(gseaGO, geneSetID = 1:3)


LFC.sig = LFC.gene[padj < 0.05 & !is.na(padj),]#subset the significant genes  7.9e-07

#We can add a column with the HGNC gene names

##### READ FROM BIOMART HUMAN DATABASES s###############################################
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

converted <- getBM(attributes=c('hgnc_symbol','ensembl_gene_id'), filters = 'ensembl_gene_id',
                   values = rownames(LFC.sig), mart = ensembl)
#Decreased expression
ql53.DEGs.down <- LFC.sig[LFC.sig$log2FoldChange < 0,]   # (-2) fold-change greater then 2
nrow(ql53.DEGs.down)

#Increased expression
ql53.DEGs.up <- LFC.sig[LFC.sig$log2FoldChange > 0,]   # 
nrow(ql53.DEGs.up); names(ql53.DEGs.up); ql53.DEGs.up$Gene_id
pwf.dn <- nullp(ql53.DEGs.up, "hg19", "Gene_id")
go.results.dn <- goseq(pwf.dn, "hg19", "ensGene")

attach(as.data.frame(LFC))


sum(!is.na(padj) & padj < 0.05 & log2FoldChange >2) #fold-change greater than 2

#Decreased expression
ql53.DEGs.down <- groups12.table$FDR < 0.05 & groups12.table$logFC<0
names(ql53.DEGs.down) <- rownames(groups12.table)
pwf.dn <- nullp(ql53.DEGs.up, "hg19", "ensGene")
go.results.dn <- goseq(pwf.dn, "hg19", "ensGene")

#Increased expression
ql53.DEGs.up <- groups12.table$FDR < 0.05 & groups12.table$logFC>0
names(ql53.DEGs.up) <- rownames(groups12.table)
pwf.up <- nullp(ql53.DEGs.down, "hg19","ensGene")
go.results.up <- goseq(pwf.up, "hg19","ensGene")
###########################################################################################################


## If you have few samples:




########################################
#
#                 END 
#  
#########################################


library(hdWGCNA)

#########################################################################################################################
#Extract the differential expression data, with false discovery rate correction

#### WCGA  ##############################################################################################################
######################################  WCCCGA #########################################################################
input_mat = df_norm_47  # mat.scaled  # t(expr_normalized)
input_mat[1:5,1:10]           # Look at first 5 rows and 10 columns
allowWGCNAThreads() 

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2))

# Call the network topology analysis function
sft = pickSoftThreshold(
  input_mat,             # <= Input data
  #blockSize = 30,
  powerVector = powers,
  verbose = 5
)

par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence")
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red"
)
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity")
)

text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

picked_power = 16
temp_cor <- cor    

cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor     # Return cor function to original namespace

# Convert labels to colors for plotting
mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = TRUE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

module_df[1:5,]
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")


#We have written out a tab delimited file listing the genes and their modules.
#However, we need to figure out which modules are associated with each trait/treatment group.
#WGCNA will calcuate an Eigangene (hypothetical central gene) for each module, 
#so it easier to determine if modules are associated with different treatments.


# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
MEs0$treatment = row.names(MEs0)

# tidy & plot data
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1,1)) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Module-trait Relationships", y = "Modules", fill="corr")


#Looking at the heatmap, the green seems positively associated (red shading) with the L groups, turquoise module is positive for the L but negative (blue shading) for the L_L1 groups, and tan module is positive in the S groups but negative elsewhere. There are other patterns, think about what the patterns mean. For this tutorial we will focus on the green, turquoise, and tan module genes.

#Examine Expression Profiles
#We’ll pick out a few modules of interest, and plot their expression profiles
subexpr = expr_normalized[submod$gene_id,]
# pick out a few modules of interest here
modules_of_interest = c("green", "turquoise", "tan")

# Pull out list of genes in that module
submod = module_df %>%
  subset(colors %in% modules_of_interest)

row.names(module_df) = module_df$gene_id

# Get normalized expression for those genes
expr_normalized[1:5,1:10]


submod_df = data.frame(subexpr) %>%
  mutate(
    gene_id = row.names(.)
  ) %>%
  pivot_longer(-gene_id) %>%
  mutate(
    module = module_df[gene_id,]$colors
  )

submod_df %>% ggplot(., aes(x=name, y=value, group=gene_id)) +
  geom_line(aes(color = module),
            alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90)
  ) +
  facet_grid(rows = vars(module)) +
  labs(x = "treatment",
       y = "normalized expression")


###Generate and Export Networks
###The network file can be generated for Cytoscape or as an edge/vertices file.
genes_of_interest = module_df %>%
  subset(colors %in% modules_of_interest)

expr_of_interest = expr_normalized[genes_of_interest$gene_id,]
expr_of_interest[1:5,1:5]


# Only recalculate TOM for modules of interest (faster, altho there's some online discussion if this will be slightly off)
TOM = TOMsimilarityFromExpr(t(expr_of_interest),
                            power = picked_power)


# Add gene names to row and columns
row.names(TOM) = row.names(expr_of_interest)
colnames(TOM) = row.names(expr_of_interest)

edge_list = data.frame(TOM) %>%
  mutate(
    gene1 = row.names(.)
  ) %>%
  pivot_longer(-gene1) %>%
  dplyr::rename(gene2 = name, correlation = value) %>%
  unique() %>%
  subset(!(gene1==gene2)) %>%
  mutate(
    module1 = module_df[gene1,]$colors,
    module2 = module_df[gene2,]$colors
  )

head(edge_list)

# Export Network file to be read into Cytoscape, VisANT, etc
write_delim(edge_list,
            file = "edgelist.tsv",
            delim = "\t")



