rm(list = ls())
library(dplyr)
library(tibble)
library(data.table)
library(openxlsx)
library(limma)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggsci)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(VennDiagram)
library(grid)


setwd('Huaxi_Pancreatitis')
data_exp <-  read.xlsx('raw_data/data all gene.xlsx',rowNames = T) %>% as.data.frame()
data_infomation <- read.xlsx('raw_data/data info_new.xlsx',rowNames = F)
data_ap <- data_exp[order(data_infomation$stat)]
data_info <- as.factor(data_infomation$stat)
data_info1 <- data_infomation[order(data_infomation$stat),]
design <- model.matrix(~0 + data_info) 
colnames(design) <- levels(data_info)
rownames(design) <- colnames(data_exp)

# limma
contr.matrix <- makeContrasts(
  NonsevereVsHealth = Nonsevere - Health,
  SevereVsHealth = Severe - Health,
  SevereVsNonsevere = Severe - Nonsevere,
  levels  = colnames(design)
)
fit <- lmFit(data_exp, design)
fit <- contrasts.fit(fit, contrasts = contr.matrix)
efit <- eBayes(fit, trend = T)
plotSA(efit)
summary(decideTests(efit))
NonsevereVsHealth <- topTable(efit, coef = 1, n=Inf)
NonsevereVsHealth = NonsevereVsHealth[(NonsevereVsHealth$adj.P.Val)<0.05 & abs(NonsevereVsHealth$logFC) > 1.5,]

SevereVsHealth <- topTable(efit, coef = 2, n=Inf)
severe.volcano <- SevereVsHealth
SevereVsHealth = SevereVsHealth[(SevereVsHealth$adj.P.Val)<0.05 & abs(SevereVsHealth$logFC) > 1.5,]

SevereVsNonsevere <- topTable(efit, coef = 3, n=Inf)
SevereVsNonsevere = SevereVsNonsevere[(SevereVsNonsevere$adj.P.Val)<0.05 & abs(SevereVsNonsevere$logFC) > 1.5,]

logfcfilterprotein <- Reduce(intersect, list(
  rownames(NonsevereVsHealth),
  rownames(SevereVsHealth),
  rownames(SevereVsNonsevere)
))
logfc.sum <- data.frame(NonsevereVsHealth[logfcfilterprotein,],SevereVsHealth[logfcfilterprotein,],SevereVsNonsevere[logfcfilterprotein,])
logfc.sum.logfc.p <- dplyr::select(logfc.sum,c('logFC','adj.P.Val','logFC.1','adj.P.Val.1','logFC.2','adj.P.Val.2'))

#-----volcano visual-----
# volcano visual
severe.volcano <- NonsevereVsHealth
severe.volcano <- SevereVsHealth
severe.volcano$changed <- factor(ifelse(severe.volcano$adj.P.Val < 0.01 & abs(severe.volcano$logFC) > 1, ifelse(severe.volcano$logFC > 1, 'Up','Down'),'Not sig.'))
severe.volcano$selectedgene <- ifelse(severe.volcano$adj.P.Val < 0.001 & abs(severe.volcano$logFC) > 1.5, rownames(severe.volcano),NA)
summary(severe.volcano$changed)
ggplot(severe.volcano, aes(logFC, -log10(adj.P.Val),
                           color = factor(changed)
))+
  geom_point(size=3)+
  labs(x = expression(Log[2]*' Fold Change'),
       y = expression(-Log[10]*' (adj.P value)')) +
  theme_grey(base_size = 15) +
  scale_color_manual(values = c('deepskyblue', 'grey', 'tomato')) +
  scale_size_manual(values = c(2,1,2)) +
  theme(legend.title = element_blank(),
        legend.position = c(0.5,0.9),
        legend.background = element_rect(fill='transparent'))+
  geom_text_repel(aes(label = selectedgene), color = 'black', size = 4,
                  #box.padding = unit(0.5, 'lines'),
                  point.padding = NA,
                  segment.color = 'black') +
  geom_hline(yintercept = -log10(0.01),linetype=2,cex=1)+ 
  geom_vline(xintercept = c(-1,1),linetype=2,cex=1)+

  theme(axis.text= element_text(colour = "black"),
        panel.border = element_rect(size=1,fill='transparent'))
 #-----------venn plot
genes_NonsevereVsHealth <- rownames(NonsevereVsHealth)
genes_SevereVsHealth <- rownames(SevereVsHealth)
my_colors <- c("#EE2C2C", "#DAA520")

venn.plot <- venn.diagram(
  x = list(
    "Nonsevere vs Health" = genes_NonsevereVsHealth,
    "Severe vs Health"    = genes_SevereVsHealth
  ),
  filename = NULL,         
  fill = my_colors,        
  alpha = c(0.6, 0.6),    
  cex = 2,                 
  fontface = "bold",       
  fontfamily = "sans",     
  cat.col = my_colors,     
  cat.cex = 1.5,        
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.pos = c(-20, 20),    
  cat.dist = c(0.05, 0.05),
  margin = 0.1,         
  main = "Differential Genes Overlap",  
  main.cex = 2                    
)

grid.newpage()
grid.draw(venn.plot)

#-----pheatmap---------
deg_genes <- union(rownames(NonsevereVsHealth), rownames(SevereVsHealth))

data_heatmap <- data_ap[rownames(data_ap) %in% deg_genes, ]
annotation_col <- data.frame(SampleType = factor(data_info1$stat))
rownames(annotation_col) <- colnames(data_heatmap)
an_colors <- list(SampleType = c(Health = '#7CA878', Nonsevere = '#F4A2A3', Severe = '#EE7072'))

pheatmap(data_heatmap,
         cluster_cols = FALSE,
         scale = 'row',
         border_color = NA,
         show_rownames = TRUE,
         show_colnames = FALSE,
         cutree_rows = 5,
         annotation_col = annotation_col,
         annotation_colors = an_colors,
         color = colorRampPalette(c('#4552A0', "white", "firebrick3"))(50))