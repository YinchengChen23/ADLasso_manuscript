setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/crc_zeller')
#KOmap <- read.table('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/ycc_crc1/mapedPW.txt', stringsAsFactors = F, sep = '\t', header = T)
#KO <- read.table('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/ycc_crc1/GSEA_w_normal.txt', stringsAsFactors = F, sep = '\t', header = T)
#KO <- read.table('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/ycc_crc1/GSEA_w_cancer.txt', stringsAsFactors = F, sep = '\t', header = T)
#KO <- read.table('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/ycc_crc1/GSEA_s_normal.txt', stringsAsFactors = F, sep = '\t', header = T)
#KO <- read.table('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/ycc_crc1/GSEA_s_cancer.txt', stringsAsFactors = F, sep = '\t', header = T)

KOmap <- read.table('/home/yincheng23/ADlasso_manuscript/Script/GSEA/submapedPW.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- read.table('/home/yincheng23/ADlasso_manuscript/Script/GSEA/GSEA_subs_cancer_II.txt', stringsAsFactors = F, sep = '\t', header = T)

KO <- KO[KO$z_score > 2,]
KO <- KO[KO$p_value < 0.05, ]
KO <- KO[KO$KO %in% KOmap$KO,]

enrichPW <- data.frame()
pathways <- unique(KOmap$pathway[KOmap$KO %in% KO$KO])
for(i in 1:length(pathways)){
  tag_p <- pathways[i]
  tag_pID <- unique(KOmap$mapid[KOmap$pathway == tag_p])
  KOinpathway <- unique(KOmap$KO[KOmap$pathway == tag_p])
  KOinpathway[KOinpathway %in% KO$KO]
  KOinpathway[!KOinpathway %in% KO$KO]
  KOnotinpathway <- unique(KOmap$KO[KOmap$pathway != tag_p])
  KOnotinpathway[KOnotinpathway %in% KO$KO]
  KOnotinpathway[!KOnotinpathway %in% KO$KO]
  r1 <- length(KOinpathway[KOinpathway %in% KO$KO])
  r2 <- length(KOinpathway[!KOinpathway %in% KO$KO])
  l1 <- length(KOnotinpathway[KOnotinpathway %in% KO$KO])
  l2 <- length(KOnotinpathway[!KOnotinpathway %in% KO$KO])
  m <- matrix(c(r1, r2, l1, l2), nrow = 2)
  res <- fisher.test(m)
  tmp <- data.frame('pathway' = tag_p, 'id' = tag_pID,'count' = r1, 'p' = res$p.value, 'oddratio' = as.numeric(res$estimate))
  enrichPW <- rbind(enrichPW, tmp)
}
enrichPW$adj.p <- p.adjust(enrichPW$p , 'fdr')
enrichPWE <- enrichPW[enrichPW$adj.p  < 0.05, ]
enrichPWE$pathway <- factor(enrichPWE$pathway, levels = enrichPWE$pathway[order(enrichPWE$count)])


KOmap <- read.table('/home/yincheng23/ADlasso_manuscript/Script/GSEA/submapedPW.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- read.table('/home/yincheng23/ADlasso_manuscript/Script/GSEA/GSEA_subs_normal_II.txt', stringsAsFactors = F, sep = '\t', header = T)

KO <- KO[KO$z_score > 2,]
KO <- KO[KO$p_value < 0.05, ]
KO <- KO[KO$KO %in% KOmap$KO,]

enrichPW <- data.frame()
pathways <- unique(KOmap$pathway[KOmap$KO %in% KO$KO])
for(i in 1:length(pathways)){
  tag_p <- pathways[i]
  tag_pID <- unique(KOmap$mapid[KOmap$pathway == tag_p])
  KOinpathway <- unique(KOmap$KO[KOmap$pathway == tag_p])
  KOinpathway[KOinpathway %in% KO$KO]
  KOinpathway[!KOinpathway %in% KO$KO]
  KOnotinpathway <- unique(KOmap$KO[KOmap$pathway != tag_p])
  KOnotinpathway[KOnotinpathway %in% KO$KO]
  KOnotinpathway[!KOnotinpathway %in% KO$KO]
  r1 <- length(KOinpathway[KOinpathway %in% KO$KO])
  r2 <- length(KOinpathway[!KOinpathway %in% KO$KO])
  l1 <- length(KOnotinpathway[KOnotinpathway %in% KO$KO])
  l2 <- length(KOnotinpathway[!KOnotinpathway %in% KO$KO])
  m <- matrix(c(r1, r2, l1, l2), nrow = 2)
  res <- fisher.test(m)
  tmp <- data.frame('pathway' = tag_p, 'id' = tag_pID,'count' = r1, 'p' = res$p.value, 'oddratio' = as.numeric(res$estimate))
  enrichPW <- rbind(enrichPW, tmp)
}
enrichPW$adj.p <- p.adjust(enrichPW$p , 'fdr')
#enrichPW <- enrichPW[enrichPW$p  < 0.05, ]
enrichPWS <- enrichPW[enrichPW$adj.p  < 0.05, ]
enrichPWS$pathway <- factor(enrichPWS$pathway, levels = enrichPW$pathway[order(enrichPW$count)])


p1 <- ggplot(enrichPWE, aes(x = pathway, y = count, fill = p)) + ggtitle('Enhanced in colorectal cancer') +
  geom_bar(stat="identity",colour = "black",size = 0.1) +  coord_flip() + ylab("KO Hits") + labs(fill = "q-value") +
  scale_fill_continuous(low = "red", high = "blue", limits = range(c(enrichPWE$p, enrichPWS$p))) +
  theme(axis.text.y = element_text(color='black', size='14'),
        axis.line = element_line(linetype = 1,colour = 'black'),
        axis.title.y = element_blank(),
        panel.background = element_rect(I(0)),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)) + coord_flip()

p2 <- ggplot(enrichPWS, aes(x = pathway, y = count, fill = p)) + ggtitle('Suppressed in colorectal cancer') +
  geom_bar(stat="identity",colour = "black",size = 0.1) +  coord_flip() + ylab("KO Hits") + labs(fill = "q-value") +
  scale_fill_continuous(low = "red", high = "blue", limits = range(c(enrichPWE$p, enrichPWS$p))) +
  theme(axis.text.y = element_text(color='black', size='14'),
        axis.title.y = element_blank(),
        axis.line = element_line(linetype = 1,colour = 'black'),
        panel.background = element_rect(I(0)),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)) + coord_flip()

p <- p1 + p2 + plot_layout(guides = "collect")

pdf("/home/yincheng23/ADlasso_manuscript/Figure/Figure5_2.pdf", width=14, height=4)
p
dev.off()



tiff("/home/yincheng23/ADlasso_manuscript/Figure/Figure5_2.tiff", units="in", width=14, height=4, res=300)
p
dev.off()



#===============================================================================
library(KEGGREST)
KO <- read.table('/home/yincheng23/ADlasso_manuscript/Script/GSEA/GSEA_subs_cancer_II.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- KO[KO$z_score > 2,]
KO <- KO[KO$p_value < 0.05, ]
KO$name = ""
KO$pathway = ""
for(i in 1:nrow(KO)){
  id <- KO$KO[i]
  query <- tryCatch(keggGet(id),
                    error=function(e) e)
  if(inherits(query, "error")) next
  KO$name[i] <- query[[1]]$NAME
  KO$pathway[i] <- paste(as.character(query[[1]]$PATHWAY), collapse = '; ')
}
write.csv(KO, '/home/yincheng23/ADlasso_manuscript/Script/GSEA/KO_sig_cancer.txt', sep = '\t')




KO <- read.table('/home/yincheng23/ADlasso_manuscript/Script/GSEA/GSEA_subs_normal_II.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- KO[KO$z_score > 2,]
KO <- KO[KO$p_value < 0.05, ]
KO$name = ""
KO$pathway = ""
for(i in 1:nrow(KO)){
  id <- KO$KO[i]
  query <- tryCatch(keggGet(id),
                    error=function(e) e)
  if(inherits(query, "error")) next
  KO$name[i] <- query[[1]]$NAME
  KO$pathway[i] <- paste(as.character(query[[1]]$PATHWAY), collapse = '; ')
}
write.csv(KO, '/home/yincheng23/ADlasso_manuscript/Script/GSEA/KO_sig_normal.txt', sep = '\t')








































setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies/ycc_crc1')
Descrip <- read.csv('KO/pred_metagenome_unstrat_descrip.tsv', as.is = TRUE, stringsAsFactors = F, sep = '\t', header = T)
KOmap <- read.table('mapedPW_cancer.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- read.table('GSEAn_cancer.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- KO[KO$z_score > 2, ]
KO <- KO[KO$p_value < 0.05, ]
KO <- KO[KO$KO %in% KOmap$KO,]

enrichPW <- data.frame()
pathways <- unique(KOmap$pathway[KOmap$KO %in% KO$KO])
for(i in 1:length(pathways)){
  tag_p <- pathways[i]
  tag_pID <- unique(KOmap$mapid[KOmap$pathway == tag_p])
  KOinpathway <- unique(KOmap$KO[KOmap$pathway == tag_p])
  KOinpathway[KOinpathway %in% KO$KO]
  KOinpathway[!KOinpathway %in% KO$KO]
  KOnotinpathway <- unique(KOmap$KO[KOmap$pathway != tag_p])
  KOnotinpathway[KOnotinpathway %in% KO$KO]
  KOnotinpathway[!KOnotinpathway %in% KO$KO]
  r1 <- length(KOinpathway[KOinpathway %in% KO$KO])
  r2 <- length(KOinpathway[!KOinpathway %in% KO$KO])
  l1 <- length(KOnotinpathway[KOnotinpathway %in% KO$KO])
  l2 <- length(KOnotinpathway[!KOnotinpathway %in% KO$KO])
  m <- matrix(c(r1, r2, l1, l2), nrow = 2)
  res <- fisher.test(m)
  tmp <- data.frame('pathway' = tag_p, 'id' = tag_pID,'count' = r1, 'p' = res$p.value, 'oddratio' = as.numeric(res$estimate))
  enrichPW <- rbind(enrichPW, tmp)
}
enrichPW <- enrichPW[enrichPW$p < 0.05, ]
enrichPW$pathway <- factor(enrichPW$pathway, levels = enrichPW$pathway[order(enrichPW$count)])
p1 <- ggplot(enrichPW, aes(x = pathway, y = count, fill = p)) +
  geom_bar(stat="identity",colour = "black",size = 0.1) +  coord_flip() + ylab("KO Hits") + labs(fill = "p-value") +
  xlab("pathway") + scale_fill_gradient(low = "red", high = "blue") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(face = "italic"),
        axis.line = element_line(linetype = 1,colour = 'black'),
        panel.background = element_rect(I(0)),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)) + coord_flip()



tiff("/home/yincheng23/ADlasso_manuscript/Figure/Figure6.tiff", units="in", width=6, height=6, res=300)
p1
dev.off()

unique(enrichPW$pathway)
KOmap$mapid[KOmap$pathway == 'Oxidative phosphorylation']





Descrip <- read.csv('KO/pred_metagenome_unstrat_descrip.tsv', as.is = TRUE, stringsAsFactors = F, sep = '\t', header = T)
KOmap <- read.table('mapedPW_normal.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- read.table('GSEAn_normal.txt', stringsAsFactors = F, sep = '\t', header = T)
KO <- KO[KO$z_score > 2, ]
KO <- KO[KO$p_value < 0.05, ]
KO <- KO[KO$KO %in% KOmap$KO,]

enrichPW <- data.frame()
pathways <- unique(KOmap$pathway[KOmap$KO %in% KO$KO])
for(i in 1:length(pathways)){
  tag_p <- pathways[i]
  tag_pID <- unique(KOmap$mapid[KOmap$pathway == tag_p])
  KOinpathway <- unique(KOmap$KO[KOmap$pathway == tag_p])
  KOinpathway[KOinpathway %in% KO$KO]
  KOinpathway[!KOinpathway %in% KO$KO]
  KOnotinpathway <- unique(KOmap$KO[KOmap$pathway != tag_p])
  KOnotinpathway[KOnotinpathway %in% KO$KO]
  KOnotinpathway[!KOnotinpathway %in% KO$KO]
  r1 <- length(KOinpathway[KOinpathway %in% KO$KO])
  r2 <- length(KOinpathway[!KOinpathway %in% KO$KO])
  l1 <- length(KOnotinpathway[KOnotinpathway %in% KO$KO])
  l2 <- length(KOnotinpathway[!KOnotinpathway %in% KO$KO])
  m <- matrix(c(r1, r2, l1, l2), nrow = 2)
  res <- fisher.test(m)
  tmp <- data.frame('pathway' = tag_p, 'id' = tag_pID,'count' = r1, 'p' = res$p.value, 'oddratio' = as.numeric(res$estimate))
  enrichPW <- rbind(enrichPW, tmp)
}
enrichPW <- enrichPW[enrichPW$p < 0.05, ]
enrichPW$pathway <- factor(enrichPW$pathway, levels = enrichPW$pathway[order(enrichPW$count)])
ggplot(enrichPW, aes(x = pathway, y = count, fill = p)) +
  geom_bar(stat="identity",colour = "black",size = 0.1) +  coord_flip() + ylab("KO Hits") + labs(fill = "p-value") +
  xlab("pathway") + scale_fill_gradient(low = "red", high = "blue") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(face = "italic"),
        axis.line = element_line(linetype = 1,colour = 'black'),
        panel.background = element_rect(I(0)),
        panel.grid.major = element_line(colour = NA),
        panel.grid.minor = element_line(colour = NA)) + coord_flip()

