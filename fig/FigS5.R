pvlTable <-  data.frame()
refpath <- '/home/yincheng23/miRNA/TCGA/selected'
selectList <- '/home/yincheng23/miRNA/TCGA/ML_para'
setwd(selectList)
dataSets <- dir()
for(cancerType in dataSets){
  setwd(refpath)
  setwd(cancerType)
  selectedList_AD <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F)
  setwd(selectList)
  setwd(cancerType)
  pvllist <- read.table('prevalence_list.txt', sep = "\t", stringsAsFactors = F)
  pvlAD <- as.numeric(pvllist$V3[pvllist$V2 %in% selectedList_AD$V1])
  print(cancerType)
  tmp <- data.frame('cancerType'=cancerType,'type'='ALL','prevalence'=as.numeric(pvllist$V3))
  pvlTable <- rbind(pvlTable, tmp)
  tmp <- data.frame('cancerType'=cancerType,'type'='ADlasso','prevalence'=pvlAD)
  pvlTable <- rbind(pvlTable, tmp)
}

p <- ggplot(pvlTable, aes(x = prevalence, fill = type)) +
     geom_density(alpha=.7, color='black') +
     facet_wrap(vars(cancerType), scales = "free", ncol=4) +
     scale_fill_manual(values = c('#FFAA33','#AAAAAA')) + 
     labs(x='Prevalence (%)', y='Density', fill='Features') + theme_bw(14) +
     theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"))

pdf("/home/yincheng23/ADlasso_manuscript/Figure/FigureS5.pdf", width=10, height=9)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/FigureS5.tiff", units="in", width=10, height=9, res=300)
p
dev.off()


