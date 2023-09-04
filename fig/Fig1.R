library(patchwork)
library(ggbeeswarm)
library(scales)
library(effsize)
#===============================================#
#                 limited model                 #
#===============================================#

pvlTable <- read.table("/home/yincheng23/ADlasso_manuscript/Script/pvl_distributed.txt", sep = "\t", stringsAsFactors = F, header = T)
pvlTable <- pvlTable[pvlTable$Method %in% c('ADlasso','LASSO_limited','SVMLASSO_limited','EN_limited','RF_limited','XGBoost_limited','MI'), ]
for(i in 1:nrow(pvlTable)){
  pvlTable$Method[i] <- gsub('_limited','',pvlTable$Method[i])
}
pvlTable$Method <- factor(pvlTable$Method, levels = c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI"))
cl <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628")
ll <- c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI")

la <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'LASSO'])
sl <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'SVMLASSO'])
en <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'EN'])
rf <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'RF'])
xg <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'XGBoost'])
mi <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'MI'])

text <- data.frame(Method = c("LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI"),
                   pos = c(1.5, 1.5, 1, 6, 5.7, 6.8),
                   cohenD = c(la$estimate,sl$estimate,en$estimate,rf$estimate,xg$estimate,mi$estimate))
text$Method <- factor(text$Method, c("LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI"))
text$cohenD <- round(text$cohenD,2)
p1 <- ggplot(pvlTable, aes(x = Method, y = pvl, color = Method)) + geom_boxplot(outlier.shape = NA) + ylim(c(0,8)) +
  scale_color_manual(values = cl) + ylab("Prevalence (%)") + ggtitle('Equivalent size model') + 
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(colour='black', size=12),
        axis.title.x = element_blank(),
        axis.text = element_text(colour='black',size=10),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = 'none') + 
        geom_text(data = text, aes(x=Method, y=pos, label = cohenD), position = position_dodge(width=0.9), fontface = "bold")


PCdf <- read.table("/home/yincheng23/ADlasso_manuscript/Script/performance_LIBSVM.txt", sep = "\t", stringsAsFactors = F, header = T)
PCdf <- PCdf[PCdf$Matrix == "AUC",]
PCdf <- PCdf[PCdf$Method %in% c('ADlasso','LASSO_limited','SVMLASSO_limited','EN_limited','RF_limited','XGBoost_limited','MI'), ]
for(i in 1:nrow(PCdf)){
  PCdf$Method[i] <- gsub('_limited','',PCdf$Method[i])
}
PCdf$Method <- factor(PCdf$Method, levels = c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI"))
cl <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628")
ll <- c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI")

p2 <- ggplot(PCdf, aes(x = Method, y = value, fill = Method)) +
  scale_fill_manual(values = cl) + scale_x_discrete(labels= ll) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(colour='black',size=10),
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none')
p2 <- p2 + ylab("AUC") + geom_bar(stat = "identity")

p1 / p2



#=================================================#
#                 unlimited model                 #
#=================================================#

pvlTable <- read.table("/home/yincheng23/ADlasso_manuscript/Script/pvl_distributed.txt", sep = "\t", stringsAsFactors = F, header = T)
pvlTable_ <- pvlTable[pvlTable$Method %in% c('ADlasso','LASSO','SVMLASSO','EN','RF','XGBoost'), ]
pvlTable_$Method <- factor(pvlTable_$Method, levels = c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost"))
cl <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33")
ll <- c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost")

la <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'LASSO'])
sl <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'SVMLASSO'])
en <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'EN'])
rf <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'RF'])
xg <- cohen.d(pvlTable$pvl[pvlTable$Method == 'ADlasso'], pvlTable$pvl[pvlTable$Method == 'XGBoost'])

text <- data.frame(Method = c("LASSO", "SVMLASSO", "EN", "RF", "XGBoost"),
                   pos = c(1.5, 1.5, 1, 4, 6),
                   cohenD = c(la$estimate,sl$estimate,en$estimate,rf$estimate,xg$estimate))
text$Method <- factor(text$Method, c("LASSO", "SVMLASSO", "EN", "RF", "XGBoost"))
text$cohenD <- round(text$cohenD,2)

p3 <- ggplot(pvlTable_, aes(x = Method, y = pvl, color = Method)) + geom_boxplot(outlier.shape = NA) + ylim(c(0,7.5)) +
  scale_color_manual(values = cl) + ylab("Prevalence (%)") + ggtitle('Full feature size model') +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(colour='black', size=12),
        axis.title.x = element_blank(),
        axis.text = element_text(colour='black',size=10),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 18),
        legend.position = 'none') +
  geom_text(data = text, aes(x=Method, y=pos, label = cohenD), position = position_dodge(width=0.9), fontface = "bold")

PCdf <- read.table("/home/yincheng23/ADlasso_manuscript/Script/performance_LIBSVM.txt", sep = "\t", stringsAsFactors = F, header = T)
PCdf <- PCdf[PCdf$Matrix == "AUC",]
PCdf <- PCdf[PCdf$Method %in% c('ADlasso','LASSO','SVMLASSO','EN','RF','XGBoost'), ]
PCdf$Method <- factor(PCdf$Method, levels = c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost"))
cl <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628")
ll <- c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI")

p4 <- ggplot(PCdf, aes(x = Method, y = value, fill = Method)) +
  scale_fill_manual(values = cl) + scale_x_discrete(labels= ll) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(colour='black',size=10),
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = 'none')
p4 <- p4 + ylab("AUC") + geom_bar(stat = "identity")


numtable <- data.frame(table(pvlTable_$Method))
colnames(numtable) <- c('Method', 'number')
numtable$Method <- factor(numtable$Method, levels = c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost"))
cl <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628")
ll <- c("ADlasso", "LASSO", "SVMLASSO", "EN", "RF", "XGBoost", "MI")

p5 <- ggplot(numtable, aes(x = Method, y = number, fill = Method)) +
  scale_fill_manual(values = cl) + scale_x_discrete(labels= ll) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_blank(),
        axis.text = element_text(colour='black',size=10),
        legend.title=element_blank(),
        axis.text.x = element_blank(),
        legend.position = 'none')
p5 <- p5 + ylab("Number") + geom_bar(stat = "identity") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))

pp1 <- p1 / p2 
pp2 <- p3 / p5 / p4
P <- pp1 | pp2


pdf("/home/yincheng23/ADlasso_manuscript/Figure/Figure1.pdf", width=8, height=8)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/Figure1.tiff", units="in", width=8, height=8, res=300)
P
dev.off()