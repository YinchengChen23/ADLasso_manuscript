library(effsize)
library(patchwork)
library(RColorBrewer)

getCD <- function(DataSet, Method, ADp, Metp, addDF){
  cohD <- cohen.d(ADp, Metp)
  tmp <- data.frame("Data" = DataSet, "Method" = Method, "cohensD" = cohD$estimate, "number" = length(Metp))
  addDF <- rbind(addDF, tmp)
  return(addDF)
}

#===========================================#
#             limitedModel                  #
#===========================================#

limitedModel = c("LASSO_limited.txt", "EN_limited.txt", "RF_limited.txt", "XGBoost_limited.txt", "MI.txt", "mRMR.txt", "ReliefF.txt", "fisher.txt", "FDC.txt")
pvlTable <-  data.frame()
abTable <-  data.frame()
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
  
  ablist <- read.table('meanAB_list.txt', sep = "\t", stringsAsFactors = F)
  abAD <- as.numeric(ablist$V3[ablist$V2 %in% selectedList_AD$V1])
  
  print(cancerType)
  for(selectFile in limitedModel){
    if(selectFile %in% c("LASSO_limited.txt", "EN_limited.txt", "RF_limited.txt", "XGBoost_limited.txt")){
      method <- strsplit(selectFile, '_limited')[[1]][1]
    } else {
      method <- strsplit(selectFile, '.txt')[[1]][1]
    }
    if(method %in% c('mRMR','fisher','FDC')){
      selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F,skipNul = T)
      pvlB <- c()
      for(ix in selectedList$V1){
        pvlB <- c(pvlB, pvllist$V3[pvllist$V2 == ix])
      }
      pvlB <- as.numeric(pvlB)
      
    }else{
      selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F,skipNul = T)
      pvlB <- as.numeric(pvllist$V3[pvllist$V2 %in% selectedList$V1])
      abB <- as.numeric(ablist$V3[ablist$V2 %in% selectedList$V1])
    }
    pvlTable <- getCD(cancerType, method, pvlAD, pvlB, pvlTable)
    abTable <- getCD(cancerType, method, abAD, abB, abTable)
  }
}

pvlTable$Data <- factor(pvlTable$Data, levels = sort(unique(pvlTable$Data)))
pvlTable$Method <- factor(pvlTable$Method, levels =  c("LASSO", "EN", "RF", "XGBoost", "MI", "mRMR",
                                                       "ReliefF", "fisher", "FDC"))
abTable$Data <- factor(abTable$Data, levels = sort(unique(abTable$Data)))
abTable$Method <- factor(abTable$Method, levels =  c("LASSO", "EN", "RF", "XGBoost", "MI", "mRMR",
                                                       "ReliefF", "fisher", "FDC"))
mycolor <- c("ADlasso" = "#E41A1C", "LASSO" = "#377EB8", "EN" = "#4DAF4A", "RF" = "#984EA3", "XGBoost" = "#FF7F00", "MI" = "#FFFF33",
             "mRMR" = "#A65628", "ReliefF" = "#F781BF", "fisher" = "#999999", "FDC" = "#9F88FF")
p1 <- ggplot(pvlTable, aes(x = Data, y = cohensD, color = Method)) + 
  ylab("Cohen's d \n (Prevalence)") + scale_color_manual(values = mycolor) +
  scale_x_discrete(labels= as.character(sort(unique(pvlTable$Data)))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title=element_blank(),
        legend.position = 'none')
sing = 1
for (i in 1:(length(unique(pvlTable$Data))-1)){
  sing = sing * -1
  p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, fill = ifelse(sing > 0, 'white', 'gray95'))
}
p1 <- p1 + geom_jitter(width = 0.4) + geom_hline(aes(yintercept = 0, lty = '0')) + geom_hline(aes(yintercept = 0.8, lty ='0.8')) +
      scale_linetype_manual(name="Cohen's d", values = c(2,3))

p2 <- ggplot(abTable, aes(x = Data, y = cohensD, color = Method)) + 
  ylab("Cohen's d\n (Abundance)") + scale_color_manual(values = mycolor) +
  scale_x_discrete(labels= as.character(sort(unique(pvlTable$Data)))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
sing = 1
for (i in 1:(length(unique(pvlTable$Data))-1)){
  sing = sing * -1
  p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, fill = ifelse(sing > 0, 'white', 'gray95'))
}
p2 <- p2 + geom_jitter(width = 0.4) + geom_hline(aes(yintercept = 0, lty = '0')) + geom_hline(aes(yintercept = 0.8, lty ='0.8')) +
      scale_linetype_manual(name="Cohen's d", values = c(2,3)) + guides(color = "none")




PCdf <- read.table('/home/yincheng23/miRNA/TCGA/script/performance_limitedModel.txt', sep = "\t", stringsAsFactors = F, header = T)
PCdf <- PCdf[PCdf$Matrix == "AUC", ]
PCdf <- PCdf[PCdf$Method != "SVMLASSO", ]
PCdf$Data <- factor(PCdf$Data, levels = sort(unique(pvlTable$Data)))
PCdf$Method <- factor(PCdf$Method, levels =  c("ADlasso","LASSO", "EN", "RF", "XGBoost", "MI", "mRMR",
                                               "ReliefF", "fisher", "FDC"))
PCdf$trans <- ifelse(PCdf$Method == 'ADlasso',2,1)
p3 <- ggplot(PCdf, aes(x = Data, y = value, color = Method)) + 
  ylab("AUC") + scale_color_manual(values = mycolor) + guides(alpha = "none") +
  scale_x_discrete(labels= as.character(sort(unique(pvlTable$Data)))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())
sing = 1
for (i in 1:(length(unique(pvlTable$Data))-1)){
  sing = sing * -1
  p3 <- p3 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, fill = ifelse(sing > 0, 'white', 'gray95'))
}
p3 <- p3 + geom_jitter(width = 0.25)

p <- p1/p2/p3 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')

pdf("/home/yincheng23/ADlasso_manuscript/Figure/FigureS6.pdf", width=8, height=6)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/FigureS6.tiff", units="in", width=8, height=6, res=300)
p
dev.off()