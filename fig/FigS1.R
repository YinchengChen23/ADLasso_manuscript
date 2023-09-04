library(effsize)
library(patchwork)
library(RColorBrewer)

getPVL <- function(data, rawpvl, List){
  pvlL <- rawpvl[rownames(data) %in% List$V1]
  return(pvlL)
}
getCD <- function(data, rawpvl, List, refList, DataSet, Method, addDF){
  pvl <- getPVL(data, rawpvl, List)
  refpvl <- getPVL(data, rawpvl, refList)
  if(length(pvl)>1){
    cohD <- cohen.d(refpvl, pvl)
    cohD <- cohD$estimate
  } else if(length(pvl) == 1) {
    stds = ((length(rawpvl)-1)* (sd(rawpvl)**2))/(length(rawpvl)-1)
    stds = stds**(0.5)
    cohD <- (pvl - mean(refpvl))/stds
  } else {
    cohD <- NA
  }
  tmp <- data.frame("Data" = DataSet, "Method" = Method, "cohensD" = cohD, "number" = nrow(List))
  addDF <- rbind(addDF, tmp)
  return(addDF)
}

#===========================================#
#           unlimitedModel                  #
#===========================================#

pvlTable <-  data.frame()
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results'
selectList <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info'
setwd(selectList)
unlimitedModel = c("LASSO.txt", "EN.txt", "RF.txt", "XGBoost.txt")
dataSets <- dir()
for(dataSet in dataSets){
  if(dataSet == 'LIBSVM_rs'){next}
  #----------- read raw data --------------
  dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/'
  dataPath_ <- paste0(dataPath, dataSet, '/ASV_table.txt')
  data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)
  data[data > 0] <- 1
  prevalence <- rowSums(data)/ncol(data)*100
  #----------- including different method --------------
  setwd(refpath)
  setwd(dataSet)
  selectedList_AD <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F)
  setwd(selectList)
  setwd(dataSet)
  print(dataSet)
  for(selectFile in unlimitedModel){
    method <- strsplit(selectFile, '.txt')[[1]][1]
    check <- readLines(selectFile)
    if(length(check) == 0){next}
    selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F,skipNul = T)
    pvlTable <- getCD(data, prevalence, selectedList, selectedList_AD, dataSet, method, pvlTable)
  }
}
pvlTable$Data <- factor(pvlTable$Data, levels = sort(unique(pvlTable$Data)))
pvlTable$Method <- factor(pvlTable$Method, levels =  c("LASSO", "EN", "RF", "XGBoost"))
mycolor <- c("ADlasso" = "#E41A1C", "LASSO" = "#377EB8", "EN" = "#4DAF4A", "RF" = "#984EA3", "XGBoost" = "#FF7F00")
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
  scale_linetype_manual(values = c(2,3))

abTable <-  data.frame()
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results'
selectList <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info'
setwd(selectList)
dataSets <- dir()
for(dataSet in dataSets){
  #----------- read raw data --------------
  if(dataSet == 'LIBSVM_rs'){next}
  dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/'
  dataPath_ <- paste0(dataPath, dataSet, '/ASV_table.txt')
  data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)
  ab <- rowSums(data)/sum(data)*100
  #----------- including different method --------------
  setwd(refpath)
  setwd(dataSet)
  selectedList_AD <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F)
  setwd(selectList)
  setwd(dataSet)
  print(dataSet)
  for(selectFile in unlimitedModel){
    method <- strsplit(selectFile, '.txt')[[1]][1]
    check <- readLines(selectFile)
    if(length(check) == 0){next}
    selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F,skipNul = T)
    abTable <- getCD(data, ab, selectedList, selectedList_AD, dataSet, method, abTable)
  }
}
abTable$Data <- factor(abTable$Data, levels = sort(unique(abTable$Data)))
abTable$Method <- factor(abTable$Method, levels =  c("LASSO", "EN", "RF", "XGBoost"))
p2 <- ggplot(abTable, aes(x = Data, y = cohensD, color = Method)) + 
  ylab("Cohen's d\n (Abundance)") + scale_color_manual(values = mycolor) +
  scale_x_discrete(labels= as.character(sort(unique(abTable$Data)))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        axis.text.x = element_blank(),
        axis.title.x = element_blank())
sing = 1
for (i in 1:(length(unique(abTable$Data))-1)){
  sing = sing * -1
  p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, fill = ifelse(sing > 0, 'white', 'gray95'))
}
p2 <- p2 + geom_jitter(width = 0.4) + geom_hline(aes(yintercept = 0, lty = '0')) + geom_hline(aes(yintercept = 0.8, lty ='0.8')) +
  scale_linetype_manual(name="Cohen's d",values = c(2,3)) + guides(color = "none")



PCdf <- read.table("/home/yincheng23/ADlasso_manuscript/Script/performance_ML.txt", sep = "\t", stringsAsFactors = F, header = T)
PCdf <- PCdf[PCdf$Matrix == "AUC", ]
PCdf <- PCdf[PCdf$Method != "SVMLASSO", ]
PCdf$Data <- factor(PCdf$Data, levels = sort(unique(pvlTable$Data)))
PCdf$Method <- factor(PCdf$Method, levels =  c("ADlasso","LASSO", "EN", "RF", "XGBoost"))
PCdf$trans <- ifelse(PCdf$Method == 'ADlasso',2,1.5)
p3 <- ggplot(PCdf, aes(x = Data, y = value, color = Method)) + 
  ylab("AUC") + scale_color_manual(values = mycolor) +
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
  p3 <- p3 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, fill = ifelse(sing > 0, 'white', 'gray95'))
}
p3 <- p3 + geom_jitter(width = 0.4)

numTable <-  data.frame()
setwd(selectList)
unlimitedModel = c("LASSO.txt", "EN.txt", "RF.txt", "XGBoost.txt")
dataSets <- dir()
for(dataSet in dataSets){
  if(dataSet == 'LIBSVM_rs'){next}
  setwd(refpath)
  setwd(dataSet)
  selectedList_AD <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F)
  tmp <- data.frame("Data" = dataSet, "Method" = 'ADlasso', "number" = nrow(selectedList_AD))
  numTable <- rbind(numTable, tmp)
  setwd(selectList)
  setwd(dataSet)
  print(dataSet)
  for(selectFile in unlimitedModel){
    method <- strsplit(selectFile, '.txt')[[1]][1]
    check <- readLines(selectFile)
    if(length(check) == 0){next}
    selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F)
    tmp <- data.frame("Data" = dataSet, "Method" = method, "number" = nrow(selectedList))
    numTable <- rbind(numTable, tmp)
  }
}
numTable$Data <- factor(numTable$Data, levels = sort(unique(numTable$Data)))
numTable$Method <- factor(numTable$Method, levels =  c("ADlasso","LASSO", "EN", "RF", "XGBoost"))
numTable$number <- log(numTable$number,10)
p4 <- ggplot(numTable, aes(x = Data, y = number, color = Method)) + 
  ylab("log10(#Features)") + scale_color_manual(values = mycolor) +
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
for (i in 1:(length(unique(numTable$Data))-1)){
  sing = sing * -1
  p4 <- p4 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf, fill = ifelse(sing > 0, 'white', 'gray95'))
}
p4 <- p4 + geom_jitter(width = 0.4) + geom_hline(aes(yintercept = 1, lty = '10')) +
  scale_linetype_manual(name='Number of feature',values = 4)

p <- p1/p2/p3/p4 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')


pdf("/home/yincheng23/ADlasso_manuscript/Figure/FigureS1.pdf", width=10, height=9)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/FigureS1.tiff", units="in", width=10, height=9, res=300)
p
dev.off()
