library(effsize)
library(patchwork)
library(RColorBrewer)

getPVL <- function(data, rawpvl, List){
  pvlL <- rawpvl[rownames(data) %in% List$V1]
  return(as.vector(pvlL))
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

pvlTable <-  data.frame()
selectList <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results'
setwd(selectList)
dataSets <- dir()
for(dataSet in dataSets){
  setwd(selectList)
  setwd(dataSet)
  selectFiles <- dir()
  if(!'ADlasso.txt' %in% selectFiles){next}
  if(dataSet == 'LIBSVM_rs'){next}
  #----------- read raw data --------------
  dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/'
  dataPath_ <- paste0(dataPath, dataSet, '/ASV_table.txt')
  data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)  # data <- read.table(paste0(dataPath, dataSet, '/metadata.txt'), sep = "\t", stringsAsFactors = F)
  data[data > 0] <- 1
  prevalence <- rowSums(data)/ncol(data)*100
  #----------- including different method --------------
  selectedList_AD <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F)
  for(selectFile in selectFiles){
    con <- file(selectFile,open="r")
    check_tag <- readLines(con)
    close(con)
    if(length(check_tag) == 0){print(dataSet);print(selectFile);next}
    if(selectFile == 'ADlasso.txt'){next}
    method <- strsplit(selectFile, '.txt')[[1]][1]
    selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F,skipNul = T)
    pvlTable <- getCD(data, prevalence, selectedList, selectedList_AD, dataSet, method, pvlTable)
  }
}
pvlTable$Data <- factor(pvlTable$Data, levels = sort(unique(pvlTable$Data)))
pvlTable$Method <- factor(pvlTable$Method, levels = c("ALDEx2", "ANCOM2", "LEfSe", "edgeR"))
mycolor <- c("ADlasso" = "#E41A1C", "ALDEx2" = "#377EB8", "ANCOM2" = "#4DAF4A", "LEfSe" = "#984EA3", "edgeR" = "#FF7F00")
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
p1 <- p1 + geom_jitter(width = 0.25) + geom_hline(aes(yintercept = 0, lty = '0')) + geom_hline(aes(yintercept = 0.8, lty ='0.8')) +
      scale_linetype_manual(values = c(2,3))



abTable <-  data.frame()
setwd(selectList)
dataSets <- dir()
for(dataSet in dataSets){
  setwd(selectList)
  setwd(dataSet)
  selectFiles <- dir()
  if(!'ADlasso.txt' %in% selectFiles){next}
  if(dataSet == 'LIBSVM_rs'){next}
  #----------- read raw data --------------
  dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/'
  dataPath_ <- paste0(dataPath, dataSet, '/ASV_table.txt')
  data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)
  ab <- rowSums(data)/sum(data)*100
  #----------- including different method --------------
  selectedList_AD <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F)
  for(selectFile in selectFiles){
    con <- file(selectFile,open="r")
    check_tag <- readLines(con)
    close(con)
    if(length(check_tag) == 0){print(dataSet);print(selectFile);next}
    if(selectFile == 'ADlasso.txt'){next}
    method <- strsplit(selectFile, '.txt')[[1]][1]
    selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F,skipNul = T)
    abTable <- getCD(data, ab, selectedList, selectedList_AD, dataSet, method, abTable)
  }
  print(dataSet)
}
abTable$Data <- factor(abTable$Data, levels = sort(unique(abTable$Data)))
abTable$Method <- factor(abTable$Method, levels = c("ALDEx2", "ANCOM2", "LEfSe", "edgeR"))
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
p2 <- p2 + geom_jitter(width = 0.25) + geom_hline(aes(yintercept = 0, lty = '0')) + geom_hline(aes(yintercept = 0.8, lty ='0.8')) +
      scale_linetype_manual(name="Cohen's d",values = c(2,3)) + guides(color = "none")

PCdf <- read.table("/home/yincheng23/ADlasso_manuscript/Script/performance_statistics.txt", sep = "\t", stringsAsFactors = F, header = T)
PCdf <- PCdf[PCdf$Matrix == "AUC", ]
PCdf$Data <- factor(PCdf$Data, levels = sort(unique(pvlTable$Data)))
PCdf$Method <- factor(PCdf$Method, levels = c("ADlasso", "ALDEx2", "ANCOM2", "LEfSe", "edgeR"))
PCdf$trans <- ifelse(PCdf$Method == 'ADlasso',2,1)

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
p3 <- p3 + geom_jitter(width = 0.25)


# sw_plastic_frere ANCOM2  
# ob_zupancic ANCOM2

numTable <-  data.frame()
setwd(selectList)
dataSets <- dir()
for(dataSet in dataSets){
  setwd(selectList)
  setwd(dataSet)
  selectFiles <- dir()
  if(!'ADlasso.txt' %in% selectFiles){next}
  if(dataSet == 'LIBSVM_rs'){next}
  for(selectFile in selectFiles){
    con <- file(selectFile,open="r")
    check_tag <- readLines(con)
    close(con)
    if(length(check_tag) == 0){
      print(dataSet)
      print(selectFile)
      method <- strsplit(selectFile, '.txt')[[1]][1]
      tmp <- data.frame("Data" = dataSet, "Method" = method, "number" = 0)
      numTable <- rbind(numTable, tmp)
      print('-------------')
    }else{
      method <- strsplit(selectFile, '.txt')[[1]][1]
      selectedList <- read.table(selectFile, sep = "\t", stringsAsFactors = F)
      tmp <- data.frame("Data" = dataSet, "Method" = method, "number" = nrow(selectedList))
      numTable <- rbind(numTable, tmp)
    }
  }
}
numTable$Data <- factor(numTable$Data, levels = sort(unique(numTable$Data)))
numTable$Method <- factor(numTable$Method, levels = c("ADlasso", "ALDEx2", "ANCOM2", "LEfSe", "edgeR"))
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
p4 <- p4 + geom_jitter(width = 0.25) + geom_hline(aes(yintercept = 1, lty = '10')) +
      scale_linetype_manual(name='Number of feature',values = 4)

p <- p1/p2/p3/p4 + plot_layout(guides = 'collect') + plot_annotation(tag_levels = 'A')
p

pdf("/home/yincheng23/ADlasso_manuscript/Figure/Figure2.pdf", width=10, height=9)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/Figure2.tiff", units="in", width=10, height=9, res=300)
p
dev.off()




# +===================================================+
# |                Diarrhea Testing                   |
# +===================================================+

pvlCTable <-  data.frame()
RawDataSetsPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies'
RBDDataSetsPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness/Diarrhea'
setwd(RBDDataSetsPath)
RBDDataSets <- dir()
RawDataSets <- dir(RawDataSetsPath)
tagList <- RawDataSets[RawDataSets %in% RBDDataSets]

par(mfrow = c(1, 3))
singleData = tagList[3]
for(singleData in tagList){
    setwd(RBDDataSetsPath)
    setwd(singleData)
    
    dataPath_ <- paste0(RBDDataSetsPath, '/', singleData, '/', singleData, '_genus_table.tsv')
    data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F, row.names = 1)
    data[data > 0] <- 1
    RBDprevalence <- rowSums(data)/ncol(data)*100
    
    dataPath_ <- paste0(RawDataSetsPath, '/', singleData, '/', 'ASV_table.txt')
    data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)
    data[data > 0] <- 1
    Rawprevalence <- rowSums(data)/ncol(data)*100
      
    cohD <- cohensD(Rawprevalence, RBDprevalence)
    tmp <- data.frame("Data" = singleData, "cohensD" = cohD,
                      'Raw' = median(Rawprevalence), "RBD" = median(RBDprevalence))
    pvlCTable <- rbind(pvlCTable, tmp)
}


kd1 <- density(Rawprevalence)
plot(kd1, col='green', lwd=2,  main = paste("Data :", tagList[1]), xlab = 'prevalence')
kd2 <- density(RBDprevalence)
lines(kd2, col='red', lwd=2)

# +======================================================================+
# |                                                                      |
# |   Density curves with gradient fill colors along the x axis          |
# |                                                                      |
# +======================================================================+

selectList <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info'
setwd(selectList)
dataSets <- dir()
pvlTable <-  data.frame()
for(dataSet in dataSets){
  #----------- read raw data --------------
  dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/'
  dataPath_ <- paste0(dataPath, dataSet, '/ASV_table.txt')
  data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)
  data[data > 0] <- 1
  prevalence <- rowSums(data)/ncol(data)*100
  tmp <- data.frame('Dataset'=dataSet, 'prevalence'=prevalence)
  pvlTable <- rbind(pvlTable, tmp)
  print(dataSet)
}

PCdf <- read.table("/home/yincheng23/ADlasso_manuscript/Script/performance_MLlimited.txt", sep = "\t", stringsAsFactors = F, header = T)
PCdf <- PCdf[PCdf$Matrix == "AUC", ]
PCdf <- PCdf[PCdf$Method == "ADlasso", ]

ggplot(pvlTable, aes(x = prevalence, y = Dataset, fill = stat(x))) +
  geom_density_ridges_gradient(color = 'white', scale = 3, size = 0.5, rel_min_height = 0.01)+
  scale_fill_viridis_c(option = "plasma", name = "Prevalence (%)")



ggplot(pvlTable, aes(x = `prevalence`, y = `Dataset`, height = stat(density))) +
  geom_density_ridges(stat = "binline", bins = 50,draw_baseline = T) 
  geom_density_ridges_gradient(PCdf, aes(y = `Dataset`, fill = value),scale = 3, size = 0.3, rel_min_height = 0.01)



# +===================================================+
# |               Obesity Testing                     |
# +===================================================+ 


RawDataSetsPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies'
RBDDataSetsPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness/Obesity'
RawDataSets <- dir(RawDataSetsPath)
RBDDataSets <- dir(RBDDataSetsPath)
par(mfrow = c(2, 4))

Rawvec <- c()
RBDvec <- c()
for(x in RBDDataSets){
  taget <- grep(x, RawDataSets, value = T)
  if(length(taget) == 0){next}
  Rawvec <- c(Rawvec, taget)
  RBDvec <- c(RBDvec, x)
}
pvlCTable <- data.frame()
for(i in 1:length(RBDvec)){
  dataPath_ <- paste0(RBDDataSetsPath, '/', RBDvec[i])
  setwd(dataPath_)
  fileslist <- dir()
  get <- grep('metadata.tsv', fileslist, value = T)
  data <- read.table(get, sep = "\t", stringsAsFactors = F, row.names = 1, header = T)
  data[data > 0] <- 1
  RBDprevalence <- rowSums(data)/ncol(data)*100
  
  dataPath_ <- paste0(RawDataSetsPath, '/', Rawvec[i], '/', 'ASV_table.txt')
  data <- read.table(dataPath_, sep = "\t", stringsAsFactors = F)
  data[data > 0] <- 1
  Rawprevalence <- rowSums(data)/ncol(data)*100
  
  cohD <- cohensD(Rawprevalence, RBDprevalence)
  tmp <- data.frame("Data" = Rawvec[i], "cohensD" = cohD,
                    'Raw' = median(Rawprevalence), "RBD" = median(RBDprevalence))
  pvlCTable <- rbind(pvlCTable, tmp)
}

kd1 <- density(Rawprevalence)
plot(kd1, col='green', lwd=2,  main = paste("Data :", RBDvec[i]), xlab = 'prevalence')
kd2 <- density(RBDprevalence)
lines(kd2, col='red', lwd=2)
