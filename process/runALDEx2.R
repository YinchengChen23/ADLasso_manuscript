library(ALDEx2)
setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies')
outPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results/'
bigDataSets <- c('e')
alreadyDown <- c('n')
dataSets <- c('ArcticTransects')
for(dataSet in dataSets){
  if(dataSet %in% bigDataSets){next}
  if(dataSet %in% alreadyDown){next}
  print(dataSet)
  setwd(dataSet)
  t1=proc.time()
  
  ASV_table <- read.table("ASV_table.txt", sep = "\t", stringsAsFactors = F)
  groupings <- read.table("metadata.txt", sep = "\t", stringsAsFactors = F)
  
  results <- aldex(reads=ASV_table, conditions = groupings$Class, mc.samples = 128, test="t", effect=TRUE,
                   include.sample.summary = FALSE, verbose=T, denom="all")
  sublist <- results[results$wi.ep < 0.05,]
  sublist <- cbind(rownames(sublist),sublist)
  
  outPath_dataSet <- paste0(outPath, dataSet)
  if(! file.exists(outPath_dataSet)){
    dir.create(outPath_dataSet)
  }
  outPath_dataSet_selectedList <- paste0(outPath, dataSet, '/','ALDEx2', '.txt')
  write.table(sublist, outPath_dataSet_selectedList, quote=FALSE, sep="\t", col.names = F, row.names = F)
  
  t2=proc.time()
  tCost=t2-t1
  print(paste0('The cost：',tCost[3][[1]]/60,'min'))
  print('-----------------------------------------------------')
  setwd('../')
}

#Rscript --vanilla runALDEx2.R
