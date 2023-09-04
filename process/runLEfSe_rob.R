dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness'
target <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness/'
outPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/LEfse_tmp/'
dataSets <- dir(target)
for(dataSet in dataSets){
  tag <- strsplit(dataSet, '_')[[1]]
  tag1 <- tag[1]; tag2 <- paste0(tag[2:length(tag)], collapse = '_')
  print(dataSet)
  setwd(dataPath); setwd(tag1); setwd(tag2)

  t1=proc.time()
  
  ASV_table <- read.table("ASV_table.txt", sep = "\t", stringsAsFactors = F)
  groupings <- read.table("metadata.txt", sep = "\t", stringsAsFactors = F)
  
  
  flip_ASV_table <- data.frame(t(ASV_table), check.names = F)
  flip_ASV_table <- cbind(groupings$Class, flip_ASV_table)
  colnames(flip_ASV_table)[1] <- "Class"
  flip_ASV_table <- cbind(rownames(groupings), flip_ASV_table)
  colnames(flip_ASV_table)[1] <- "id"
  ret_tab <- data.frame(t(flip_ASV_table), check.names = F)
  
  
  outPath_dataSet <- paste0(outPath, dataSet)
  if(! file.exists(outPath_dataSet)){
    dir.create(outPath_dataSet)
  }
  outPath_dataSet_selectedList <- paste0(outPath, dataSet, '/','LEfSe.tmp')
  write.table(ret_tab, outPath_dataSet_selectedList, quote=FALSE, sep="\t", col.names = F)
  
  t2=proc.time()
  tCost=t2-t1
  print(paste0('The cost：',tCost[3][[1]]/60,'min'))
  print('-----------------------------------------------------')
  setwd('../')
}


#docker run -it --rm -v /home/yincheng23/ADlasso_manuscript/Datasets/data:/tmp yincheng23/lefse:0.0.4
#cd /tmp/LEfse_tmp
#DataSets=$(ls -l ./ |awk '/^d/ {print $NF}')
#for datsSet in $DataSets
#do
#   echo $datsSet
#   cd $datsSet
#   format_input.py LEfSe.tmp data_in -c 2 -u 1 -o 1000000
#   run_lefse.py data_in LEfSe_raw.txt
#   cat LEfSe_raw.txt | awk '{if($3>2){print $0}}' > LEfSe.txt
#   cd ..
#done




#cd /home/yincheng23/ADlasso_manuscript/Datasets/data/LEfse_tmp
#Folder=$(ls -l ./ |awk '/^d/ {print $NF}')
#for data in $Folder
#do
#  cd $data
#  echo $data
#  if [[ -f LEfSe.txt ]]; then
#  cp LEfSe.txt /home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results/$data/LEfSe.txt
#  fi
#  cd ..
#done




