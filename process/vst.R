# +===================================================+
# |          check property of each dataset           |
# +===================================================+

wd <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies'
setwd(wd)
folder <- dir()
n_samples <- c()
n_features <- c()
datasets <- c()
for(dataset in folder){
  setwd(wd)
  print(dataset)
  if(substr(dataset,1,3) == 'ycc'){next}
  
  setwd(dataset)
  files <- dir()
  ASVsfilename <- grep('ASVs_table.tsv',files, value = T)
  metafilename <- grep('meta.tsv',files, value = T)
  
  con <- file(ASVsfilename)
  file_1_line1 <- readLines(con,n=1)
  close(con)
  
  if(grepl("Constructed from biom file", file_1_line1)){
    ASV_table <- read.table(ASVsfilename, sep="\t", skip=1, header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }else{
    ASV_table <- read.table(ASVsfilename, sep="\t", header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }
  
  groupings <- read.table(metafilename, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
  
  #number of samples
  AS <- (colnames(ASV_table))
  MS <- (rownames(groupings))
  inter <- intersect(AS,MS)
  ASu <- setdiff(AS,MS)
  MSu <- setdiff(MS,AS)
  print(paste(c(length(ASu), length(inter), length(MSu))))
  n_samples <- c(n_samples, length(inter))
  n_features <- c(n_features, nrow(ASV_table))
  datasets <- c(datasets, dataset)
  print('-----------------------------------------------------')
}

data_property <- data.frame(featureSize = n_features,  sampleSize = n_samples, name = datasets)
refer <- read.csv('/home/yincheng23/ADlasso_manuscript/Script/38list.csv')
setdiff(data_property$name, refer$Dataset.Name)
data_property <- data_property[data_property$name %in% refer$Dataset.Name, ]
data_property$source <- ""
data_property$platform <- ""
data_property$region <- ""
for(i in 1:nrow(data_property)){
  data_property$source[i] <- strsplit(refer$X[refer$Dataset.Name == data_property$name[i]],' ')[[1]][1]
  data_property$platform[i] <- refer$Seq..Tech.[refer$Dataset.Name == data_property$name[i]]
  data_property$region[i] <- refer$X16S.variable.region[refer$Dataset.Name == data_property$name[i]]
}

# +===================================================+
# |                VST transformation                 |
# +===================================================+

library(DESeq2)
setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies')
folder <- dir()
for(dataset in folder){
  print(dataset)
  if(substr(dataset,1,3) == 'ycc'){next}
  
  setwd(dataset)
  files <- dir()
  ASVsfilename <- grep('ASVs_table.tsv',files, value = T)
  metafilename <- grep('meta',files, value = T)
  
  con <- file(ASVsfilename)
  file_1_line1 <- readLines(con,n=1)
  close(con)
  
  if(grepl("Constructed from biom file", file_1_line1)){
    ASV_table <- read.table(ASVsfilename, sep="\t", skip=1, header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }else{
    ASV_table <- read.table(ASVsfilename, sep="\t", header=T, row.names = 1, 
                            comment.char = "", quote="", check.names = F)
  }
  
  groupings <- read.table(metafilename, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
  groupings$Class <- groupings[,1]
  ASV_sample <- (colnames(ASV_table))
  Meta_sample <- (rownames(groupings))
  intersect_sample <- intersect(ASV_sample,Meta_sample)
  
  ASV_table <- ASV_table[,intersect_sample]
  groupings <- groupings[intersect_sample,]
  
  n_feature <- nrow(ASV_table)
  ASV_table <- ASV_table[rownames(ASV_table)[rowSums(ASV_table) > 0],]
  n_feature_filt <- nrow(ASV_table)
  
  print(paste('sample :', length(ASV_sample), '->', length(intersect_sample)))
  print(paste('feature :', n_feature, '->', n_feature_filt))
  
  
  ASV_pseudo_count <- ASV_table + 1
  groupings$Class <- factor(groupings$Class)
  
  print('strat to DESeq2')
  dds <- DESeqDataSetFromMatrix(countData = ASV_pseudo_count,
                                colData = groupings,
                                ~ Class)
  print('strat to VST')
  vsd <- varianceStabilizingTransformation(dds, fitType="mean", blind= T)
  print('VST done')
  vsd_table <- assay(vsd)
  write.table(vsd_table, file="ASV_vst.txt", quote=F, sep="\t", row.names=T)
  write.table(ASV_table, file="ASV_table.txt", quote=F, sep="\t", row.names=T)
  write.table(groupings, file="metadata.txt", quote=F, sep="\t", row.names=T)
  print('-----------------------------------------------------')
  setwd('../')
}

# +===================================================+
# |          VST transformation in my data            |
# +===================================================+

library(DESeq2)
setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies')
folder <- dir()
for(dataset in folder){
  print(dataset)
  if(substr(dataset,1,3) != 'ycc'){next}
  
  setwd(dataset)
  files <- dir()
  ASVsfilename <- grep('ASV_table.txt',files, value = T)
  metafilename <- grep('metadata.txt',files, value = T)
  
  ASV_table <- read.table(ASVsfilename, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
  groupings <- read.table(metafilename, sep="\t", row.names = 1, header=T, comment.char = "", quote="", check.names = F)
  ASV_sample <- (colnames(ASV_table))
  Meta_sample <- (rownames(groupings))
  intersect_sample <- intersect(ASV_sample,Meta_sample)
  
  ASV_table <- ASV_table[,intersect_sample]
  groupings <- groupings[intersect_sample,]
  ASV_table <- ASV_table[rownames(ASV_table)[rowSums(ASV_table) > 0],]
  ASV_pseudo_count <- ASV_table + 1
  groupings$Class <- factor(groupings$Class)
  
  print('strat to DESeq2')
  dds <- DESeqDataSetFromMatrix(countData = ASV_pseudo_count,
                                colData = groupings,
                                ~ Class)
  print('strat to VST')
  vsd <- varianceStabilizingTransformation(dds, fitType="mean", blind= T)
  print('VST done')
  vsd_table <- assay(vsd)
  write.table(vsd_table, file="ASV_vst.txt", quote=F, sep="\t", row.names=T)
  write.table(ASV_table, file="ASV_table.txt", quote=F, sep="\t", row.names=T)
  write.table(groupings, file="metadata.txt", quote=F, sep="\t", row.names=T)
  print('-----------------------------------------------------')
  setwd('../')
}

# +===================================================+
# |                    Check taxa                     |
# +===================================================+

setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies')
folder <- dir()
i <- 40
dataset <- folder[i]
setwd(dataset)
files <- dir()
ASVsfilename <- grep('ASVs_table.tsv',files, value = T)
con <- file(ASVsfilename)
file_1_line1 <- readLines(con,n=1)
close(con)
if(grepl("Constructed from biom file", file_1_line1)){
  ASV_table <- read.table(ASVsfilename, sep="\t", skip=1, header=T, row.names = 1, 
                          comment.char = "", quote="", check.names = F)
}else{
  ASV_table <- read.table(ASVsfilename, sep="\t", header=T, row.names = 1, 
                          comment.char = "", quote="", check.names = F)
}
print(dataset)
print(rownames(ASV_table)[1:10])
print(dir())



need_annaotation_rowname <- c('ArcticFireSoils', 'ArcticFreshwaters','ArcticTransects','Ji_WTP_DS','Office')

need_annaotation_fasta   <- c('glass_plastic_oberbeckmann_ASVs','GWMC_ASIA_NA','GWMC_HOT_COLD','sed_plastic_hoellein',
                              'sed_plastic_rosato','seston_plastic_mccormick','sw_plastic_frere','sw_sed_detender',
                              'wood_plastic_kesy')

cannot_annaotation       <- c('art_scher','asd_son','BISCUIT','Blueberry','cdi_schubert',
                              'cdi_vincent','Chemerin','crc_baxter','crc_zeller','edd_singh',
                              'Exercise','hiv_dinh','hiv_lozupone','hiv_noguerajulian','ibd_gevers',
                              'ibd_papa','MALL','ob_goodrich','ob_ross','ob_turnbaugh',
                              'ob_zhu','ob_zupancic','par_scheperjans','t1d_alkanani','t1d_mejialeon')