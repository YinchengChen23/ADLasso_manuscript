library(edgeR)
library(phyloseq)
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

setwd('/home/yincheng23/ADlasso_manuscript/Datasets/Studies')
dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness'
outPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness/'
dataSets <- dir(outPath)
for(dataSet in dataSets){
  tag <- strsplit(dataSet, '_')[[1]]
  tag1 <- tag[1]; tag2 <- paste0(tag[2:length(tag)], collapse = '_')
  print(dataSet)
  setwd(dataPath); setwd(tag1); setwd(tag2)
  
  t1=proc.time()
  
  ASV_table <- read.table("ASV_table.txt", sep = "\t", stringsAsFactors = F, check.names = F)
  groupings <- read.table("metadata.txt", sep = "\t", stringsAsFactors = F)
  
  OTU <- phyloseq::otu_table(ASV_table, taxa_are_rows = T)
  sampledata <- phyloseq::sample_data(groupings, errorIfNULL = T)
  phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
  test <- phyloseq_to_edgeR(physeq = phylo, group="Class")
  et = exactTest(test)
  tt = topTags(et, n=nrow(test$table), adjust.method="BH", sort.by="PValue")
  res <- tt@.Data[[1]]
  sublist <- res[res$FDR < 0.001,]
  sublist <- cbind(rownames(sublist),sublist)
  
  outPath_dataSet <- paste0(outPath, dataSet)
  if(! file.exists(outPath_dataSet)){
    dir.create(outPath_dataSet)
  }
  outPath_dataSet_selectedList <- paste0(outPath, dataSet, '/','edgeR', '.txt')
  write.table(sublist, outPath_dataSet_selectedList, quote=FALSE, sep="\t", col.names = F, row.names = F)
  
  t2=proc.time()
  tCost=t2-t1
  print(paste0('The cost：',tCost[3][[1]]/60,'min'))
  print('-----------------------------------------------------')
  setwd('../')
}

#Rscript --vanilla runedgeR.R
