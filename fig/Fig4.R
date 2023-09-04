library(patchwork)
checkpoint <- function(selectFile){
  con <- file(selectFile,open="r")
  check_tag <- readLines(con)
  close(con)
  if(length(check_tag) == 0){
    return(FALSE)
  } else {
    return(TRUE)
  }
}

extract_genus  <- function(taxalist){
  taxalist <- unique(taxalist)
  if(grepl('__', taxalist[1])){
    taxa_genus <- c()
    for(t in taxalist){
      splitlist <- strsplit(t, 'g__')[[1]]
      if(length(splitlist) > 1){
        subg <- splitlist[2]
        taxa_genus <- c(taxa_genus, subg)
      }
    }
    return(taxa_genus)
  } else {
    return(taxalist)
  }
}

fisher <- function(d1gs, d2gs, d1g, d2g){
  t1 <- intersect(d1gs, d2gs)
  t2 <- intersect(d1gs ,setdiff(d2g, d2gs))
  t3 <- intersect(d2gs ,setdiff(d1g, d1gs))
  t4 <- unique(c(setdiff(d1g, d1gs), setdiff(d2g, d2gs)))
  m = matrix(c(length(t1),length(t2),length(t3),length(t4)), nrow = 2)
  res <- fisher.test(m)
  return(res)
}

begdf <- data.frame()
datapath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Testing_Bias_robustness'
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness'
setwd(refpath)
for(x in dir()){
  dis <- strsplit(x, '_')[[1]][1]
  dataset <- strsplit(x, '_')[[1]]
  dataset <- paste(dataset[2:length(dataset)], collapse = "_")
  setwd(datapath); setwd(dis); setwd(dataset)
  data <- read.table("ASV_table.txt", sep = "\t", stringsAsFactors = F)
  gtotal <- extract_genus(rownames(data))
  for(z in gtotal){
    tmp <- data.frame(disease = dis, dataset= dataset, genus = z)
    begdf <- rbind(begdf, tmp)
  }
}


selectdf <- data.frame()
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness'
setwd(refpath)
for(x in dir()){
  dis <- strsplit(x, '_')[[1]][1]
  dataset <- strsplit(x, '_')[[1]]
  dataset <- paste(dataset[2:length(dataset)], collapse = "_")
  setwd(refpath); setwd(x)
  for(y in dir()){
    method <- strsplit(y, '.txt')[[1]][1]
    if(checkpoint(y)){
      selectedList <- read.table(y, sep = "\t", stringsAsFactors = F,skipNul = T)
      glist <- extract_genus(selectedList[,1])
      for(z in glist){
        tmp <- data.frame(disease = dis, dataset= dataset, method = method, genus = z)
        selectdf <- rbind(selectdf, tmp)
      }
    }
  }
}

subselectdf <- selectdf[selectdf$disease == 'Diarrhea',]
methods <- unique(subselectdf$method)
datasets <- unique(subselectdf$dataset)
pldf <- data.frame()
for(m in methods){
  for(d1 in  datasets){
    for(d2 in datasets){
      if(d1 != d2){
        d1g <- unique(begdf$genus[begdf$dataset == d1])
        d2g <- unique(begdf$genus[begdf$dataset == d2])
        d1gs <- unique(subselectdf$genus[subselectdf$dataset == d1 & subselectdf$method == m])
        d2gs <- unique(subselectdf$genus[subselectdf$dataset == d2 & subselectdf$method == m])
        FS <- fisher(d1gs, d2gs, d1g, d2g)
        JC <- length(intersect(d1gs, d2gs))/length(union(d1gs, d2gs))
        #JC <- length(intersect(d1gs, d2gs))/max(length(unique(d1gs)), length(unique(d2gs)))
        tmp <- data.frame(method = m, number = length(intersect(d1gs, d2gs)), JC = JC, odds = as.numeric(FS$estimate), p = FS$p.value)
        pldf <- rbind(pldf, tmp)
      }
    }
  }
}

subselectdf <- selectdf[selectdf$disease == 'Obesity',]
methods <- unique(subselectdf$method)
datasets <- unique(subselectdf$dataset)
pldf1 <- data.frame()
for(m in methods){
  for(d1 in  datasets){
    for(d2 in datasets){
      if(d1 != d2){
        d1g <- unique(begdf$genus[begdf$dataset == d1])
        d2g <- unique(begdf$genus[begdf$dataset == d2])
        d1gs <- unique(subselectdf$genus[subselectdf$dataset == d1 & subselectdf$method == m])
        d2gs <- unique(subselectdf$genus[subselectdf$dataset == d2 & subselectdf$method == m])
        FS <- fisher(d1gs, d2gs, d1g, d2g)
        JC <- length(intersect(d1gs, d2gs))/length(union(d1gs, d2gs))
        #JC <- length(intersect(d1gs, d2gs))/max(length(unique(d1gs)), length(unique(d2gs)))
        tmp <- data.frame(method = m, number = length(intersect(d1gs, d2gs)), JC = JC, odds = as.numeric(FS$estimate), p = FS$p.value)
        pldf1 <- rbind(pldf1, tmp)
      }
    }
  }
}
mycolor <- c("ADlasso" = "#E41A1C", "ALDEx2" = "#377EB8", "ANCOM2" = "#4DAF4A", "LEfSe" = "#984EA3", "edgeR" = "#FF7F00")
p1 <- ggplot(pldf, aes(x = JC, y = odds, color = method, size = number)) + geom_point(alpha=0.5) +
  labs(x = 'Jaccard similarity', y = 'Odds ratio', color = 'Method', size = 'Overlapping genera') +
  scale_color_manual(values = mycolor) + ggtitle('Diarrhea') +
  scale_size(breaks = range(c(pldf$number, pldf1$number))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        legend.position = 'none')
p2 <- ggplot(pldf1, aes(x = JC, y = odds, color = method, size = number)) + geom_point(alpha=0.5) +
  labs(x = 'Jaccard similarity', y = 'Odds ratio', color = 'Method', size = 'Overlapping genera') +
  scale_color_manual(values = mycolor) + ggtitle('Obesity') + 
  scale_size_continuous(breaks = range(c(pldf$number, pldf1$number))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10))
p1 + p2
#ggplot(pldf, aes(x = number, y = JC, color = method, size = JC)) + geom_point(alpha=0.3) +
#  labs(x = 'Overlapping genera', y = 'Jaccard index', color = 'Method', size = 'Jaccard index') +
#  scale_color_manual(values = mycolor)

selectdf <- data.frame()
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info_rob'
refadpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness'
limitedModel = c("LASSO_limited.txt", "EN_limited.txt", "RF_limited.txt", "XGBoost_limited.txt", "MI.txt", "mRMR.txt", "ReliefF.txt", "fisher.txt", "FDC.txt")
setwd(refpath)
for(x in dir()){
  dis <- strsplit(x, '_')[[1]][1]
  dataset <- strsplit(x, '_')[[1]]
  dataset <- paste(dataset[2:length(dataset)], collapse = "_")
  #------------- loading ADlasso -------------
  setwd(refadpath); setwd(x)
  selectedList <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F,skipNul = T)
  glist <- extract_genus(selectedList[,1])
  for(z in glist){
    tmp <- data.frame(disease = dis, dataset= dataset, method = 'ADlasso', genus = z)
    selectdf <- rbind(selectdf, tmp)
  }
  #------------- loading other -------------
  setwd(refpath); setwd(x)
  for(y in dir()){
    if(y %in% limitedModel){
      if(y %in% c("LASSO_limited.txt", "EN_limited.txt", "RF_limited.txt", "XGBoost_limited.txt")){
        method <- strsplit(y, '_limited')[[1]][1]
      } else {
        method <- strsplit(y, '.txt')[[1]][1]
      }
      selectedList <- read.table(y, sep = "\t", stringsAsFactors = F,skipNul = T)
      glist <- extract_genus(selectedList[,1])
      for(z in glist){
        tmp <- data.frame(disease = dis, dataset= dataset, method = method, genus = z)
        selectdf <- rbind(selectdf, tmp)
      }
    }
  }
}
subselectdf <- selectdf[selectdf$disease == 'Diarrhea',]
methods <- unique(subselectdf$method)
datasets <- unique(subselectdf$dataset)
pldf <- data.frame()
for(m in methods){
  for(d1 in  datasets){
    for(d2 in datasets){
      if(d1 != d2){
        d1g <- unique(begdf$genus[begdf$dataset == d1])
        d2g <- unique(begdf$genus[begdf$dataset == d2])
        d1gs <- unique(subselectdf$genus[subselectdf$dataset == d1 & subselectdf$method == m])
        d2gs <- unique(subselectdf$genus[subselectdf$dataset == d2 & subselectdf$method == m])
        FS <- fisher(d1gs, d2gs, d1g, d2g)
        JC <- length(intersect(d1gs, d2gs))/length(union(d1gs, d2gs))
        #JC <- length(intersect(d1gs, d2gs))/max(length(unique(d1gs)), length(unique(d2gs)))
        tmp <- data.frame(method = m, number = length(intersect(d1gs, d2gs)), JC = JC, odds = as.numeric(FS$estimate), p = FS$p.value)
        pldf <- rbind(pldf, tmp)
      }
    }
  }
}

subselectdf <- selectdf[selectdf$disease == 'Obesity',]
methods <- unique(subselectdf$method)
datasets <- unique(subselectdf$dataset)
pldf1 <- data.frame()
for(m in methods){
  for(d1 in  datasets){
    for(d2 in datasets){
      if(d1 != d2){
        d1g <- unique(begdf$genus[begdf$dataset == d1])
        d2g <- unique(begdf$genus[begdf$dataset == d2])
        d1gs <- unique(subselectdf$genus[subselectdf$dataset == d1 & subselectdf$method == m])
        d2gs <- unique(subselectdf$genus[subselectdf$dataset == d2 & subselectdf$method == m])
        FS <- fisher(d1gs, d2gs, d1g, d2g)
        JC <- length(intersect(d1gs, d2gs))/length(union(d1gs, d2gs))
        #JC <- length(intersect(d1gs, d2gs))/max(length(unique(d1gs)), length(unique(d2gs)))
        tmp <- data.frame(method = m, number = length(intersect(d1gs, d2gs)), JC = JC, odds = as.numeric(FS$estimate), p = FS$p.value)
        pldf1 <- rbind(pldf1, tmp)
      }
    }
  }
}
mycolor <- c("ADlasso" = "#E41A1C", "LASSO" = "#377EB8", "EN" = "#4DAF4A", "RF" = "#984EA3", "XGBoost" = "#FF7F00", "MI" = "#FFFF33",
             "mRMR" = "#A65628", "ReliefF" = "#F781BF", "fisher" = "#999999", "FDC" = "#9F88FF")
p3 <- ggplot(pldf, aes(x = JC, y = odds, color = method, size = number)) + geom_point(alpha=0.5) +
  labs(x = 'Jaccard similarity', y = 'Odds ratio', color = 'Method', size = 'Overlapping genera') +
  scale_color_manual(values = mycolor) + ggtitle('Diarrhea') +
  scale_size(breaks = range(c(pldf$number, pldf1$number))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        legend.position = 'none')
p4 <- ggplot(pldf1, aes(x = JC, y = odds, color = method, size = number)) + geom_point(alpha=0.5) +
  labs(x = 'Jaccard similarity', y = 'Odds ratio', color = 'Method', size = 'Overlapping genera') +
  scale_color_manual(values = mycolor) + ggtitle('Obesity') + 
  scale_size_continuous(breaks = range(c(pldf$number, pldf1$number))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10))

selectdf <- data.frame()
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/ML_parm_info_rob'
refadpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results_robustness'
unlimitedModel = c("LASSO.txt", "EN.txt", "RF.txt", "XGBoost.txt")
setwd(refpath)
for(x in dir()){
  dis <- strsplit(x, '_')[[1]][1]
  dataset <- strsplit(x, '_')[[1]]
  dataset <- paste(dataset[2:length(dataset)], collapse = "_")
  #------------- loading ADlasso -------------
  setwd(refadpath); setwd(x)
  selectedList <- read.table('ADlasso.txt', sep = "\t", stringsAsFactors = F,skipNul = T)
  glist <- extract_genus(selectedList[,1])
  for(z in glist){
    tmp <- data.frame(disease = dis, dataset= dataset, method = 'ADlasso', genus = z)
    selectdf <- rbind(selectdf, tmp)
  }
  #------------- loading other -------------
  setwd(refpath); setwd(x)
  for(y in dir()){
    if(y %in% unlimitedModel){
      if(y %in% c("LASSO_limited.txt", "EN_limited.txt", "RF_limited.txt", "XGBoost_limited.txt")){
        method <- strsplit(y, '_limited')[[1]][1]
      } else {
        method <- strsplit(y, '.txt')[[1]][1]
      }
      selectedList <- read.table(y, sep = "\t", stringsAsFactors = F,skipNul = T)
      glist <- extract_genus(selectedList[,1])
      for(z in glist){
        tmp <- data.frame(disease = dis, dataset= dataset, method = method, genus = z)
        selectdf <- rbind(selectdf, tmp)
      }
    }
  }
}
subselectdf <- selectdf[selectdf$disease == 'Diarrhea',]
methods <- unique(subselectdf$method)
datasets <- unique(subselectdf$dataset)
pldf <- data.frame()
for(m in methods){
  for(d1 in  datasets){
    for(d2 in datasets){
      if(d1 != d2){
        d1g <- unique(begdf$genus[begdf$dataset == d1])
        d2g <- unique(begdf$genus[begdf$dataset == d2])
        d1gs <- unique(subselectdf$genus[subselectdf$dataset == d1 & subselectdf$method == m])
        d2gs <- unique(subselectdf$genus[subselectdf$dataset == d2 & subselectdf$method == m])
        FS <- fisher(d1gs, d2gs, d1g, d2g)
        JC <- length(intersect(d1gs, d2gs))/length(union(d1gs, d2gs))
        #JC <- length(intersect(d1gs, d2gs))/max(length(unique(d1gs)), length(unique(d2gs)))
        tmp <- data.frame(method = m, number = length(intersect(d1gs, d2gs)), JC = JC, odds = as.numeric(FS$estimate), p = FS$p.value)
        pldf <- rbind(pldf, tmp)
      }
    }
  }
}

subselectdf <- selectdf[selectdf$disease == 'Obesity',]
methods <- unique(subselectdf$method)
datasets <- unique(subselectdf$dataset)
pldf1 <- data.frame()
for(m in methods){
  for(d1 in  datasets){
    for(d2 in datasets){
      if(d1 != d2){
        d1g <- unique(begdf$genus[begdf$dataset == d1])
        d2g <- unique(begdf$genus[begdf$dataset == d2])
        d1gs <- unique(subselectdf$genus[subselectdf$dataset == d1 & subselectdf$method == m])
        d2gs <- unique(subselectdf$genus[subselectdf$dataset == d2 & subselectdf$method == m])
        FS <- fisher(d1gs, d2gs, d1g, d2g)
        JC <- length(intersect(d1gs, d2gs))/length(union(d1gs, d2gs))
        #JC <- length(intersect(d1gs, d2gs))/max(length(unique(d1gs)), length(unique(d2gs)))
        tmp <- data.frame(method = m, number = length(intersect(d1gs, d2gs)), JC = JC, odds = as.numeric(FS$estimate), p = FS$p.value)
        pldf1 <- rbind(pldf1, tmp)
      }
    }
  }
}
p5 <- ggplot(pldf, aes(x = JC, y = odds, color = method, size = number)) + geom_point(alpha=0.5) +
  labs(x = 'Jaccard similarity', y = 'Odds ratio', color = 'Method', size = 'Overlapping genera') +
  scale_color_manual(values = mycolor) + ggtitle('Diarrhea') +
  scale_size(breaks = range(c(pldf$number, pldf1$number))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        legend.position = 'none')
p6 <- ggplot(pldf1, aes(x = JC, y = odds, color = method, size = number)) + geom_point(alpha=0.5) +
  labs(x = 'Jaccard similarity', y = 'Odds ratio', color = 'Method', size = 'Overlapping genera') +
  scale_color_manual(values = mycolor) + ggtitle('Obesity') + 
  scale_size_continuous(breaks = range(c(pldf$number, pldf1$number))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=10),
        legend.position = 'none')

p <- (p1 + p2) / (p3 + p4) + plot_annotation(tag_levels = 'A')
pdf("/home/yincheng23/ADlasso_manuscript/Figure/Figure4.pdf", width=10, height=8)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/Figure4.tiff", units="in", width=10, height=8, res=300)
p
dev.off()