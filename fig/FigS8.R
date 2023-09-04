library(ggplot2)
library(ggridges)
library(ggtext)
library(glue)

pvlTable <-  data.frame()
dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/'
dataSets <- dir(dataPath)
for(dataSet in dataSets){
  #----------- read raw data --------------
  if(dataSet == 'LIBSVM_rs'){next}
  setwd(dataPath)
  setwd(dataSet)
  data <- read.table('ASV_table.txt', sep = "\t", stringsAsFactors = F)
  data[data > 0] <- 1
  prevalence <- rowSums(data)/ncol(data)*100
  pvlTable <- rbind(pvlTable, data.frame('Data'=dataSet,'pvl'=prevalence))
}
pvlTable$Data <- factor(pvlTable$Data, levels = sort(unique(pvlTable$Data)))

highlight = function(x, pat, color="black", family="") {
  ifelse(grepl(pat, x), glue("<b style='font-family:{family}; color:{color}'>{x}</b>"), x)
}

p1 <- ggplot(pvlTable, aes(x = pvl, y = Data, fill = stat(x))) +
             geom_density_ridges_gradient(scale = 3, size = 0.3, rel_min_height = 0.01) +
             scale_fill_viridis_c(name = "", option = "C", limits = c(0, 100),
                                  trans = scales::pseudo_log_trans(sigma = 0.1)) +
             scale_x_continuous(trans='log10') + theme_minimal() +
             theme(legend.position = 'none') + xlab('Prevalence (%)') +
             scale_y_discrete(labels= function(x) highlight(x, "ob_zupancic", "red")) +
             theme(axis.text.y=element_markdown())


PCdf <- read.table("/home/yincheng23/ADlasso_manuscript/Script/performance_statistics.txt", sep = "\t", stringsAsFactors = F, header = T, row.names = 1)
PCdf <- PCdf[PCdf$Matrix == "AUC", ]
PCdf <- PCdf[PCdf$Method == "ADlasso", ]
PCdf <- rbind(PCdf, data.frame('Data'='block', 'Method'='ADlasso', 'Matrix'='AUC', 'value'=NA))
PCdf$Data <- factor(PCdf$Data, levels = c(levels(pvlTable$Data),'block'))
PCdf$AUC <- ifelse(PCdf$value > 0.99, '> 0.99', '< 0.99')

p2 <- ggplot(PCdf, aes(x = Data, y = value, fill = AUC)) +
             geom_bar(stat="identity") + theme_minimal() + coord_flip() + ylab('AUC') +
             scale_fill_manual(values = c('#FFAA33','#AAAAAA')) +
             theme(axis.text.y = element_blank(),
                   axis.title.y = element_blank())

p <- p1 + p2 + plot_layout(widths= c(7,3))


pdf("/home/yincheng23/ADlasso_manuscript/Figure/FigureS8.pdf", width=8, height=8)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/FigureS8.tiff", units="in", width=8, height=8, res=300)
p
dev.off()