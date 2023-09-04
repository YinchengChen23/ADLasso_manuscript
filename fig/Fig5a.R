pvlTable <-  data.frame()
refpath <- '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results'

dataPath <- '/home/yincheng23/ADlasso_manuscript/Datasets/Studies/crc_zeller/'
data <- read.table(paste0(dataPath,'/ASV_table.txt'), sep = "\t", stringsAsFactors = F)
taxa <- read.table(paste0(dataPath,'/taxa_table.txt'), sep = "\t", stringsAsFactors = F)
meta <- read.table(paste0(dataPath,'/metadata.txt'), sep = "\t", stringsAsFactors = F)

ASV <- data
for(i in 1:ncol(ASV)){
  ASV[,i] <- ASV[,i]/ colSums(ASV)[i]
}

pdata <- data
pdata[pdata > 0] <- 1
idxNon <- rownames(meta)[meta$Class == 'Normal']
idxCan <- rownames(meta)[meta$Class == 'Cancer']

data_path = '/home/yincheng23/ADlasso_manuscript/Datasets/data/selected_results'
selected <- read.table(paste0(data_path,'/crc_zeller/ADlasso.txt'), sep = "\t", stringsAsFactors = F)
boxplot <- data.frame()
barplot <- data.frame()
wigplot <- data.frame()
xlabe <- c()
SigID <- c()
for(i in 1:nrow(selected)){
  tag <- selected$V1[i]
  tmp <- data.frame('group'= c(rep('Normal',length(idxNon)),rep('Cancer',length(idxCan))),
                    'ab' = as.numeric(c(ASV[tag,idxNon],ASV[tag,idxCan])))
  tmp$ASV <- tag
  tmp$Genus <- taxa$Genus[rownames(taxa) == tag]
  boxplot <- rbind(boxplot, tmp)
  tmp <- data.frame('group'= c('Normal','Cancer'),
                    'ab' = c(sum(pdata[tag,idxNon])/length(idxNon)*100, sum(pdata[tag,idxCan])/length(idxCan)*100))
  tmp$ASV <- tag
  tmp$Genus <- taxa$Genus[rownames(taxa) == tag]
  barplot <- rbind(barplot, tmp)
  wigplot <- rbind(wigplot,  data.frame('ASV'= tag, 'w' = selected$V2[i], 'Genus' = taxa$Genus[rownames(taxa) == tag]))
  
  xlabe <- c(xlabe, taxa$Genus[rownames(taxa) == tag])
  res <- wilcox.test(as.numeric(ASV[tag,idxNon]), as.numeric(ASV[tag,idxCan]))
  if(res$p.value < 0.05){
    SigID <- c(SigID, tag)
    print(paste(tag, taxa$Genus[rownames(taxa) == tag],res$p.value))
  }
}
keepID <- c()
for(x in unique(barplot$ASV)){
  if(unique(barplot$Genus[barplot$ASV == x]) == 'unidentified'){next}
  if(sum(barplot$ab[barplot$ASV == x]) > 140){
    keepID <- c(keepID, x)
  }
}
keepID <- c()
for(x in unique(barplot$Genus)){
  if(x == 'unidentified'){next}
  tagID <- barplot$ASV[barplot$Genus == x]
  t1 <- mean(barplot$ab[barplot$ASV %in% tagID & barplot$group == 'Normal'])
  t2 <- mean(barplot$ab[barplot$ASV %in% tagID & barplot$group == 'Cancer'])
  if(t1 > 50 & t2 > 50){
    keepID <- c(keepID, tagID)
  }
}

boxplot <- boxplot[boxplot$ASV %in% keepID,]
barplot <- barplot[barplot$ASV %in% keepID,]
wigplot <- wigplot[wigplot$ASV %in% keepID,]

wigplot$w <- wigplot$w * -1
wigplot$state <- 'neutral'
wm <- c()
for(x in unique(wigplot$Genus)){
  achor <- median(wigplot$w[wigplot$Genus == x])
  wm <- c(wm, achor)
  if(achor > 0){
    wigplot$state[wigplot$Genus == x] <- 'Cancer'
  }else if(achor < 0){
    wigplot$state[wigplot$Genus == x] <- 'Normal'
  }
}
names(wm) <- unique(wigplot$Genus)
boxplot$Genus <- factor(boxplot$Genus, levels = names(sort(wm)))
wigplot$Genus <- factor(wigplot$Genus, levels = names(sort(wm)))


taxalen <- length(unique(barplot$Genus))
mycolor <- c("Normal" = "#FFCC22", "Cancer" = "#EE7700")
p1 <- ggplot(boxplot, aes(x = Genus, y = ab, color = group)) +  coord_flip() +
  ylab("Relative abundance (%)") + scale_color_manual(values = mycolor) + 
  scale_x_discrete() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black', size=12),
        legend.title=element_blank(),
        legend.position = 'top')

for (i in 1:(taxalen-1)){
  p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = 10^-Inf, ymax = 10^Inf,
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
}
p1 <- p1 +  geom_boxplot(width = 0.5, outlier.shape = NA)

mycolor <- c("Normal" = "#FFCC22","Cancer" = "#EE7700")
p3 <- ggplot(wigplot, aes(x = Genus, y = w)) +  coord_flip() + 
  scale_color_manual(values = mycolor) +
  scale_fill_manual(values = mycolor) + 
  ylab("weight of ADlasso") + scale_x_discrete() +
  theme(panel.background = element_rect(fill = 'transparent'),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color='black'),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(colour='black',size=14),
        #axis.text.y = element_blank(),
        #axis.title.y = element_blank(),
        legend.title=element_blank(),
        legend.position = 'top')
for (i in 1:(taxalen-1)){
  p3 <- p3 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                      fill = ifelse(i %% 2 == 0, 'white', 'gray95'))
}
p3 <- p3 +  geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(width = 0.5,aes(fill = state))

pdf("/home/yincheng23/ADlasso_manuscript/Figure/Figure5_1.pdf", width=8, height=12)
p3
dev.off()


tiff("/home/yincheng23/ADlasso_manuscript/Figure/Figure5_1.tiff", units="in", width=8, height=12, res=300)
p3
dev.off()
