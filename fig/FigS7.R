library(Rmisc)
df1 <- read.csv('/home/yincheng23/ADlasso_manuscript/Script/mult_result.csv',stringsAsFactors = F)
df2 <- summarySE(df1, measurevar="AUC", groupvars=c("dataset","feature","class"))
df2$class[df2$class == 'Small adenoma'] <- 'Small\nadenoma'
df2$class[df2$class == 'Large adenoma'] <- 'Large\nadenoma'
df2$class[df2$class == 'High risk normal'] <- 'High risk\nnormal'
df2$class[df2$class == 'Adv. adenoma'] <- 'Advanced\nadenoma'
df2$class <- factor(df2$class, levels = c('Normal','High risk\nnormal','Adenoma','Advanced\nadenoma',
                                          'Small\nadenoma','Large\nadenoma','Cancer',
                                          'CD','UC'))

p <- ggplot(df2, aes(x = class, y = AUC, fill = feature)) +
     geom_bar(stat="identity", color=NA, position=position_dodge()) +
     geom_errorbar(aes(ymin = AUC - sd, ymax = AUC + sd),
                   position = position_dodge(0.9), width = .2) +
     facet_wrap(vars(dataset), scales = "free", nrow=2) +
     scale_fill_manual(values = c('#FFCC22','#0000FF')) + 
     labs(x='', y='AUC', fill='') + theme_bw(14) +
     theme(strip.background = element_blank(),
           text = element_text(size=14),
           axis.text.x = element_text(size = 10,colour = 'black'))

pdf("/home/yincheng23/ADlasso_manuscript/Figure/FigureS7.pdf", width=10, height=6)
p
dev.off()

tiff("/home/yincheng23/ADlasso_manuscript/Figure/FigureS7.tiff", units="in", width=10, height=6, res=300)
p
dev.off()