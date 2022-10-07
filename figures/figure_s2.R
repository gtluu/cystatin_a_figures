library(ggplot2)

ivis <- read.csv('900_mice_ivis.csv')
colnames(ivis)[1] <- 'Mouse'
ivis$Mouse <- factor(ivis$Mouse)
ivis$Day <- factor(ivis$Day)

svg('Figure_S2.svg', width=16, height=9)
ggplot(ivis, aes(x=Day, y=Avg_Rad_Eff, group=Mouse, colour=Mouse)) +
  geom_line(size=1.5) +
  geom_point(size=3) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Increase in Tumor Buden Over Time',
       x='Day',
       y='Average Radiant Efficiency',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))
dev.off()
