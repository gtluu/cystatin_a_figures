library(ggplot2)
library(RColorBrewer)

initial_slopes <- read.csv('flower_calibrated_initial_slopes.csv')
colnames(initial_slopes)[1] <- 'mouse'
initial_slopes$mouse <- factor(initial_slopes$mouse)
initial_slopes$day <- factor(initial_slopes$day)

svg('nat_comm_fig4a.svg', width=12, height=12)
ggplot(initial_slopes, aes(x=day, y=initial_slopes, group=mouse, colour=mouse)) +
  geom_line(size=1.5) +
  geom_point(size=3) +
  geom_hline(yintercept=40, linetype='dashed', size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='FLOWER Calibrated Initial Slopes Over Time',
       x='Day',
       y='Initial Slope (fm/s)',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))
dev.off()

svg('nat_comm_fig4b.svg', width=12, height=12)
ggplot(initial_slopes, aes(x=day, y=initial_slopes, group=mouse, colour=mouse)) +
  geom_smooth(method='lm', aes(color=mouse), se=FALSE, size=1.5) +
  geom_hline(yintercept=40, linetype='dashed', size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Linear Regression for Calibrated Initial Slopes Over Time',
       x='Day',
       y='Initial Slope (fm/s)',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))
dev.off()
