library(ggplot2)
library(RColorBrewer)

setwd('path/to/data')

initial_slopes <- read.csv('flower_calibrated_initial_slopes.csv')
colnames(initial_slopes)[1] <- 'mouse'
initial_slopes$mouse <- factor(initial_slopes$mouse)
initial_slopes$day <- factor(initial_slopes$day)

svg('Figure_4A.svg', width=12, height=12)
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

svg('Figure_4B.svg', width=12, height=12)
ggplot(initial_slopes, aes(x=day, y=initial_slopes, group=mouse, colour=mouse)) +
  geom_smooth(method='lm', aes(color=mouse), se=FALSE, size=1.5) +
  geom_hline(yintercept=40, linetype='dashed', size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Linear Regression for Calibrated \nInitial Slopes Over Time',
       x='Day',
       y='Initial Slope (fm/s)',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))
dev.off()

ggplot(initial_slopes, aes(x=day, y=initial_slopes, group=mouse, colour=mouse)) +
  geom_smooth(method='lm', aes(color=mouse), se=TRUE, size=1.5) +
  geom_hline(yintercept=40, linetype='dashed', size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Linear Regression for Calibrated \nInitial Slopes Over Time',
       x='Day',
       y='Initial Slope (fm/s)',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))

days <- c(0, 14, 28, 42, 56)
means <- c()
for (i in days) {
  means <- c(means,
             mean(initial_slopes[which(initial_slopes$day==i),]$initial_slopes))
}

mean_initial_slopes <- data.frame(day=days,
                                  initial_slopes=means)

ggplot(initial_slopes, aes(x=day, y=initial_slopes)) +
  geom_boxplot(outlier.shape=NA) +
  geom_hline(yintercept=40, linetype='dashed', size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='FLOWER Mean Calibrated Initial Slopes Over Time',
       x='Day',
       y='Initial Slope (fm/s)',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))

ggplot(initial_slopes, aes(x=day, y=initial_slopes)) +
  geom_smooth(method='lm', se=TRUE, size=1.5) +
  geom_hline(yintercept=40, linetype='dashed', size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Linear Regression FLOWER \nMean Calibrated Initial Slopes Over Time',
       x='Day',
       y='Initial Slope (fm/s)',
       colour='Mouse') +
  scale_color_manual(values=brewer.pal(5, 'Dark2'))
