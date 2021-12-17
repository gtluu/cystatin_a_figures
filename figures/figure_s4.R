library(ggplot2)
library(RColorBrewer)
library(stringr)

setwd('F:/Backups/F Drive NVME 12-2021 2/melissa/r')

data <- read.csv('Figure_S4.csv')
colnames(data)[1] <- 'time'
data$Sample <- str_replace(data$Sample, 'buffer', 'Buffer Subtraction')
data$Sample <- str_replace(data$Sample, '100pM', '100 pM Cystatin A')
data$Sample <- str_replace(data$Sample, '10nM', '10 nM Cystatin A')
data$Sample <- str_replace(data$Sample, '100nM', '100 nM Cystatin A')
data$Sample <- str_replace(data$Sample, '200nM', '200 nM Cystatin A')
data$Sample <- factor(data$Sample, levels=c('Buffer Subtraction',
                                            '100 pM Cystatin A',
                                            '10 nM Cystatin A',
                                            '100 nM Cystatin A',
                                            '200 nM Cystatin A'))

svg('nat_comm_figs4.svg', width=16, height=9)
ggplot(data, aes(x=time, y=wavelength, colour=Sample)) +
  geom_line(size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='FLOWER Relative Shifts of Cystatin A Over Time',
       x='Time (sec)',
       y='Wavelength (fm)',
       colour='Concentration') +
  scale_fill_brewer(palette='Paired')
dev.off()
