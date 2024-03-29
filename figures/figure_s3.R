library(stringr)
library(ggplot2)
library(RColorBrewer)

m901 <- read.csv('M901_updated.csv')
colnames(m901)[1] <- 'time'
m901$sample <- str_replace(m901$sample, '100pM', '100 pM Cystatin A')
m901$sample <- str_replace(m901$sample, 'day0', 'Day 0')
m901$sample <- str_replace(m901$sample, 'day14', 'Day 14')
m901$sample <- str_replace(m901$sample, 'day28', 'Day 28')
m901$sample <- str_replace(m901$sample, 'day42', 'Day 42')

m902 <- read.csv('M902_updated.csv')
colnames(m902)[1] <- 'time'
m902$sample <- str_replace(m902$sample, '100pM', '100 pM Cystatin A')
m902$sample <- str_replace(m902$sample, 'day0', 'Day 0')
m902$sample <- str_replace(m902$sample, 'day14', 'Day 14')
m902$sample <- str_replace(m902$sample, 'day28', 'Day 28')
m902$sample <- str_replace(m902$sample, 'day42', 'Day 42')

m903 <- read.csv('M903_updated.csv')
colnames(m903)[1] <- 'time'
m903$sample <- str_replace(m903$sample, '100pM', '100 pM Cystatin A')
m903$sample <- str_replace(m903$sample, 'day0', 'Day 0')
m903$sample <- str_replace(m903$sample, 'day14', 'Day 14')
m903$sample <- str_replace(m903$sample, 'day28', 'Day 28')
m903$sample <- str_replace(m903$sample, 'day42', 'Day 42')
m903$sample <- str_replace(m903$sample, 'day56', 'Day 56')

m904 <- read.csv('M904_updated.csv')
colnames(m904)[1] <- 'time'
m904$sample <- str_replace(m904$sample, '100pM', '100 pM Cystatin A')
m904$sample <- str_replace(m904$sample, 'day0', 'Day 0')
m904$sample <- str_replace(m904$sample, 'day14', 'Day 14')
m904$sample <- str_replace(m904$sample, 'day28', 'Day 28')
m904$sample <- str_replace(m904$sample, 'day42', 'Day 42')

m905 <- read.csv('M905_updated.csv')
colnames(m905)[1] <- 'time'
m905$sample <- str_replace(m905$sample, '100pM', '100 pM Cystatin A')
m905$sample <- str_replace(m905$sample, 'day0', 'Day 0')
m905$sample <- str_replace(m905$sample, 'day14', 'Day 14')
m905$sample <- str_replace(m905$sample, 'day28', 'Day 28')
m905$sample <- str_replace(m905$sample, 'day42', 'Day 42')
m905$sample <- str_replace(m905$sample, 'day56', 'Day 56')

svg('Figure_S3A.svg', width=16, height=9)
ggplot(m901, aes(x=time)) +
  geom_line(aes(y=shift, colour=sample), size=1.5) +
  geom_line(aes(y=fitting, colour=sample), size=1.5, linetype='dashed') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mouse 901 FLOWER Relative Shifts Over Time',
       x='Time (sec)',
       y='Relative Shift (fm)',
       colour='Mouse/Sample') +
  scale_fill_brewer(palette='Paired')
dev.off()

svg('Figure_S3B.svg', width=16, height=9)
ggplot(m902, aes(x=time)) +
  geom_line(aes(y=shift, colour=sample), size=1.5) +
  geom_line(aes(y=fitting, colour=sample), size=1.5, linetype='dashed') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mouse 902 FLOWER Relative Shifts Over Time',
       x='Time (sec)',
       y='Relative Shift (fm)',
       colour='Mouse/Sample') +
  scale_fill_brewer(palette='Paired')
dev.off()

svg('Figure_S3C.svg', width=16, height=9)
ggplot(m903, aes(x=time)) +
  geom_line(aes(y=shift, colour=sample), size=1.5) +
  geom_line(aes(y=fitting, colour=sample), size=1.5, linetype='dashed') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mouse 903 FLOWER Relative Shifts Over Time',
       x='Time (sec)',
       y='Relative Shift (fm)',
       colour='Mouse/Sample') +
  scale_fill_brewer(palette='Paired')
dev.off()

svg('Figure_S3D.svg', width=16, height=9)
ggplot(m904, aes(x=time)) +
  geom_line(aes(y=shift, colour=sample), size=1.5) +
  geom_line(aes(y=fitting, colour=sample), size=1.5, linetype='dashed') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mouse 904 FLOWER Relative Shifts Over Time',
       x='Time (sec)',
       y='Relative Shift (fm)',
       colour='Mouse/Sample') +
  scale_fill_brewer(palette='Paired')
dev.off()

svg('Figure_S3E.svg', width=16, height=9)
ggplot(m905, aes(x=time)) +
  geom_line(aes(y=shift, colour=sample), size=1.5) +
  geom_line(aes(y=fitting, colour=sample), size=1.5, linetype='dashed') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mouse 905 FLOWER Relative Shifts Over Time',
       x='Time (sec)',
       y='Relative Shift (fm)',
       colour='Mouse/Sample') +
  scale_fill_brewer(palette='Paired')
dev.off()
