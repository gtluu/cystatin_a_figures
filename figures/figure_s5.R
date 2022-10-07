library(ggplot2)

# standard concentration units == pM
# initial slope units == fm/s
calibration_curve <- read.csv('flower_cystatin_a_calibration_curve.csv')
colnames(calibration_curve) <- c('standard_concentration', 'initial_slope')
calibration_curve$standard_concentration <- log(calibration_curve$standard_concentration)

svg('Figure_S6.svg', width=12, height=12)
ggplot(calibration_curve, aes(x=standard_concentration, y=initial_slope)) +
  geom_point(size=3, colour='black') +
  geom_smooth(method='lm', se=FALSE, size=1.5, colour='black') +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Cystatin A Calibration Curve',
       x='log(Cystatin A Concentration) (pM)',
       y='Initial Slope (fm/s)')
dev.off()
