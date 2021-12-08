library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)
library(ggplot2)
library(RColorBrewer)

dir <- file.path('path/to/data')
spectra_list_raw <- import(dir)

preprocessing_workflow <- function(spectra_list, cores, snr, tol) {
  spectra_list <- trim(spectra_list, c(4000, 20000))
  spectra_list <- transformIntensity(spectra_list, method='sqrt', mc.cores=cores)
  spectra_list <- smoothIntensity(spectra_list, method='SavitzkyGolay', halfWindowSize=10, mc.cores=cores)
  spectra_list <- removeBaseline(spectra_list, method='TopHat', mc.cores=cores)
  spectra_list <- calibrateIntensity(spectra_list, method='TIC', mc.cores=cores)
  peaks <- detectPeaks(spectra_list, method='MAD', halfWindowSize=20, SNR=snr, mc.cores=cores)
  peaks <- alignPeaks(peaks, minFreq=0.75, tolerance=tol, SNR=snr, mc.cores=cores)
  
  featureMatrix <- intensityMatrix(peaks, spectra_list)
  filenames <- basename(unlist(lapply(peaks, function(x) attributes(x)$metaData$file)))
  attributes(featureMatrix)$dimnames[[1]] <- basename(filenames)
  
  featureDf <- as.data.frame(featureMatrix)
  
  drop_col_index <- which(sapply(seq(1, ncol(featureDf)), function(x) length(unique(featureDf[,x])) >= (nrow(featureDf) * 0.5)))
  featureDf <- featureDf[, drop_col_index]
  
  return(featureDf)
}

processed_spectra_list <- preprocessing_workflow(spectra_list_raw, 1, 3, 0.2)

# 1406 == indice for cystatin A
prot_hits_indices <- c(1406)

cystatin <- processed_spectra_list[prot_hits_indices[1]]

# split feature df row names
split_rownames <- str_split(rownames(cystatin), '_')
mouse_num_vector <- c()
day_vector <- c()
for (i in split_rownames) {
  mouse <- as.character(i[1])
  if (nchar(mouse) == 8) {
    mouse_num_vector <- c(mouse_num_vector, substr(mouse, 6, nchar(mouse)))
  } else if (nchar(mouse) == 4) {
    mouse_num_vector <- c(mouse_num_vector, substr(mouse, 2, nchar(mouse)))
  }
  if (tolower(i[2]) == tolower('Day')) {
    if (nchar(i[3]) == 2) {
      day_vector <- c(day_vector, as.numeric(i[3]))
    } else if (nchar(i[3] == 3)) {
      day_vector <- c(day_vector, as.numeric(substr(i[3], 1, 2)))
    }
  } else if (tolower(i[2]) == tolower('Week')) {
    if (nchar(i[3]) == 2) {
      day_vector <- c(day_vector, as.numeric(as.numeric(i[3]) * 7))
    } else if (nchar(i[3] == 3)) {
      day_vector <- c(day_vector, as.numeric(as.numeric(substr(i[3], 1, 2)) * 7))
    }
  }
}
rm(i, mouse)

# add m/z 11007, mouse, and day columns to feature df
cystatin$mz <- rep(as.numeric(colnames(cystatin)), nrow(cystatin))
cystatin$Mouse <- as.factor(mouse_num_vector)
cystatin$day <- as.factor(day_vector)
# change 11007.2 colname to intensity
colnames(cystatin)[1] <- 'intensity'

# remove outliers
q1 <- quantile(cystatin$intensity, .25)
q3 <- quantile(cystatin$intensity, .75)
iqr <- IQR(cystatin$intensity)

cystatin <- subset(cystatin, cystatin$intensity > (q1 - 1.5 * iqr) & cystatin$intensity < (q3 + 1.5 * iqr))

# get mean and median intensity per mouse per day
cystatin_mean <- data.frame(Mouse=factor(), day=factor(), intensity=numeric())
for (m in c(901, 902, 903, 904, 905)) {
  for (d in c(0, 7, 14, 21, 28, 35, 42, 56)) {
    df <- cystatin[which(cystatin$Mouse == m & cystatin$day == d),]
    cystatin_mean <- rbind(cystatin_mean, c(m, d, mean(df$intensity)))
  }
}
rm(m, d, df)
colnames(cystatin_mean) <- c('Mouse', 'day', 'mean_intensity')
cystatin_mean$Mouse <- as.factor(cystatin_mean$Mouse)
cystatin_mean$day <- as.factor(cystatin_mean$day)

svg('nat_comm_fig2a.svg', width=16, height=9)
ggplot(cystatin, aes(x=day, y=intensity, fill=Mouse)) +
  geom_boxplot(outlier.shape=NA) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Intensity of Cystatin A (m/z 11007 +/- 30 ppm) Over Time',
       x='Day',
       y='Processed Peak Intensity') +
  scale_fill_brewer(palette='Dark2')
dev.off()

svg('nat_comm_fig2b.svg', width=16, height=9)
ggplot(cystatin_mean, aes(x=day, y=mean_intensity, group=Mouse)) +
  geom_smooth(method='lm', aes(color=Mouse), se=FALSE, size=2.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mean Intensity of Cystatin A (m/z 11007 +/- 30 ppm) Over Time',
       x='Day',
       y='Processed Mean Peak Intensity') +
  scale_colour_brewer(palette='Dark2')
dev.off()
