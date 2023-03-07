library(MALDIquant)
library(MALDIquantForeign)
library(MALDIrppa)
library(stringr)
library(ggplot2)
library(RColorBrewer)

# set working dir and import data
dir <- file.path('path/to/data')  # data from MSV000090494
spectra_list_raw <- import(dir)

# wrapper for preprocessing steps that outputs dataframe with feature matrix
preprocessing_workflow <- function(spectra_list, cores, snr, tol) {
  spectra_list <- trim(spectra_list, c(4000, 20000))
  spectra_list <- transformIntensity(spectra_list, method='sqrt', mc.cores=cores)
  spectra_list <- smoothIntensity(spectra_list, method='SavitzkyGolay', halfWindowSize=10, mc.cores=cores)
  spectra_list <- removeBaseline(spectra_list, method='TopHat', mc.cores=cores)
  spectra_list <- calibrateIntensity(spectra_list, method='TIC', mc.cores=cores)
  peaks <- detectPeaks(spectra_list, method='MAD', halfWindowSize=20, SNR=snr, mc.cores=cores)
  peaks <- alignPeaks(peaks, minFreq=0.25, tolerance=tol, SNR=snr, mc.cores=cores)
  
  featureMatrix <- intensityMatrix(peaks, spectra_list)
  filenames <- basename(unlist(lapply(peaks, function(x) attributes(x)$metaData$file)))
  attributes(featureMatrix)$dimnames[[1]] <- basename(filenames)
  
  featureDf <- as.data.frame(featureMatrix)
  
  drop_col_index <- which(sapply(seq(1, ncol(featureDf)), function(x) length(unique(featureDf[,x])) >= (nrow(featureDf) * 0.5)))
  featureDf <- featureDf[, drop_col_index]
  
  return(featureDf)
}

# processed data to feature matrix df
processed_spectra_list <- preprocessing_workflow(spectra_list_raw, 1, 3, 0.2)

# NOTE: mouse cystatin A was 11007.2
mz <- as.numeric(colnames(processed_spectra_list))
mz[1050:1100]
mz[1099]  # 11006.6 +/- 55 ppm
prot_hits_indices <- c(1099)

# get df with intensities for cystatin signal
cystatin <- processed_spectra_list[prot_hits_indices[1]]

# split feature df row names
split_rownames <- str_split(rownames(cystatin), '_')
tampon_vector <- c()
for (i in split_rownames) {
  tampon_vector <- c(tampon_vector, as.character(i[1]))
}
rm(i)

# modify feature df
cystatin$mz <- rep(as.numeric(colnames(cystatin)), nrow(cystatin))
cystatin$Tampon <- as.factor(tampon_vector)
# change 11006 colname to intensity
colnames(cystatin)[1] <- 'intensity'

# read in annotation data
# add condition via annotations
annotations <- read.csv('tampon_annotations.csv', header=TRUE, colClasses='character')
colnames(annotations) <- c('id', 'date_consented', 'diagnosis', 'full_diagnosis', 'comments')
weights <- read.csv('tampon_weights_2020NOV.csv', header=TRUE, colClasses='character', nrows=54)
colnames(weights) <- c('sample', 'id', 'initial_wt', 'final_wt', 'wt')

simple_annot <- merge(weights, annotations, by='id', all=TRUE)
simple_annot <- simple_annot[,c('id', 'sample', 'diagnosis')]
simple_annot[is.na(simple_annot)] <- 'Unknown'
simple_annot[3,]$diagnosis <- 'Unknown'
simple_annot <- simple_annot[order(simple_annot$sample),]
row.names(simple_annot) <- NULL
simple_annot <- simple_annot[which(simple_annot$sample %in% c('TP10', 'TP24', 'TP34', 'TP42', 'TP48', 'TP50')),]

conditions <- simple_annot$diagnosis

cystatin <- cbind(cystatin, condition=c(rep(simple_annot$diagnosis[1:6], each=48)))

# remove outliers
q1 <- quantile(cystatin$intensity, .25)
q3 <- quantile(cystatin$intensity, .75)
iqr <- IQR(cystatin$intensity)

cystatin <- subset(cystatin, cystatin$intensity > (q1 - 1.5 * iqr) & cystatin$intensity < (q3 + 1.5 * iqr))

# get mean and median intensity per mouse per day
cystatin_mean_1 <- data.frame(Tampon=factor(), Condition=factor(), mean_intensity=numeric())
for (t in c('TP10', 'TP24', 'TP34', 'TP42', 'TP48', 'TP50')) {
  df <- cystatin[which(cystatin$Tampon == t),]
  cystatin_mean_1 <- rbind(cystatin_mean_1, c(t, unique(df$condition), mean(df$intensity)))
}
colnames(cystatin_mean_1) <- c('tampon', 'condition', 'intensity')
cystatin_mean_1$tampon <- as.factor(cystatin_mean_1$tampon)
cystatin_mean_1$condition <- as.factor(cystatin_mean_1$condition)
cystatin_mean_1$intensity <- as.numeric(cystatin_mean_1$intensity)

# benign vs ovarian cancer t-test
t.test(intensity ~ condition, data=cystatin)

# bar plot with MALDI intensities
tampon_subset <- cystatin_mean_1[which(cystatin_mean_1$tampon %in% c('TP10', 'TP24', 'TP34', 'TP42', 'TP48', 'TP50')),]
colnames(tampon_subset)[1] <- 'Tampon'
svg('Figure_5A.svg', width=16, height=9)
ggplot(cystatin, aes(x=condition, y=intensity, fill=Tampon)) +
  geom_boxplot(outlier.shape=NA) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Mean Intensity of Cystatin A (m/z 11007 +/- 55 ppm)',
       x='Condition',
       y='Processed Peak Intensity')
dev.off()

data <- read.csv('flower_calibrated_initial_slopes_tampons_long.csv')
colnames(data) <- c('time', 'shift', 'Tampon')

svg('Figure_5B.svg', width=16, height=9)
ggplot(data, aes(x=time)) +
  geom_line(aes(y=shift, colour=Tampon), size=1.5) +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.background=element_blank(),
        axis.line=element_line(colour='black', size=1.5),
        text=element_text(size=32)) +
  labs(title='Patient Tampon Extracts - FLOWER Relative Shifts Over Time',
       x='Time (sec)',
       y='Relative Shift (fm)',
       colour='Tampon') +
  scale_fill_brewer(palette='Paired')
dev.off()
