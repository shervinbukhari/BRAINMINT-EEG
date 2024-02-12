# Creates an overview based on summary files from EEG preprocessing for all participants split by study arm.


# Load packages and data
library(pastecs)
library(readr)
library(ggplot2)
library(scico)
library(ggtext)
library(tidyverse)


# BRAINMINT pipeline
#summary_files <- list.files("/tsd/p33/scratch/no-backup/projects/BRAINMINT/eeg/MMN/BRAINMINT_preprocessed/", 
 #                                 pattern = "*_summary.txt", full.names = TRUE)

# MINIMAL
summary_files <- list.files("/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/analysis_MMN/MMN_ERP_QC_output/preprocessed_data/pregnancy/", 
                           pattern = "*_summary.txt", full.names = TRUE)

#MINIMAL with line noise
#summary_files <- list.files("/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/analysis_MMN/MMN_ERP_QC_output/minimal_line_noise/", 
 #                          pattern = "*_summary.txt", full.names = TRUE)

# PREP Pipeline
#summary_files <- list.files("/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/analysis_MMN/MMN_ERP_QC_output/full_preprocessing/full_preprocessed_data/pregnancy/", 
 #                         pattern = "*_summary.txt", full.names = TRUE)


# Create one dataset
data <- data.frame()

for (i in 1:length(summary_files)) {
  filename <- basename(summary_files[i])  # Get the filename
  id <- strsplit(filename, "_")[[1]][1:2]  # Extract ID from filename
  id <- paste(id, collapse = "_")
  temp <- read.table(summary_files[i], header=TRUE, stringsAsFactors=FALSE)
#    read_delim(summary_files[i], "\t", escape_double = FALSE)
#  temp[1,2] <- as.character(temp[1,2])
  temp[1, 1] <- id  
  data <- rbind(data, temp)
  rm(temp)
}

data$group <- NA
data$group[data$ID %in% 65000:66999] <- "adolescence"
data$group[data$ID %in% 68000:68999] <- "pregnancy"

data$group <- ifelse(substr(data$ID, 1, 2) %in% c("65","66"), "adolescence", 
                   ifelse(substr(data$ID, 1, 2) == "68", "pregnancy", NA))
data$group <- as.factor(data$group)

data = data %>%
  filter(group == 'pregnancy')

data$Num_InterpolatedChan <- as.numeric(data$Num_InterpolatedChan)
data$NumComp <- as.numeric(data$NumComp)
data$RejComp <- as.numeric(data$RejComp)
data$PctTimeRetained.ASR. <- as.numeric(data$PctTimeRetained.ASR.)*100
data$PctRejectedEpochs <- as.numeric(data$PctRejectedEpochs)*100
data$PctRejComp <- data$RejComp/data$NumComp*100



# Summary stats
summary(data)
stat.desc(data)


# Plots
# as histogram
hist_rejComp <- ggplot(data, aes(RejComp), fill = group) +
  geom_histogram(binwidth = 1, aes(fill = group), alpha = 0.4) +
  ggtitle("Rejected Components") +
  scale_fill_manual(values = c("dodgerblue3", "forestgreen")) +
  xlab("Number of Rejected Components") +
  ylab("Count") +
  scale_color_manual(values = c("dodgerblue3","forestgreen")) +
  theme_minimal()

hist_prosRejComp <- ggplot(data, aes(PctRejComp), fill= group) +
  geom_histogram(binwidth = 1, aes(fill = group), alpha = 0.4) +
  ggtitle("Pct. Rejected Components") +
  xlab("% Rejected Components") +
  ylab("Count") +
  scale_color_manual(values = c("dodgerblue3","forestgreen")) +
  scale_fill_manual(values = c("dodgerblue3", "forestgreen")) +
  theme_minimal()

hist_prosTimeRetained <-ggplot(data, aes(PctTimeRetained.ASR.), fill= group) +
  geom_histogram(binwidth = 1, aes(fill = group), alpha = 0.4) +
  ggtitle("Data Retained after ASR") +
  xlab("Pct. Time Retained") +
  ylab("Count") +
  scale_color_manual(values = c("dodgerblue3","forestgreen")) +
  scale_fill_manual(values = c("dodgerblue3", "forestgreen")) +
  theme_minimal()

hist_prosEpochRejected <-ggplot(data, aes(PctRejectedEpochs), fill= group) +
  geom_histogram(binwidth = 1, aes(fill = group), alpha = 0.4) +
  ggtitle("Percentage of Epochs Rejected") +
  xlab("Pct. Rejected Epochs ") +
  ylab("Count") +
  scale_color_manual(values = c("dodgerblue3","forestgreen")) +
  scale_fill_manual(values = c("dodgerblue3", "forestgreen")) +
  theme_minimal()


# hist_rejEpochs <- ggplot(data, aes(PctTimeRetained), fill= group) +
#   geom_histogram(binwidth = 1, aes(fill = group), alpha = 0.4) +
#   ggtitle("Rejected Epochs") +
#   xlab("Number of Epochs") +
#   ylab("Count") +
#   scale_color_manual(values = c("dodgerblue3","forestgreen")) +
#   scale_fill_manual(values = c("dodgerblue3", "forestgreen")) +
#   theme_minimal()

hist_intChan <- ggplot(data, aes(Num_InterpolatedChan), fill= group) +
  geom_histogram(binwidth = 1, aes(fill = group), alpha = 0.4) +
  ggtitle("Interpolated Channels") +
  xlab("Number of Channels") +
  ylab("Count") +
  scale_color_manual(values = c("dodgerblue3","forestgreen")) +
  scale_fill_manual(values = c("dodgerblue3", "forestgreen")) +
  theme_minimal()

# save plots
for (plot in ls(pattern = "hist_")) {
  ggsave(paste(plot, ".pdf", sep = ""), get(plot), width = 20, heigh = 12.5, bg = "white")
}
