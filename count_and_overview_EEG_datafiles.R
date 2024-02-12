# Problem: We want to have an easily accessible overview of EEG data that shows 
# 1. which participants we have raw data files for 2. which files are available for each participant

#Solution: This BRAINMINT script can be used to 
#A) count the number of (correctly named) raw EEG files and behavioral files available on cluster for each task.
#B) Compare available data to what has been added to the logs
#C) How many participants have a folder and finally outputs a .xlsx table where each participant is listed and which EEG and behavioral files are available for them. 

# Adapted for linux and slurm compatibility from count_EEG_data.Rmd on windows by Shervin Bukhari 27.09.2022

# !adapt adolescence code for multiple timepoints if needed! #change line 53
## COUNTS DATA ON CLUSTER AND OUTPUTS TABLES TO CLUSTER DIRECTORY
## 
#install.packages("openxlsx","tidyverse")
library(openxlsx)
library(tidyverse)
library(readxl)
#############
#ADOLESCENCE#
#############
# Read Datacollection logg and remove unnecessary columns
adolescence_logg<-  read_xlsx(path ="/tsd/p33/data/durable/groups/imaging/BRAINMINT/project_organisation/data_collection/Ungdom_logg061121.xlsm", col_names = TRUE, range = "EEG!C11:H999")#load data collection logg
# Remove empty rows and unnecessary columns
#adolescence_logg$`Participant(ID)` <- paste(adolescence_logg$`Participant(ID)`,adolescence_logg$`Time point`,sep="_")

adolescence_logg <- adolescence_logg[!(is.na(adolescence_logg$Date)),c("Date","Participant(ID)")]
names(adolescence_logg) <- c("logg_date","logg_ID")
adolescence_logg <- adolescence_logg[order(adolescence_logg$logg_ID),] 
adolescence_logg$logg_date <- format(as.Date(adolescence_logg$logg_date, format = "%d/%m/%Y"), "%Y-%m-%d")

#Load EEG data, listing the folders for each participant etc.
eegfolder_adolesc <- "/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_adolescents/"
eegdirs_adolesc <- list.dirs(path = eegfolder_adolesc, full.names = TRUE, recursive = FALSE) #list of all directories. 1 per participant
eegdirs_adolesc<-eegdirs_adolesc[ grepl("65*|66*", eegdirs_adolesc)]
length(eegdirs_adolesc) #number of directories, or participants
eegfiles_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="*.bdf") #list all files 
length(eegfiles_adolesc) #total number of files
eegrest_closed_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="REST_closed", ignore.case = TRUE) #list all files containing task name for example"REST_closed", case insensitive
eegrest_open_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="REST_open", ignore.case =TRUE)
eegRLWM_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="RLWM", ignore.case = TRUE)
eegMMN_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="MMN", ignore.case = TRUE)
eegSST_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="SST", ignore.case = TRUE)
eegEmopics_adolesc <- list.files(eegdirs_adolesc, all.files=TRUE, pattern="Emopics", ignore.case = TRUE)

#Load 3D scans
folder_adolesc_3D <- "/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_head_scans/3Dscans_adolescents/"
files_adolesc_3D <- list.files(folder_adolesc_3D, pattern="*.ply", all.files=TRUE)

#Get Subject ID from EEG folder names
IDs_adolesc<-gsub("/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_adolescents//", "", eegdirs_adolesc) #creates a vector with subject ID names from each folder in /egg_raw_adolescents/
Timepoint_adolesc <- sapply(strsplit(IDs_adolesc, split='_', fixed=TRUE), function(x) (x[2])) #Not used but can be implemented if there are multiple timepoints
IDs_adolesc<-str_remove(IDs_adolesc,"_01")
eegTasks<-c("Rest_closed","Rest_open","RLWM","MMN","SST","Emopics","3D_Scan") #EEGTasks # In addition to one EEG file per task and one 3D scan, there should be one file with behavioral data for RLWM, SST and emopics respectively
EEGdata_adolesc<- data.frame(matrix(nrow=length(IDs_adolesc),ncol=length(eegTasks))) #creates data frame the length of number of participants x (tasks + 3D Scan column)
dimnames(EEGdata_adolesc)<-list(IDs_adolesc,eegTasks) #name the columns after task and rows after IDs

#For each participant with a file, if there is a file matching .bdf-file for each of the 6 tasks, put a Y in the data frame, otherwise a N
for (i in seq_along(eegdirs_adolesc)) 
{
  if (0 != length(list.files(eegdirs_adolesc[i], all.files=TRUE, pattern="REST_closed",ignore.case = TRUE))){
    EEGdata_adolesc[i,1]<-"Yes"
  } else {
    EEGdata_adolesc[i,1]<-"X"
  }
}

for (i in seq_along(eegdirs_adolesc)) 
{
  if (0!=length(list.files(eegdirs_adolesc[i], all.files=TRUE, pattern="REST_open",ignore.case = TRUE))){
    EEGdata_adolesc[i,2]<-"Yes"
  } else {
    EEGdata_adolesc[i,2]<-"X"
  }
}

for (i in seq_along(eegdirs_adolesc)) 
{
  if (0!=length(list.files(eegdirs_adolesc[i], all.files=TRUE, pattern="RLWM",ignore.case = TRUE))){
    EEGdata_adolesc[i,3]<-"Yes"
  } else {
    EEGdata_adolesc[i,3]<-"X"
  }
}
for (i in seq_along(eegdirs_adolesc)) 
{
  if (0!=length(list.files(eegdirs_adolesc[i], all.files=TRUE, pattern="MMN",ignore.case = TRUE))){
    EEGdata_adolesc[i,4]<-"Yes"
  } else {
    EEGdata_adolesc[i,4]<-"X"
  }
}
for (i in seq_along(eegdirs_adolesc)) 
{
  if (0!=length(list.files(eegdirs_adolesc[i], all.files=TRUE, pattern="SST",ignore.case = TRUE))){
    EEGdata_adolesc[i,5]<-"Yes"
  } else {
    EEGdata_adolesc[i,5]<-"X"
  }
}
for (i in seq_along(eegdirs_adolesc)) 
{
  if (0!=length(list.files(eegdirs_adolesc[i], all.files=TRUE, pattern="Emopics",ignore.case = TRUE))){
    EEGdata_adolesc[i,6]<-"Yes"
  } else {
    EEGdata_adolesc[i,6]<-"X"
  }
}

# See if there is a file in files_adolesc_3D matching the participant id, insert Y if TRUE or N if false
EEGdata_adolesc <- EEGdata_adolesc %>% 
  rownames_to_column("id") %>% 
  rowwise() %>% 
  mutate(`3D_Scan` = ifelse(any(str_detect(files_adolesc_3D, id)), "Yes", "X"))

### NOW WE WANT TO REPEAT A SIMILAR APPROACH BUT FOR BEHAVIORAL DATA FILES, ADDING THEM TO THE DATAFRAME ## 
behaviorfolder_adolesc <-"/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_behavior/eeg_behavior_adolescents/"
behaviordir_adolesc <- list.dirs(path = behaviorfolder_adolesc, full.names = TRUE, recursive = FALSE) #list of all directories. 1 per participant
behaviordir_adolesc<-behaviordir_adolesc[ grepl("65*|66*", behaviordir_adolesc)]
behaviorfiles_adolesc <- list.files(behaviorfolder_adolesc, all.files=TRUE, pattern="*.csv|*.mat") #list all files 
length(behaviorfiles_adolesc) #total number of files

behaviorRLWM_adolesc <- list.files(behaviordir_adolesc, all.files=TRUE, pattern="rlwmpst",ignore.case = TRUE)
behaviorSST_adolesc <- list.files(behaviordir_adolesc, all.files=TRUE, pattern="SST", ignore.case = TRUE)
behaviorEmopics_adolesc <- list.files(behaviordir_adolesc, all.files=TRUE, pattern="Emopics", ignore.case = TRUE)

#See if there is a file in behaviorRLWM, behaviorSST and behaviorEmopics matching the participant id, insert Y if TRUE or N if false
EEGdata_adolesc <- EEGdata_adolesc %>% 
  #  rownames_to_column("id") %>% 
  rowwise() %>% 
  mutate(`behaviorRLWM` = ifelse(any(str_detect(behaviorRLWM_adolesc, id)), "Yes", "X"))
EEGdata_adolesc <- EEGdata_adolesc %>% 
  #  rownames_to_column("id") %>% 
  rowwise() %>% 
  mutate(`behaviorSST` = ifelse(any(str_detect(behaviorSST_adolesc, id)), "Yes", "X"))
EEGdata_adolesc <- EEGdata_adolesc %>% 
  # rownames_to_column("id") %>% 
  rowwise() %>% 
  mutate(`behaviorEmopics` = ifelse(any(str_detect(behaviorEmopics_adolesc, id)), "Yes", "X"))

EEGdata_adolesc <- merge(adolescence_logg,EEGdata_adolesc, by.x = "logg_ID", by.y = "id", all.x = TRUE)
EEGdata_adolesc <- EEGdata_adolesc %>% arrange(desc(logg_date))
#Write the table as a xlsx file with current date added
currentDate <- Sys.Date()
adol_filename <- paste("/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/EEG_DataLists/Adolescence_EEGdatalist",currentDate,".xlsx",sep="")
write.xlsx(EEGdata_adolesc,file=adol_filename, rowNames=TRUE)

#############
#PREGNANCY###
#############
##Since there are several sessions per participant, we will need to make some adjustments so that different time points are accounted for. 
rm(list = ls()) 

pregnancy_logg<-  read_xlsx(path = "/tsd/p33/data/durable/groups/imaging/BRAINMINT/project_organisation/data_collection/Gravid_logg020721.xlsm", col_names = TRUE, range = "EEG!C7:H999")#load data collection logg
# Remove empty rows and unnecessary columns
pregnancy_logg$`Participant(ID)` <- paste(pregnancy_logg$`Participant(ID)`,pregnancy_logg$`Time point`,sep="_")
pregnancy_logg <- pregnancy_logg[!(is.na(pregnancy_logg$Date)),c("Date","Participant(ID)")] #remove empty cells
names(pregnancy_logg) <- c("logg_date","logg_ID")
pregnancy_logg <- pregnancy_logg[order(pregnancy_logg$logg_ID),]
pregnancy_logg$logg_date <- format(as.Date(pregnancy_logg$logg_date, format = "%d/%m/%Y"), "%Y-%m-%d")

#Read folders and files
eegfolder_preg <- "/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_pregnancy/"
eegdirs_preg = list.dirs(path = eegfolder_preg, full.names = TRUE, recursive = FALSE) #list of all directories. 1 per participant
eegdirs_preg<-eegdirs_preg[ grepl("68*", eegdirs_preg)] #removes folders not correctly named.
length(eegdirs_preg) #number of directories, or participants
eegfiles_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="*.bdf") #list all files 
length(eegfiles_preg) #total number of files
eegrest_closed_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="REST_closed", ignore.case = TRUE) #list all files containing task name for example"REST_closed", case insensitive
eegrest_open_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="REST_open", ignore.case =TRUE)
eegRLWM_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="RLWM", ignore.case = TRUE)
eegMMN_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="MMN", ignore.case = TRUE)
eegSST_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="SST", ignore.case = TRUE)
eegEmopics_preg <- list.files(eegdirs_preg, all.files=TRUE, pattern="Emopics", ignore.case = TRUE)

folder_preg_3D <- '/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_head_scans/3Dscans_pregnancy'
files_preg_3D <- list.files(folder_preg_3D, pattern="*.ply", all.files=TRUE)

# Because behavioral files do not follow the same naming convention as the manually named files above, we need to convert the time point names to their corresponding integer and the ID without their timepoint suffix 
IDs_preg<-gsub("/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_pregnancy//", "", eegdirs_preg) #creates a vector with subject ID names from each folder
Timepoint_preg <- sapply(strsplit(IDs_preg, split='_', fixed=TRUE), function(x) (x[2]))
Timepoint_integer_preg <- Timepoint_preg
Timepoint_integer_preg <- gsub(pattern="pre", replacement="1", Timepoint_integer_preg) 
Timepoint_integer_preg <- gsub(pattern="mid2", replacement="2", Timepoint_integer_preg)
Timepoint_integer_preg <- gsub(pattern="post1", replacement="3",Timepoint_integer_preg) 
Timepoint_integer_preg <- gsub(pattern="post2", replacement="4",Timepoint_integer_preg)
Timepoint_integer_preg <- gsub(pattern="control18m", replacement="5",Timepoint_integer_preg)
Timepoint_integer_preg<-data.frame(Timepoint_integer_preg)
eegTasks<-c("stripped_ID","Timepoint","Timepoint_integer","ID_Timepoint", "Rest_closed","Rest_open","RLWM","MMN","SST","Emopics","3D_Scan","behaviorRLWM","behaviorSST","behaviorEmopics") #Tasks + 3D scan

IDs_preg_stripped<-sapply(strsplit(IDs_preg, split='_', fixed=TRUE), function(x) (x[1]))
IDs_preg_stripped<-data.frame(IDs_preg_stripped)

EEGdata_preg<- data.frame(matrix(nrow=length(IDs_preg),ncol=length(eegTasks))) #creates data frame the length of number of participants x tasks
dimnames(EEGdata_preg)<-list(IDs_preg,eegTasks) #name the columns after task and rows after IDs
EEGdata_preg[,1]<-IDs_preg_stripped
EEGdata_preg[,2]<-Timepoint_preg
EEGdata_preg[,3]<-Timepoint_integer_preg
EEGdata_preg <- EEGdata_preg %>% 
  mutate(Timepoint_integer = as.numeric(Timepoint_integer),
         ID_Timepoint = str_glue("{stripped_ID}_{Timepoint_integer}"))

#For each participant with a file, if there is a file matching the task name, put a Y in the data frame, otherwise a N
for (i in seq_along(eegdirs_preg)) 
{
  if (0 != length(list.files(eegdirs_preg[i], all.files=TRUE, pattern="REST_closed",ignore.case = TRUE))){
    EEGdata_preg[i,5]<-"Yes"
  } else {
    EEGdata_preg[i,5]<-"X"
  }
}

for (i in seq_along(eegdirs_preg)) 
{
  if (0!=length(list.files(eegdirs_preg[i], all.files=TRUE, pattern="REST_open",ignore.case = TRUE))){
    EEGdata_preg[i,6]<-"Yes"
  } else {
    EEGdata_preg[i,6]<-"X"
  }
}

for (i in seq_along(eegdirs_preg)) 
{
  if (0!=length(list.files(eegdirs_preg[i], all.files=TRUE, pattern="RLWM",ignore.case = TRUE))){
    EEGdata_preg[i,7]<-"Yes"
  } else {
    EEGdata_preg[i,7]<-"X"
  }
}
for (i in seq_along(eegdirs_preg)) 
{
  if (0!=length(list.files(eegdirs_preg[i], all.files=TRUE, pattern="MMN",ignore.case = TRUE))){
    EEGdata_preg[i,8]<-"Yes"
  } else {
    EEGdata_preg[i,8]<-"X"
  }
}
for (i in seq_along(eegdirs_preg)) 
{
  if (0!=length(list.files(eegdirs_preg[i], all.files=TRUE, pattern="SST",ignore.case = TRUE))){
    EEGdata_preg[i,9]<-"Yes"
  } else {
    EEGdata_preg[i,9]<-"X"
  }
}
for (i in seq_along(eegdirs_preg)) 
{
  if (0!=length(list.files(eegdirs_preg[i], all.files=TRUE, pattern="Emopics",ignore.case = TRUE))){
    EEGdata_preg[i,10]<-"Yes"
  } else {
    EEGdata_preg[i,10]<-"X"
  }
}

# See if there is a file in files_preg_3D matching the participant id, insert Yes if TRUE or X if false
EEGdata_preg <- EEGdata_preg %>% 
  rownames_to_column("id") %>% 
  rowwise() %>% 
  mutate(`3D_Scan` = ifelse(any(str_detect(files_preg_3D, id)), "Yes", "X"))

### NOW WE WANT TO REPEAT A SIMILAR APPROACH BUT FOR BEHAVIORAL DATA FILES, ADDING THEM TO THE DATAFRAME ## 
behaviorfolder_preg <- "/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_behavior/eeg_behavior_pregnancy"
behaviordir_preg <- list.dirs(path = behaviorfolder_preg, full.names = TRUE, recursive = FALSE) #list of all directories. 1 per participant
behaviordir_preg <- behaviordir_preg[ grepl("68*", behaviordir_preg)]
behaviorfiles_preg <- list.files(behaviordir_preg, all.files=TRUE, pattern="*.csv|*.mat") #list all files 

behaviorRLWM_preg <- list.files(behaviordir_preg, all.files=TRUE, pattern="rlwmpst", ignore.case = TRUE)
behaviorSST_preg <- list.files(behaviordir_preg, all.files=TRUE, pattern="SST", ignore.case = TRUE)
behaviorEmopics_preg <- list.files(behaviordir_preg, all.files=TRUE, pattern="Emopics", ignore.case = TRUE)

behaviorRLWM_preg<-gsub("ID","", behaviorRLWM_preg)
behaviorRLWM_preg_df <- behaviorRLWM_preg %>% 
  as_tibble() %>% 
  separate(value, into = c(NA,"stripped_ID",NA,"timepoint",NA), sep = "_") %>% 
  mutate(timepoint = parse_number(timepoint),
         id_tp = str_glue("{stripped_ID}_{timepoint}")) #

behaviorSST_preg_df <- behaviorSST_preg %>% 
  as_tibble() %>% 
  separate(value, into = c(NA,NA,"stripped_ID","timepoint",NA), sep = "_") %>% 
  mutate(timepoint = parse_number(timepoint),
         id_tp = str_glue("{stripped_ID}_{timepoint}")) #
behaviorEmopics_preg <-gsub("sub","", behaviorEmopics_preg)
behaviorEmopics_preg_df <- behaviorEmopics_preg %>% 
  as_tibble() %>% 
  separate(value, into = c(NA,NA,"stripped_ID",NA,"timepoint", NA), sep = "_") %>% 
  mutate(timepoint = parse_number(timepoint),
         id_tp = str_glue("{stripped_ID}_{timepoint}")) #


# See if there is a file of RLWM, SST, and Emopics  matching the participant id and time point, insert Y if TRUE or X if false

EEGdata_preg <- EEGdata_preg %>% 
  mutate(behaviorRLWM = ifelse(ID_Timepoint %in% behaviorRLWM_preg_df$id_tp, "Yes", "X"))

EEGdata_preg <- EEGdata_preg %>% 
  mutate(behaviorSST = ifelse(ID_Timepoint %in% behaviorSST_preg_df$id_tp, "Yes", "X"))
#         id_tp = str_glue("{stripped_ID}_{Timepoint_integer}"),
#         id_tp = str_glue("{stripped_ID}_{Timepoint_integer}")
#Timepoint_integer = as.numeric(Timepoint_integer)
EEGdata_preg <- EEGdata_preg %>% 
  mutate(behaviorEmopics = ifelse(ID_Timepoint %in% behaviorEmopics_preg_df$id_tp, "Yes", "X"))

EEGdata_preg <- merge(pregnancy_logg,EEGdata_preg, by.x = "logg_ID", by.y = "id", all.x = TRUE)
EEGdata_preg <- EEGdata_preg %>% arrange(desc(logg_date))
#Write the table as a .xlsx file with current date added
currentDate <- Sys.Date()
preg_filename <- paste("/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/EEG_DataLists/Pregnancy_EEGdatalist",currentDate,".xlsx",sep="")
write.xlsx(EEGdata_preg,file=preg_filename, rowNames=TRUE)

