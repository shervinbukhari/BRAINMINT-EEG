% Standard template for script to be used to preprocessing of BRAINMINT
% data. The script has many parameters that can be tweaked and adjusted
% depending on the specific paradigm or anaylsis. This provides a basic
% structure.
% By Shervin H. Bukhari, November 2023
% Please Note:
% Create your own copy of this function before you start editing or changing paths.
%
% Code for including and dealing with EOGs in ICA decomposition is unfinished
% ICA is performed on continuous data.
% Summary files are created. One file contains names of interpolated
% electrodes
%
function Preprocess(subject_session_ID)
%% Define paths, paradigm and load EEGlab
addpath('/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/.../scripts/...') % Must be edited
rmpath(genpath('/gpfs/colossus01/software/matlab/R2017a/toolbox/eeglab13_4_4b/'))
EEGLABfolder = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/toolboxes/eeglab2022.0/';
addpath(EEGLABfolder);

rawfile = char(subject_session_ID);
parts = strsplit(subject_session_ID, '_');
subject = parts{1};
TP = parts{2};
paradigm = 'InsertName'; %REST_closed, REST_open, RLWM, MMN, SST or Emopics

% check if adolescent or pregnancy
if strcmp(subject(1:2),'65')
    % adolescent
    group = 'adolescence';
    rawdir = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_adolescents';
%Define your output folder for adolescence
    outdir = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/.../adolescence';
    subject_session = [subject, '_', TP];
elseif strcmp(subject(1:2),'66')
    % adolescent
    group = 'adolescence';
    rawdir = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_adolescents';
%Define your output folder for adolescence
    outdir = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/.../adolescence';
    subject_session = [subject, '_', TP];
elseif strcmp(subject(1:2),'68')
    % pregnancy
    group = 'pregnancy';
    rawdir = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_pregnancy';
    
%Define your output folder for pregnancy
    outdir = '/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/.../pregnancy';
    subject_session = [subject, '_', TP];
end

% Get a list of files that match the pattern
fileList = ls([rawdir, '/', subject_session, '/', subject_session, '_' paradigm '*.bdf']);%insert paradigm name before asterix

% Convert the file list to a cell array of strings
fileList = strsplit(fileList);

% Find the index of the filename containing the string 'merged.set'
mergedIdx = find(contains(fileList, 'merged.set'));

% If there is a file containing 'merged', select it
if ~isempty(mergedIdx)
    rawfile = fileList{mergedIdx};
else
% Otherwise, select the first file in the list
    rawfile = fileList{1};
end
% Trim any whitespace from the selected file name
rawfile = strtrim(rawfile);

[ALLEEG , ~, CURRENTSET ALLCOM] = eeglab;
pop_editoptions('option_single', false, 'option_savetwofiles', false);
%% PARAMETER SETTINGS%%
%Please edit these parameters instead of the main body of the code for
%readability. Do not add semi-colon, so that they are printed in the
%slurm-log
fprintf('\n\n\n**** Selected Preprocessing parameter settings for %s %s %s ****\n\n\n', subject, TP, paradigm);
TriggerLatency = 0 %.02 % seconds latency
TriggerstoMove = {'trigger1','trigger2'} % Insert the EEG.event.type code of the triggers you want to displace by TriggerLatency.
EEGchannel_indices = 1:64 % EEG channel indices, do not change
allScalpchannel_indices = 1:66 % Indices of all channels including the external channels around eyes
allelectrophyschannel_indices = 1:72
HighPassFrequency = 0.1 % Pass frequencies above
LowPassFrequency = 100 % Pass frequencies below
EMG_highpass = 28 % Highpass threshold EMG
EMG_lowpass = 500 % Lowpass threshold EMG
ECG_highpass = 0.1 % Highpass threshold ECG
Demean = 'Yes' % Demean the continuous signal or not
DownSamplingRate = 512 %Which frequency to downsample data to.
LineNoiseFrequency = [50 100 150 200 250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000] % HOW PARANOID ARE YOU ABOUT LINE NOISE HARMONICS?
IncludeEOGs ='No' %'Yes' % Whether or not to include EOGs in the ICA decomposition, No also excludes all other external sensors
% FOR CLEAN_RAWDATA()
fprintf('Parameters for clean_raw()')
FlatLineCriterion = 'off' % Max. tolerated flatline duration (seconds)
ChannelCriterion = 0.8 % Minimum channel correlation. If a channel is correlated at less than this value to a reconstruction of it based on other channels, it is considered abnormal in the given time window
LineNoiseCriterion = 4 % If a channel has more line noise relative to its signal than this value, in standard deviations based on the total channel population, it is considered abnormal.
BurstCriterion = 'off'% 20 %Standard deviation cutoff for removal of bursts (via ASR). Data portions whose variance is larger than this threshold relative to the calibration data are considered missing data and will be removed.
WindowCriterion = 'off'% 0.25 %Criterion for removing time windows that were not repaired completely. This mayhappen if the artifact in a window was composed of too many simultaneous uncorrelated sources (for example, extreme movements such as jumps). This is the maximum fraction of contaminated channels that are tolerated in the final output data for each considered window
WindowCriterionTolerances = 'off' %[-Inf 7] %These are the power tolerances outside of which a channel in the finaloutput data is considered "bad", in standard deviations relative to a robust EEG power distribution (lower and upper bound). Any time window in the final (repaired) output which has more than the tolerated fraction (set by the WindowCriterion parameter) of channel with a power outside of this range will be considered incompletely repaired and will be removed from the output.
HighPass = 'off' %Transition band for the initial high-pass filter in Hz. Only use if data is not already HP filtered.
Distance = 'Euclidian' %Unable to find documentation of this parameter, but space is Euclidian
BurstRejection = 'off' %If 'on' reject portions of data containing burst instead of correcting them using ASR.
ICflagsettings = [NaN NaN;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN;NaN NaN] 
% Rejects components marked as 90% probability of eye, muscle 
% Format is [Brain; Muscle; Eye; Heart; Line Noise; Channel Noise, Other]
% First digit indicates the probability threshold and 1 = flags if probability is atleast P| 0 = flags probability is less than P. P=probability threshold             
% [NaN NaN;0.9 1;0.9 1;0.9 1;NaN NaN;NaN NaN;NaN NaN]
trimamplitudeThreshold = 500 % µV
trimpointSpreadWidth = 500 % ms
fprintf('Epoching and Epoch Rejection settings')
window_in_seconds       = [0 0] % [pre post] stimulus in seconds
events_of_interest = (".." |"..") ; %insert trigger codes
event_list              = cell(unique({EEG.event.type}));
indx                    = contains(event_list,events_of_interest);
epoch_names             = event_list(indx);
rejection_window = [0 0] % start and end in seconds, specifies time window to mark for rejection with pop_eegthresh
rejection_deviation = [-0 0] % negative and positive potential deviation in µV to mark for rejection with pop_eegthresh
fprintf('-----------------------------------------------------------------------')
%% LOAD EEG FILE
fprintf('\n\n\n**** Loading %s dataset for participant %s %s %s ****\n\n\n', subject,TP, paradigm);
if contains([rawfile],'merged.set')
    EEG=pop_loadset([rawfile], 'ref', 38);
else
    EEG = pop_biosig([rawfile],'ref', 38);
        if length(EEG.chanlocs) == 75
        EEG = pop_chanedit(EEG, 'lookup',[EEGLABfolder 'plugins/dipfit/standard_BESA/standard-10-5-cap385.elp'],'load',{'/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/BRAINMINT.ced','filetype' 'autodetect'});
        else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG = pop_chanedit(EEG, 'lookup',[EEGLABfolder 'plugins/dipfit/standard_BESA/standard-10-5-cap385.elp'],'load',{'/ess/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/BRAINMINT_withoutGSR.ced','filetype' 'autodetect'});
        end
end

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

if strcmp(IncludeEOGs, 'Yes')
% Do nothing
else %Remove them from the channel list
%EEG = pop_select( EEG, 'channel',{'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'});    
end

EEG.group = group;
EEG.subject = subject; 
%% Create bipolar HEOG channel and run initial EOG preprocessing
if EEG.nbchan == 75
        EEG.data(76,:)=EEG.data(66,:)-EEG.data(65,:);
        EEG.chanlocs(76).labels = 'HEOG';
else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG.data(74,:)=EEG.data(66,:)-EEG.data(65,:);
        EEG.chanlocs(74).labels = 'HEOG';
end%% Create bipolar EMG channels and run initial EMG preprocessing
if EEG.nbchan == 75
        EEG.data(77,:)=EEG.data(68,:)-EEG.data(67,:);
        EEG.chanlocs(77).labels = 'Left_EMG';
        EEG.data(78,:)=EEG.data(70,:)-EEG.data(69,:);
        EEG.chanlocs(78).labels = 'Right_EMG';
else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG.data(75,:)=EEG.data(68,:)-EEG.data(67,:);
        EEG.chanlocs(75).labels = 'Left_EMG';
        EEG.data(76,:)=EEG.data(70,:)-EEG.data(69,:);
        EEG.chanlocs(76).labels = 'Right_EMG';
end   

%% Band pass filter the data between 28 and 500 Hz
EEG = pop_eegfiltnew(EEG, 'locutoff',EMG_highpass,'plotfreqz',0,'channels',{'Left_EMG','Right_EMG'});
EEG = pop_eegfiltnew(EEG, 'hicutoff',EMG_lowpass,'plotfreqz',0,'channels',{'Left_EMG','Right_EMG'});

%% Calculate the root mean square of the signal using moving 40-point (20 ms) windows (39-point overlap and zero padding) 
if EEG.nbchan == 75
        EEG.data(77,:) = rms_PPI(EEG.data(77,:), 40,39,1);
        EEG.data(78,:) = rms_PPI(EEG.data(78,:), 40,39,1);
else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG.data(75,:) = rms_PPI(EEG.data(75,:), 40,39,1);
        EEG.data(76,:) = rms_PPI(EEG.data(76,:), 40,39,1);
end   

%% Create bipolar ECG channel and run initial ECG preprocessing
if EEG.nbchan == 75
        EEG.data(79,:)=EEG.data(72,:)-EEG.data(71,:);
        EEG.chanlocs(79).labels = 'ECG';    
        EEG = pop_eegfiltnew(EEG, 'locutoff',ECG_highpass,'plotfreqz',0,'channels',{'ECG'});

else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG.data(77,:)=EEG.data(72,:)-EEG.data(71,:);
        EEG.chanlocs(77).labels = 'ECG';
        EEG = pop_eegfiltnew(EEG, 'locutoff',ECG_highpass,'plotfreqz',0,'channels',{'ECG'});
end   

%% Create bipolar GSR channel and run initial GSR preprocessing
if EEG.nbchan == 75
        EEG.data(80,:)=EEG.data(74,:)-EEG.data(73,:);
        EEG.chanlocs(80).labels = 'GSR';     
% else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
%         EEG.data(78,:)=EEG.data(74,:)-EEG.data(73,:);
%         EEG.chanlocs(78).labels = 'GSR';
end   

%% Create scaled respiration channel
if EEG.nbchan == 75
        EEG.data(81,:)=EEG.data(75,:)/100;
        EEG.chanlocs(81).labels = 'RESP_scaled';     
else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG.data(78,:)=EEG.data(73,:)/100;
        EEG.chanlocs(78).labels = 'RESP_scaled';
end   

%% Update number of channels 
if EEG.nbchan == 75
        EEG.nbchan = 81;  
else %If the recording has fewer than 75 electrodes, assume GSRs are missing and use this coordinate template instead
        EEG.nbchan = 78;  
end   
%% APPLY FILTERING - HIGH PASS, LOWPASS AND LINE NOISE REMOVAL
fprintf('\n\n\n**** Applying line noise and high pass filtering to participant %s %s %s ****\n\n\n', subject,TP, paradigm);
if strcmp(Demean, 'Yes')
EEG = pop_rmbase( EEG, [], [], [1:72]);
end
%Clean line noise
EEG = pop_cleanline(EEG,'linefreqs',LineNoiseFrequency,'ChanCompIndices', EEGchannel_indices);

% High Pass filtering
EEG = pop_eegfiltnew(EEG, 'locutoff',HighPassFrequency,'channels',{EEG.chanlocs(allelectrophyschannel_indices).labels});

% Downsample EEG files to 512 Hz
fprintf('\n\n\n**** Downsampling for participant %s %s %s ****\n\n\n', subject,TP, paradigm);
EEG = pop_resample( EEG, DownSamplingRate);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

% Select scalp electrodes
all_EEG = EEG; 
EEG = pop_select( EEG, 'channel',{'Fp1','AF7','AF3','F1','F3','F5','F7','FT7','FC5','FC3','FC1','C1','C3','C5','T7','TP7','CP5','CP3','CP1','P1','P3','P5','P7','P9','PO7','PO3','O1','Iz','Oz','POz','Pz','CPz','Fpz','Fp2','AF8','AF4','AFz','Fz','F2','F4','F6','F8','FT8','FC6','FC4','FC2','FCz','Cz','C2','C4','C6','T8','TP8','CP6','CP4','CP2','P2','P4','P6','P8','P10','PO8','PO4','O2'});


%% Remove bad channels
fprintf('\n\n\n**** Bad channel detection for participant %s %s %s ****\n\n\n', subject, TP, paradigm);
[EEG.urchanlocs] = deal(EEG.chanlocs); % Keep original channels
EEG = pop_reref(EEG,[]);
%EEG = pop_reref(EEG, {EEG.chanlocs(EEGchannel_indices).labels}, 'keepref', 'on'); % average referencing (excluding HEOGs)
orignal_timelength=EEG.xmax;
EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',FlatLineCriterion,'ChannelCriterion',ChannelCriterion,'LineNoiseCriterion',LineNoiseCriterion,'Highpass',HighPass,'BurstCriterion',BurstCriterion,...
    'WindowCriterion',WindowCriterion,'BurstRejection',BurstRejection,'Distance',Distance,'WindowCriterionTolerances',WindowCriterionTolerances);

    if(isfield(EEG.etc,'clean_channel_mask'))
        clean_channel_mask = EEG.etc.clean_channel_mask;
    else
        clean_channel_mask = ones(64,1);
    end
EEG.etc.badchannels= EEG.urchanlocs(clean_channel_mask == 0);
data_rank=length(EEG.chanlocs);

%% INTERPOLATION and REREFERENCE
% Interpolate removed channels 
EEG = pop_interp(EEG, EEG.urchanlocs);    
EEG = pop_reref(EEG, []); 

% Reference EOGs to average reference
if strcmp(IncludeEOGs, 'Yes') % Not working
%pop_select(EEG, 'channel', {'LO1', 'LO2'});
else %Do nothing
end
%% Reconstruct continuous data with selected external channels 
all_EEG.data(1:64,:)=EEG.data(1:64,:);

all_EEG.urchanlocs = EEG.urchanlocs;

if(isfield(EEG.etc,'clean_channel_mask'))
    all_EEG.etc.badchannels= EEG.etc.badchannels;
end

EEG = all_EEG;

%% Select clean continuous data to run ICA 
%% Run trimOutlier() Rejects +/- 500 ms around any datapoint exceeding 500 mV across the 64 scalp channels
fprintf('\n\n\n**** Applying trimOutlier to data of participant %s %s %s ****\n\n\n', subject, TP, paradigm);
EEG = trimOutlier_timepoints(EEG, EEGchannel_indices, trimamplitudeThreshold, trimpointSpreadWidth);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
percentage_timeretained= EEG.xmax/orignal_timelength;

%% ICA
% Run ICA and IC Label
fprintf('\n\n\n**** Running ICA on participant %s %s %s ****\n\n\n', subject,TP, paradigm);
EEG = pop_runica(EEG, 'icatype', 'picard', 'maxiter',500,'mode','standard','pca',data_rank, 'chanind', EEGchannel_indices);
%%
%% Reconstruct continuous data with ICA decomposition based on clean data
all_EEG.icaweights=EEG.icaweights;
all_EEG.icachansind=EEG.icachansind;
all_EEG.icaact=EEG.icaact;
all_EEG.icasphere=EEG.icasphere;
all_EEG.icawinv=EEG.icawinv;
all_EEG.icasplinefile=EEG.icasplinefile;

EEG = all_EEG;
EEG = pop_iclabel(EEG, 'default');
EEG = pop_icflag(EEG, ICflagsettings); 
[pr, ind] = max(EEG.etc.ic_classification.ICLabel.classifications, [], 2); %for summary file
number_ICs = length(EEG.etc.ic_classification.ICLabel.classifications); %for summary file
 IC_reject = find(EEG.etc.ic_classification.ICLabel.classifications(:,2) >= 0.90 | ...
    EEG.etc.ic_classification.ICLabel.classifications(:,3) >= 0.90); %If
%    you want more complicated criteria, such as reject everything that is
%    not x% brain y% other, you can use adjust IC_reject and replace [] with IC_reject
%    in pop_subcomp
fprintf('\n\n\n**** Rejecting the following ICs for participant %s %s %s ****\n\n\n', subject, TP, paradigm);
fprintf('%5d', IC_reject);
fprintf('\n\n\n');
EEG = pop_subcomp( EEG, [], 0);

%Low Pass filtering
EEG = pop_eegfiltnew(EEG, 'hicutoff',LowPassFrequency,'plotfreqz',0,'channels',{EEG.chanlocs(allelectrophyschannel_indices).labels});
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );

%% Epoching
fprintf('\n\n\n**** Epoching for participant %s %s %s ****\n\n\n', subject,TP, paradigm);

%move event markers 20 ms forward (to control for hardware latency)
for index = 1:length(EEG.event)
    if ismember(EEG.event(index).type, TriggerstoMove)
    EEG.event(index).latency = EEG.event(index).latency+TriggerLatency*EEG.srate; %moves them latency*EEG.srate samples forward
    end
end

EEG = eeg_checkset(EEG, 'eventconsistency');
EEG = pop_epoch( EEG, epoch_names, window_in_seconds, 'epochinfo', 'yes', 'newname', [ subject '_' TP '_' paradigm '_BRAINMINT_epoched' ]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
number_epochs = EEG.trials;



%% EPOCH REJECTION
fprintf('\n\n\n**** Running epoch rejection for participant %s %s %s ****\n\n\n', subject, TP, paradigm);

EEG = eeg_checkset( EEG );
EEG = pop_eegthresh(EEG,1,EEGchannel_indices, rejection_deviation(1), rejection_deviation(2), rejection_window(1), rejection_window(2) , 2, 0);
EEG = eeg_rejsuperpose( EEG, 1, 1, 1, 1, 1, 1, 1, 1);
rejected_epochs = sum(EEG.reject.rejglobal);
EEG = pop_rejepoch( EEG, [find(EEG.reject.rejglobal)] ,0);

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

EEG = pop_saveset( EEG, 'filename',[subject '_' TP '_' paradigm 'BRAINMINT_Preprocessed.set'],'filepath',outdir); %change name to your preference

[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
%% %% Create summary file
percentage_trialsrejected = rejected_epochs/number_epochs;
fprintf('\n\n\n**** Creating summary file for participant %s %s %s ****\n\n\n', subject, TP, paradigm);
summary = [str2num(subject) , length(EEG.etc.badchannels),...
     number_ICs, length(IC_reject), numel(find(ind(IC_reject) == 2)), numel(find(ind(IC_reject) == 3)), numel(find(ind(IC_reject) == 4))...
    , percentage_timeretained,percentage_trialsrejected, rejected epochs];
Interpolated_Channels = strjoin({EEG.etc.badchannels.labels},',');
fileID = fopen([outdir filesep subject '_' TP '_' paradigm '_BRAINMINT_Preprocessing_summary.txt'], 'w');
fprintf(fileID, 'ID\t Num_InterpolatedChan\t NumComp\t RejComp\t Muscle\t Eye\t Heart\t PctTimeRetained(trim)\t PctRejectedEpochs\t Num_RejectedEpochs\n');  
fprintf(fileID, '%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\t', summary);
fclose(fileID);
fileID = fopen([outdir filesep subject '_' TP '_' paradigm '_BRAINMINT_Preprocessing_InterpolatedChannels.txt'], 'w')
fprintf(fileID, 'Num_InterpolatedChan\n');  
fprintf(fileID, '%s\t ', Interpolated_Channels);
fclose(fileID);

fprintf('\n\n\n**** PREPROCESSING FINISHED FOR PARTICIPANT %s %s %s ****\n\n\n', subject,TP,paradigm);
clear all;
restoredefaultpath
rehash toolboxcache
end