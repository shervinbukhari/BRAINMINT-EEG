%% Script for merging data files from the same participant and session if recordings has been interrupted midway 
% 
% NB! Merged file is stored as a .set-file. Storing as eeg not possible due
% to bug in pop_writeeeg causing the eventstructure to be messed up.
% the EEG recording was accidentally stopped and resumed during a task,
% resulting in two files.
% This script will help you identify cases where there are two
% matches for the same file name based on task, paradigm and timepoint
% but requires manual steps at the end to actually merge the files.
% NB!: Files cannot be merged haphazardly. Sometimes a task has been
% restarted completely.
%
%

% Shervin Bukhari 26.09.2023
%%
restoredefaultpath;
EEGLABfolder = '/cluster/p/p33/cluster/groups/imaging/BRAINMINT/eeg/analysis/toolboxes/eeglab2022.0';
addpath(genpath(EEGLABfolder));
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
pop_editoptions('option_single', false, 'option_savetwofiles', false);

dir_eeg_pregnancy= '/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_pregnancy/';
dir_eeg_adolescence ='/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/eeg_raw_data/eeg_raw_adolescents/';
% List all eeg files within the subject directories
eeg_files_pregnancy = [dir([dir_eeg_pregnancy '6*/*.bdf']); dir([dir_eeg_pregnancy '6*/*.set'])] ;
eeg_files_adolescence = [dir([dir_eeg_adolescence '6*/*.bdf']); dir([dir_eeg_pregnancy '6*/*.set'])];

% concatenate the lists
all_eeg_files = [eeg_files_pregnancy;eeg_files_adolescence];

% Remove rest eyes/closed, since next section cannot distinguish them
eeg_files = all_eeg_files;
for i=numel(all_eeg_files):-1:1
    if contains(all_eeg_files(i).name,'REST') 
    eeg_files(i) = [];
    end
end
%Remove participants/sessions where there already is a merged file from the list
%Search for files containing 'merged'
merged_files = struct('name', {}, 'folder', {});
for i = 1:length(eeg_files)
    if contains(eeg_files(i).name, 'merged')
        merged_files(end+1).name = eeg_files(i).name;
        merged_files(end).folder = eeg_files(i).folder;
    end
end 
% Remove any file that has the same path as a previously merged from list
for i = length(eeg_files):-1:1
    if any(strcmp({merged_files.folder},eeg_files(i).folder))
    eeg_files(i) = [];
    end
end
%% Find cases where there are two files for a paradigm
matching_indices = cell(size(eeg_files));
% Look for cases in which ID, time point and pardigm occur twice
for i=1:numel(eeg_files)
    split_name = split(eeg_files(i).name, '_');
    eeg_files(i).id = char(split_name(1));
    eeg_files(i).timepoint = char(split_name(2));
    eeg_files(i).paradigm = char(split_name(3));
end
for i=1:numel(eeg_files)
    curr_file = eeg_files(i).name;
    % Loop over each file again to find matching file names
    for j = 1:numel(eeg_files)
        % Check if file names match and indices are different
        if contains(curr_file, eeg_files(j).id) && contains(curr_file, eeg_files(j).timepoint) && contains(curr_file, eeg_files(j).paradigm) && i ~= j
            % Add the index to the matching indices cell array
            matching_indices{i} = [matching_indices{i}, j];
        end
    end
end
% remove empty fields
matching_indices = matching_indices(~cellfun('isempty', matching_indices));

%% Merge the files and create README
%This part requires manual work.
% The matching_indices variable lists each file index in eeg_files which
% has a potential duplicate file that needs it needs to be merged with.
% However not all files with duplicates should be merged, sometimes the
% paradigm has been restarted or reran. Therefore, one has to manually
% check the data collection logg for each participant to see if there are any comments from
% data collection that gives additional information.
% List of files to excluded from merging: 
% '68082_mid2_SST_20210204.bdf' & '68082_mid2_SST_20210204_run2.bdf'
% '68103_pre_RLWM_20210329.bdf' & '68103_pre_RLWM_20210224.bdf'
% '68196_pre_SST_20210816.bdf' & '68196_pre_SST_20210816_run2.bdf'
% '68289_pre_RLWM_20220111.bdf' & '68289_pre_RLWM_2_20220111.bdf'
% '68444_pre_SST_20220510.bdf' & 68444_pre_SST_20220510_part2.bdf'
% '65062_01_RLWM_20211024.bdf' & '65062_01_RLWM_20211024_run2.bdf'
% '65104_01_MMN_20210811.bdf' & '65104_01_MMN_20210811_triggers_fixed.bdf'
% '65160_01_RLWM_20210908.bdf' & '65160_01_RLWM_20210908_run2.bdf'
% To merge two EEG files input the eeg_files index of the two files in x
% and y below. The code automatically organizes the data chronologically
% when merging and gives the resulting file the name of the first file with
% "_merged" added to the end.

clear EEG, clear INEEG1 INEEG2
x=3600;%index original file 
y=3599;%index second file
filestomerge{1}=[eeg_files(x).folder, '/' eeg_files(x).name]; %original file
filestomerge{2}=[eeg_files(y).folder, '/' eeg_files(y).name]; %second file

INEEG1 = pop_biosig(filestomerge{1});
INEEG2 = pop_biosig(filestomerge{2}); 
%Makes sure the files are merged in order of creation so that the data is
%chronologically ordered.
split_name = split(filestomerge{1}, '.');
if datenum(INEEG1.etc.T0) < datenum(INEEG2.etc.T0)
EEG=pop_mergeset(INEEG1,INEEG2);
elseif datenum(INEEG2.etc.T0) < datenum(INEEG1.etc.T0)
EEG=pop_mergeset(INEEG2,INEEG1);
end
merged_filename = [split_name{1},'_merged.set'] ;
%Save the merged file
pop_saveset(EEG, 'filename' ,merged_filename, 'version', '7.3');
%create README file 
%Specify the file path and name
filepath = eeg_files(x).folder; 
dateStr =datestr(date, 'yyyy-mm-dd');

%Open the file for writing
fid = fopen([filepath filesep 'README'], 'w');
%Check if the file was successfully opened
if fid == -1
    error('Unable to create file');
end

%Write the content to the file
fprintf(fid, 'The files %s and %s have been merged into %s using EEGLAB on %s. \n', ...
eeg_files(x).name, eeg_files(x).name, merged_filename, dateStr); %define the names
%Close the file
fclose(fid);
