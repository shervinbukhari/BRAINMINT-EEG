% Script for transforming BRAINMINT EEG files to BIDS format
% Shervin Bukhari 04.10.23
function convert_to_bids_eeg(dir_in)
 dir_out = [dir_in 'bids'];
    if ~exist(dir_out)
        mkdir(dir_out)
    end
% set EEGlab
EEGLABfolder = '/cluster/projects/p33/groups/imaging/BRAINMINT/eeg/analysis/toolboxes/eeglab2022.0/';

addpath(EEGLABfolder);
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

pop_editoptions('option_single', false, 'option_savetwofiles', false);

% Turn of parallel pools
ps = parallel.Settings;
ps.Pool.AutoCreate = false;

% List all bdf files within the subject directories that start with 6*
files = [];
files = dir([dir_in '6*/*.bdf']);
if isempty(files)
    error('Files is empty')
end
%% BIDS Conversion
% Loop through all files
% extract the ID and paradigm, open the bdf file and extract date

strange_files = double.empty(0,2); % safety to identify those that don't follow the naming coventions
for i = 1 :size(files,1)
    fprintf('\n Running rename_bdf pipeline on %s \n', files(i).name)
    try % safety net in case something is wrong with the file, e.g. empty file
        split_name = split(files(i).name, '_');
        
        % find ID
        files(i).ID = char(split_name(1));
        
        % find the paradigm
        if contains(files(i).name,'Emopics','IgnoreCase',1)
            files(i).paradigm = 'Emopics';
        elseif contains(files(i).name,'SST','IgnoreCase',1)
            files(i).paradigm = 'SST';
        elseif contains(files(i).name,'RLWM','IgnoreCase',1)
            files(i).paradigm = 'RLWM';
        elseif contains(files(i).name,'MMN','IgnoreCase',1)
            files(i).paradigm = 'MMN';
        elseif contains(files(i).name,'REST_closed','IgnoreCase',1)
            files(i).paradigm = 'REST_closed';
        elseif contains(files(i).name,'REST_open','IgnoreCase',1)
            files(i).paradigm = 'REST_open';
        else
	    fprintf(['\n Something went wrong here: ' files(i).name '\n']);
            strange_files = [strange_files;[i 1]];
            continue
        end
        
        % get time point
        if (str2num(files(i).ID) >= 68000)
            if contains(files(i).name,'_01_') || contains(files(i).name,'_1_') || contains(files(i).name,'_pre_')
                files(i).timepoint = 'pre';
            elseif contains(files(i).name,'_02_') || contains(files(i).name,'_2_') || contains(files(i).name,'_mid2_')
                files(i).timepoint = 'mid2';
            elseif contains(files(i).name,'_03_') || contains(files(i).name,'_3_') || contains(files(i).name,'_post1_')
                files(i).timepoint = 'post1';
            elseif contains(files(i).name,'_04_') || contains(files(i).name,'_4_') || contains(files(i).name,'_post2_')
                files(i).timepoint = 'post2';
            elseif contains(files(i).name,'control', 'IgnoreCase', true)
                files(i).timepoint = 'control18m';
            else
		fprintf(['\n Something went wrong here: ' files(i).name '\n']);
                strange_files = [strange_files;[i 2]];
            end
        elseif (str2num(files(i).ID) >= 65000 && str2num(files(i).ID) < 68000)
            % if adolescents (CHANGE WHEN MORE TIME POINTS ARE COLLECTED)
            files(i).timepoint = '01';
        else
	    fprintf(['\n Something went wrong here: ' files(i).name '\n']);
            strange_files = [strange_files; [i 3]];
            continue
        end
        
        % load file and extract date
        rawfile = [ files(i).folder filesep files(i).name ];
        EEG = pop_biosig(rawfile, 'ref',38);
        
        year = num2str(EEG.etc.T0(1));
        
        month = '';
        if length(num2str(EEG.etc.T0(2)))== 1
            month = ['0' num2str(EEG.etc.T0(2))];
        else
            month = num2str(EEG.etc.T0(2));
        end
        
        day = '';
        if length(num2str(EEG.etc.T0(3)))== 1
            day = ['0' num2str(EEG.etc.T0(3))];
        else
            day = num2str(EEG.etc.T0(3));
        end
        
        files(i).recdate = [year month day];
        clear('EEG')
    catch
        fprintf(['\n Something went wrong here: ' files(i).name '\n']);
        strange_files = [strange_files;[i 4]];
    end
    
end

% save files names of the strange files
fileID = fopen([dir_in 'recent_import/strange_files' filesep 'Strange_filesBIDS' date '.txt'], 'w');
fprintf(fileID, 'strange_files\t error\n');  
for i = 1:size(strange_files,1)
    j = strange_files(i,1);
    fprintf(fileID, '%s\t %d\n', files(j).name, strange_files(i,2)); 
end
fclose(fileID);

% move strange files into separate folder
for i = 1:size(strange_files,1)
    j = strange_files(i,1);
    fprintf('\n Moving %s to strange files folder \n', files(j).name)
    movefile([files(j).folder filesep files(j).name], [dir_in 'recent_import/strange_files']);
end

% remove strange files
files(strange_files(:,1)) = [];
%% Move files
% movefile source destination

fileoverview = struct('old_names','old', 'new_names', 'new');
for i = 1:size(files,1)
    orig_file = [files(i).folder filesep files(i).name];
    new_subject_directory = [dir_out filesep 'sub-' files(i).ID];
    new_session_directory = [dir_out filesep 'sub-' files(i).ID filesep 'ses-' files(i).timepoint];
    new_eeg_directory =  [dir_out filesep 'sub-' files(i).ID filesep 'ses-' files(i).timepoint filesep 'eeg'];
    new_filename = ['sub-' files(i).ID '_' files(i).timepoint '_' files(i).paradigm '_' files(i).recdate '_eeg.bdf'];
     
    if ~exist(new_subject_directory)
        mkdir(new_subject_directory)
    end
    
    if ~exist(new_session_directory)
        mkdir(new_session_directory)
    end
    
    if ~exist(new_eeg_directory)
    mkdir(new_eeg_directory)
    end
    
    % check if that file already exists, can handle up to two duplications
    if isfile([new_eeg_directory filesep new_filename])
        new_filename = ['sub-' files(i).ID '_' files(i).timepoint '_' files(i).paradigm '_' files(i).recdate '_dupl_eeg.bdf'];
        if isfile([new_directory filesep new_filename])
            new_filename = ['sub-' files(i).ID '_' files(i).timepoint '_' files(i).paradigm '_' files(i).recdate '_dupl2_eeg.bdf'];
        end
    end
    
    %copyfile(orig_file, [new_directory filesep new_filename])
    movefile(orig_file, [new_eeg_directory filesep new_filename]);
    
    % create overview of old and new names
    fileoverview(i).old_names = files(i).name;
    fileoverview(i).new_names = new_filename;
end

% save fileoverview
fileID = fopen([dir_in filesep 'recent_import' filesep 'BIDS_conversion_EEG_files_overview_' date '.txt'], 'w');
fprintf(fileID, 'old_names\t new_names\n');  
for i = 1:size(fileoverview,2) 
    fprintf(fileID, '%s\t %s\n', fileoverview(i).old_names, fileoverview(i).new_names); 
end
fclose(fileID);
% Generate auxiliary bids files 
channelloc_to_tsv(EEG)
eeg_writechanfile(EEG, [new_eeg_directory filesep 'sub-' files(i).ID '_' files(i).timepoint '_' files(i).paradigm '_' files(i).recdate])
eeg_writeeventsfiles(EEG, [new_eeg_directory filesep 'sub-' files(i).ID '_' files(i).timepoint '_' files(i).paradigm '_' files(i).recdate]) 
electrodes_to_tsv(EEG)

end