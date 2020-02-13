






%%

clear all

system = 1 %0 (1 = MAC, 0 = Windows)
use_costrap = 0 %1 = do
use_unfold = 0 % 1 = do unfold 0 = don't

if system == 1
 
cd('/Volumes/methlab/Neurometric/Antisaccades/ALL/eeglab14_1_2b/')
addpath('/Volumes/methlab/Neurometric/Antisaccades/ALL/')
addpath('/Volumes/methlab/Neurometric/Antisaccades/ALL/FastICA_25')
% addpath('tools');
eeglab;
close
%addpath('/Users/mplome/Desktop/Tools')
table_old =readtable("/Volumes/methlab/Neurometric/Antisaccades/ALL/old.xlsx");
table_young =readtable("/Volumes/methlab/Neurometric/Antisaccades/ALL/young.xlsx");
oldIDsA = table_old{3:end,1};
oldIDsB= table_old{3:end,2};
all_oldIDs = [oldIDsA; oldIDsB];
youngIDsA = table_young{3:end,1};
youngIDsB= table_young{3:end,2};
all_youngIDs = [youngIDsA; youngIDsB];
 
%raw= '/Users/mplome/Desktop/balanced/EEG_preprocessed'; % path to preprocessed eeg files
raw = '/Volumes/methlab/Neurometric/Antisaccades/main_analysis/EEG'
etfolder='/Volumes/methlab/Neurometric/Antisaccades/main_analysis/ET'; %MO
 
 
%WINDOWS!!!!
elseif system == 0 % Windows
%ta sciezka ma byc('/Volumes/methlab/Neurometric/Antisaccades/ALL/eeglab14_1_2b/'
cd('W:\Neurometric\Antisaccades\ALL\eeglab14_1_2b')
addpath('W:\Neurometric\Antisaccades\ALL\')
addpath('W:\Neurometric\Antisaccades\ALL\FastICA_25')
addpath('W:\Neurometric\Antisaccades\ALL\UNFOLD')

eeglab;
close
addpath('Tools') %
table_old = readtable("W:\Neurometric\Antisaccades\main_analysis\Resources\old.xlsx");
table_young =readtable("W:\Neurometric\Antisaccades\main_analysis\Resources\young.xlsx");
oldIDsA = table_old{3:end,1};
oldIDsB= table_old{3:end,2};
all_oldIDs = [oldIDsA; oldIDsB];
youngIDsA = table_young{3:end,1};
youngIDsB= table_young{3:end,2};
all_youngIDs = [youngIDsA; youngIDsB];
 
raw= 'W:\Neurometric\Antisaccades\main_analysis\EEG'; % path to preprocessed eeg files
etfolder='W:\Neurometric\Antisaccades\main_analysis\ET'; %MO
 
end
 
d=dir(raw) %what folders are in there (each folder = one subject)
 
d(1:3)=[] % get rid of the . and .. folders as well as .DS_Store on mac
 
 
OLD_OR_YOUNG = {'old', 'yng'};
data_pro_old = {};
data_pro_yng = {};
data_anti_old = {};
data_anti_yng = {};


data_pro_mpe_old = {};
data_pro_mpe_yng = {};
data_anti_mpe_old = {};
data_anti_mpe_yng = {};

 
suffix_pro_mpe_old = {};
suffix_anti_mpe_old = {};
suffix_pro_mpe_yng = {};
suffix_anti_mpe_yng = {};

all_subjects ={};
old_subjects = {};
young_subjects = {};

%% 
for i=1:3 %length(d) %loop over all subjects
    if d(i).isdir
        subjectfolder=dir([d(i).folder filesep d(i).name  ]);
        
        deleteindex=[];
        for ii=1:length(subjectfolder)
            if not(endsWith(subjectfolder(ii).name, '_EEG.mat')) || startsWith(subjectfolder(ii).name,'bip') || startsWith(subjectfolder(ii).name,'red')
                deleteindex(end+1)=ii;
            end
        end
        subjectfolder(deleteindex)=[];
        FullEEG=[];
        names ={};
        for kkk=1:length(subjectfolder)
        names{kkk,1} =    subjectfolder(kkk).name;
        end
        
        
        for ii=1:length(subjectfolder)
           clear ind_t ind
            for kk = 1:length(subjectfolder)
                ind(kk) = ~isempty(strfind(names{kk},['AS',num2str(ii)]));
            end
            ind_t = find(ind);
            if isempty(ind_t)
                continue;
            end

            load ([subjectfolder(ind_t).folder filesep subjectfolder(ind_t).name]) % gets loaded as EEG
          %  fileindex=subjectfolder(ii).name(end-8) %here you need to find the index from thefile (end-someting) indexing
            etfile=  [etfolder filesep d(i).name filesep d(i).name '_AS' num2str(ii) '_ET.mat'] %define string of the complete path to the matching ET file.
             
            %EEG = pop_reref(EEG,[47 83])
            
            EEG = pop_eegfiltnew(EEG,[],30);

            EEG = pop_reref(EEG,[]) %tu jest wazna sprawa, to reref
            
            %merge ET into EEG
            ev1=94 %first trigger of eeg and ET file
            ev2=50 % end trigger in eeg and ET file
            EEG=pop_importeyetracker(EEG, etfile,[ev1 ev2], [1:4], {'TIME' 'L_GAZE_X' 'L_GAZE_Y' 'L_AREA'},1,1,0,0,4);
           
            %% change triggers here
            
       % countblocks = 0;
        previous = '';
        for e = 1:length(EEG.event)
           % if strcmp(EEG.event(e).type,'94  ')
           %     countblocks = countblocks+1;
           % end
           
           if ~isempty(strfind(subjectfolder(ind_t).name,'AS2')) || ~isempty(strfind(subjectfolder(ind_t).name,'AS3')) || ~isempty(strfind(subjectfolder(ind_t).name,'AS4')) 
           % if countblocks == 2 || countblocks == 3 || countblocks == 4 % antisaccade blocks
                if strcmp(EEG.event(e).type,'10  ') % change 10 to 12 for AS
                    EEG.event(e).type = '12  ';
                elseif strcmp(EEG.event(e).type,'11  ')
                    EEG.event(e).type = '13  '; % change 11 to 13 for AS
               
                end
                if strcmp(EEG.event(e).type,'40  ') 
                    EEG.event(e).type = '41  ';
                end
            end

            if strcmp(EEG.event(e).type, 'L_saccade')
                if strcmp(previous, '10  ')
                    EEG.event(e).type = 'L_saccade_10';%pro left
                elseif strcmp(previous, '11  ')
                    EEG.event(e).type = 'L_saccade_11';%pro right
                elseif strcmp(previous, '12  ')
                    EEG.event(e).type = 'L_saccade_12';%anti left
                elseif strcmp(previous, '13  ')
                    EEG.event(e).type = 'L_saccade_13'; %anti right
                end             
            end
            
            if ~strcmp(EEG.event(e).type, 'L_fixation') ...
                    && ~strcmp(EEG.event(e).type, 'L_blink')
                previous = EEG.event(e).type;
            end
        end
            
            
            
            
            
            
            
            
            %sprawdz pupil size, w razie czego zmien
            %zmienilem ii = 2 na isempty.fu;lleeg, bo bylo zakladane ze przy ii = 2 fulleeg juz nie jest puste zakladamy
            if isempty(FullEEG)
                FullEEG=EEG;
            else
                FullEEG=pop_mergeset(FullEEG,EEG);
            end
        end

        if isempty(FullEEG)
            continue
        end
        
            
        
        

        
        id=d(i).name ;
        clear young old
        young=    any(contains(all_youngIDs,id));
        old =     any(contains(all_oldIDs,id));
        
        all_subjects{end+1,1} =id;
        if young == 0
        old_subjects{end+1,1} = id;
        elseif young == 1
        young_subjects{end+1,1} = id;  
        end
        %young means 1, old means 0
        all_ages(i) =    young;
        all_eeg{i} = FullEEG;
        
      
      
        
       %% COSTRAP 
       if use_costrap == 1
       Elecs2Use = [17 12 7]; % equivalent to Fp1 Fp2 E15
       TSP_per_second = 3;
       MinTSP_Freq = 25;
     % MARTYNA ELECTRODE = [104 9 19]; - frontal eye fields 
       FullEEG = costrap(FullEEG,Elecs2Use,TSP_per_second,MinTSP_Freq,0);

       end
%      
        %% UNFOLD
        if use_unfold == 1

        FullEEG_noeye = FullEEG;
        FullEEG_noeye.chanlocs = FullEEG.chanlocs(1:105);
        FullEEG_noeye.data = FullEEG.data(1:105, :);
        FullEEG_noeye.nbchan = 105;
        
        [ufresult,EEG_res,EEG_model] = unfold_martyna(FullEEG_noeye);
 
        EEG_res.data = [EEG_res.data; FullEEG.data(106:109, :)];
        FullEEG = EEG_res;

        end
     

        
        
        %% epoching
        
        BEGIN = -0.3;
        END = 0.5;
        try
        catch
            invalid(i) = true;
            fprintf('my_pop_epoch() failed for %s\n', id);
           % continue;    
        end
        %teraz laczymy 10 i 11
        data_pro_mpe{i} = EpochLatencyOf(...
            my_pop_epoch_vol2(FullEEG, '10  ', 'L_saccade_10'), ...
            'L_saccade_10');
        % mpev = EpochLatencyOf(mpev, 'L_saccade_10')
        % dziala na zwyklym pop_epoch
        data_pro_mpe{i} = ClipToBoundsExtra(data_pro_mpe{i}, -500, 100);
        clear tmpdata 
        tmpdata = EpochLatencyOf(...
            my_pop_epoch_vol2(FullEEG, '11  ', 'L_saccade_11'), ...
            'L_saccade_11');
        tmpdata = ClipToBoundsExtra(tmpdata, -100, 500);
        
        data_pro_mpe{i} = MergeSets(data_pro_mpe{i}, tmpdata)

        
        %teraz laczymy 12 i 13
        data_anti_mpe{i} = EpochLatencyOf(...
            my_pop_epoch_vol2(FullEEG, '12  ', 'L_saccade_12'), ...
            'L_saccade_12');
        data_anti_mpe{i} = ClipToBoundsExtra(data_anti_mpe{i}, -100, 500);
        clear tmpdata 
        tmpdata = EpochLatencyOf(...
            my_pop_epoch_vol2(FullEEG, '13  ', 'L_saccade_13'), ...
            'L_saccade_13'); 
        tmpdata = ClipToBoundsExtra(tmpdata, -500, 100);
                
        data_anti_mpe{i} = MergeSets(data_anti_mpe{i}, tmpdata);

    

        
        invalid(i) = false;
        
        
      %dividing - old/young
        eval(['data_pro_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_mpe{i}']);
        eval(['data_anti_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_mpe{i}']);
        
        
%         eval(['suffix_pro_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = suffix_pro_mpe{i}']);
%         eval(['suffix_anti_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = suffix_anti_mpe{i}']);         
%         
    end
        
end

%save('data_mpe_all_subjects.mat','data_anti_mpe','data_anti_mpe_old','data_anti_mpe_yng','data_pro_mpe_old','data_pro_mpe_yng','data_pro_mpe','all_subjects','young_subjects','old_subjects','-v7.3')

%save('suffix_mpe_all_subjects.mat','suffix_anti_mpe','suffix_anti_mpe_old','suffix_anti_mpe_yng','suffix_pro_mpe_old','suffix_pro_mpe_yng','suffix_pro_mpe','all_subjects','young_subjects','old_subjects','-v7.3')


 

%%
data_pro_mpe_old_100_150 = {};
data_pro_mpe_old_150_250 = {};
data_pro_mpe_old_250_350 = {};
data_pro_mpe_old_350_450 = {};
data_pro_mpe_old_450_500 = {};


data_anti_mpe_old_100_150 = {};
data_anti_mpe_old_150_250 = {};
data_anti_mpe_old_250_350 = {};
data_anti_mpe_old_350_450 = {};
data_anti_mpe_old_450_500 = {};


data_pro_mpe_yng_100_150 = {};
data_pro_mpe_yng_150_250 = {};
data_pro_mpe_yng_250_350 = {};
data_pro_mpe_yng_350_450 = {};
data_pro_mpe_yng_450_500 = {};


data_anti_mpe_yng_100_150 = {};
data_anti_mpe_yng_150_250 = {};
data_anti_mpe_yng_250_350 = {};
data_anti_mpe_yng_350_450 = {};
data_anti_mpe_yng_450_500 = {};

for i=1:length(data_anti_mpe_old)
for j=1:length(data_anti_mpe_old{i}.epoch)
    sz = data_anti_mpe_old{i}.epoch(j).latency
    if sz > 450
        data_anti_mpe_old_450_500{end+1} = data_anti_mpe_old{i}.cell_data{j};
    elseif sz > 350
        data_anti_mpe_old_350_450{end+1} = data_anti_mpe_old{i}.cell_data{j};
    elseif sz > 250
        data_anti_mpe_old_250_350{end+1} = data_anti_mpe_old{i}.cell_data{j};
    elseif sz > 150
        data_anti_mpe_old_150_250{end+1} = data_anti_mpe_old{i}.cell_data{j};
    else
        data_anti_mpe_old_100_150{end+1} = data_anti_mpe_old{i}.cell_data{j};
    end
end
end

for i=1:length(data_pro_mpe_old)
for j=1:length(data_pro_mpe_old{i}.epoch)
    sz = data_pro_mpe_old{i}.epoch(j).latency
    if sz > 450
        data_pro_mpe_old_450_500{end+1} = data_pro_mpe_old{i}.cell_data{j};
    elseif sz > 350
        data_pro_mpe_old_350_450{end+1} = data_pro_mpe_old{i}.cell_data{j};
    elseif sz > 250
        data_pro_mpe_old_250_350{end+1} = data_pro_mpe_old{i}.cell_data{j};
    elseif sz > 150
        data_pro_mpe_old_150_250{end+1} = data_pro_mpe_old{i}.cell_data{j};
    else
        data_pro_mpe_old_100_150{end+1} = data_pro_mpe_old{i}.cell_data{j};
    end
end
end

for i=1:length(data_anti_mpe_yng)
for j=1:length(data_anti_mpe_yng{i}.epoch)
    sz = data_anti_mpe_yng{i}.epoch(j).latency
    if sz > 450
        data_anti_mpe_yng_450_500{end+1} = data_anti_mpe_yng{i}.cell_data{j};
    elseif sz > 350
        data_anti_mpe_yng_350_450{end+1} = data_anti_mpe_yng{i}.cell_data{j};
    elseif sz > 250
        data_anti_mpe_yng_250_350{end+1} = data_anti_mpe_yng{i}.cell_data{j};
    elseif sz > 150
        data_anti_mpe_yng_150_250{end+1} = data_anti_mpe_yng{i}.cell_data{j};
    else
        data_anti_mpe_yng_100_150{end+1} = data_anti_mpe_yng{i}.cell_data{j};
    end
end
end

for i=1:length(data_pro_mpe_yng)
for j=1:length(data_pro_mpe_yng{i}.epoch)
    sz = data_pro_mpe_yng{i}.epoch(j).latency
    if sz > 450
        data_pro_mpe_yng_450_500{end+1} = data_pro_mpe_yng{i}.cell_data{j};
    elseif sz > 350
        data_pro_mpe_yng_350_450{end+1} = data_pro_mpe_yng{i}.cell_data{j};
    elseif sz > 250
        data_pro_mpe_yng_250_350{end+1} = data_pro_mpe_yng{i}.cell_data{j};
    elseif sz > 150
        data_pro_mpe_yng_150_250{end+1} = data_pro_mpe_yng{i}.cell_data{j};
    else
        data_pro_mpe_yng_100_150{end+1} = data_pro_mpe_yng{i}.cell_data{j};
    end
end
end

%%

%%%%%%%%%%%%%%%%%%%%FUNCTIOMS
function [X] = MergeSets(X, varargin)
    %last 50
    varmat = cell2mat(varargin);
    X.event = horzcat(X.event, varmat.event);
    X.trials = sum(vertcat(X.trials, varmat.trials));
    X.xmin = min(vertcat(X.xmin, varmat.xmin));
    X.xmax = max(vertcat(X.xmax, varmat.xmax));
    X.cell_data = cat(2, X.cell_data, varmat.cell_data);
    X.epoch = horzcat(X.epoch, varmat.epoch);
end

%here sorting
function [dataset_out] = SortByLatencyOf(dataset, event)
    dataset_out = EpochLatencyOf(dataset, event);
    [B,I] = sort([dataset_out.epoch.latency]);
    dataset_out.epoch_sorted = dataset_out.epoch(I);
    dataset_out.data_sorted = dataset_out.data(:,:,I);
end

function [dataset] = EpochLatencyOf(dataset, event)
    for i=1:length(dataset.epoch)
        sacc_idx = strcmp(dataset.epoch(i).eventtype, event);
        if sum(sacc_idx) == 1
            dataset.epoch(i).latency = dataset.epoch(i).eventlatency{sacc_idx};
        else
            dataset.epoch(i).latency = -1;
        end
    end
end




