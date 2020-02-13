%% NEw! 25.01

%%

clear all

system = 0 %0 (1 = MAC, 0 = Windows)
use_costrap = 1 %1 = do
use_unfold = 1 % 1 = do unfold 0 = don't

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
data_pro_right_old = {};
data_pro_right_yng = {};
data_anti_right_old = {};
data_anti_right_yng = {};
data_pro_left_old = {};
data_pro_left_yng = {};
data_anti_left_old = {};
data_anti_left_yng = {};


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
for i=1:length(d) %loop over all subjects
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
     
%DOTAD SPOKO
        
    %%  START
    
        BEGIN = -0.75;
        END = 0.8;
        try
        %data epoching
        data_pro_right{i} = pop_epoch(FullEEG, {'11  '}, [BEGIN,END]);
        data_pro_left{i} = pop_epoch(FullEEG, {'10  '}, [BEGIN,END]);
        data_anti_right{i} = pop_epoch(FullEEG, {'13  '}, [BEGIN,END]);
        data_anti_left{i} = pop_epoch(FullEEG, {'12  '}, [BEGIN,END]);
  
        data_pro_right{i} = pop_rmbase(data_pro_right{i},1000*[BEGIN,0]);
        data_pro_left{i} = pop_rmbase(data_pro_left{i},1000*[BEGIN,0]);
        data_anti_right{i} = pop_rmbase(data_anti_right{i},1000*[BEGIN,0]);
        data_anti_left{i} = pop_rmbase(data_anti_left{i},1000*[BEGIN,0]);
        
%          I'm interested only in saccades with latencies [100;500]ms ale nie
%         dziala
        data_pro_right{i} = FilterByEventLatency(data_pro_right{i}, 'L_saccade_11', 50, 800);
        data_pro_left{i} = FilterByEventLatency(data_pro_left{i}, 'L_saccade_10', 50, 800);
        data_anti_right{i} = FilterByEventLatency(data_anti_right{i}, 'L_saccade_13', 50, 800);
        data_anti_left{i} = FilterByEventLatency(data_anti_left{i}, 'L_saccade_12', 50, 800);
        
        catch
            invalid(i) = true;
            fprintf('my_pop_epoch() failed for %s\n', id);
            continue;    
        end
        invalid(i) = false;
    
        %dividing - old/young
        eval(['data_pro_right_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_right{i}']);
        eval(['data_pro_left_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_left{i}'])
        eval(['data_anti_right_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_right{i}']);
        eval(['data_anti_left_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_left{i}']);
    end
        
end

%%
for i=1:length(data_pro_right_yng)
    [data_pro_right_yng{i}, pri] =ClipToBounds(data_pro_right_yng{i}, -100, 500);
    [data_pro_left_yng{i}, pli] = ClipToBounds(data_pro_left_yng{i}, -500, 100);
    [data_anti_right_yng{i}, ari] = ClipToBounds(data_anti_right_yng{i}, -500, 100);
    [data_anti_left_yng{i}, ali] = ClipToBounds(data_anti_left_yng{i}, -100, 500);

    data_pro_yng{i} = MergeSets(data_pro_right_yng{i}, data_pro_left_yng{i});
    data_anti_yng{i} = MergeSets(data_anti_right_yng{i}, data_anti_left_yng{i});
end

for i=1:length(data_pro_right_old)
    [data_pro_right_old{i}, pri] =ClipToBounds(data_pro_right_old{i}, -100, 500);
    [data_pro_left_old{i}, pli] = ClipToBounds(data_pro_left_old{i}, -500, 100);
    [data_anti_right_old{i}, ari] = ClipToBounds(data_anti_right_old{i}, -500, 100);
    [data_anti_left_old{i}, ali] = ClipToBounds(data_anti_left_old{i}, -100, 500);

    data_pro_old{i} = MergeSets(data_pro_right_old{i}, data_pro_left_old{i});
    data_anti_old{i} = MergeSets(data_anti_right_old{i}, data_anti_left_old{i});
end

%%
%%

ITEMS = {
    {data_pro_yng, {'L_saccade_10', 'L_saccade_11'}, 'PRO  YOUNG GROUP'}
    {data_anti_yng, {'L_saccade_12', 'L_saccade_13'}, 'ANTI   YOUNG GROUP'}
    
	{data_pro_old, {'L_saccade_10', 'L_saccade_11'}, 'PRO OLD GROUP'}
	{data_anti_old, {'L_saccade_12', 'L_saccade_13'}, 'ANTI OLD GROUP'}
	
};
I_DATA = 1;
I_TRIGGER = 2;
I_TITLE = 3;
%%



%TUTAJ NOWE KODOWANE
ELECTRODE = [19 9 104];

VAR_NAMES = {'pro_yng', 'anti_yng', 'pro_old','anti_old'};
for I=1:length(ITEMS)
  for i=1:length(ITEMS{I}{I_DATA})
    if isempty(ITEMS{I}{I_DATA}{i}.data)
        continue;
    end
    erp(i,:,:) = mean(ITEMS{I}{I_DATA}{i}.data,3);
  end
  groupaverage_erp = squeeze(mean(mean(erp(:,ELECTRODE,:),2),1));
  eval(['erp_' VAR_NAMES{I} ' = erp;']);
  eval(['groupaverage_erp_' VAR_NAMES{I} ' = groupaverage_erp;']);
  
  clear erp  groupaverage_erp;
end
nsubj_old=size(data_anti_old,2);
nsubj_yng=size(data_anti_yng,2);



%% plots
figure
subplot(2,2,1)
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'old', 'yng'})
title ('PRO')
hold off

subplot(2,2,3)
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'old', 'yng'})
title ('ANTI')
hold off


%% to wyglada spoko - pierwszy jest stary plot 
figure;
for I=1:length(ITEMS)
    subplot(2,2,I);
    %new method starts here
    X = MergeSets(ITEMS{I}{I_DATA}{:});%tutaj z tych: PODAWANEJEST EVENTS

    X = SortByLatencyOf(X, ITEMS{I}{I_TRIGGER})
    X.data_sorted;
    image_sc_data= squeeze(mean(X.data_sorted(ELECTRODE,:,:)))';
    save([ITEMS{I}{I_TITLE},'_data_rt.mat'],'X','-v7.3')

    colormap('jet');
    caxis([-25,25]);
    
    


    new_to_image  = [];
    rts = [];
    for i = 26: size(image_sc_data,1) - 25
       new_to_image(end + 1, :) =  mean(image_sc_data(i-25:i+25,:),1);
      
       rts(end + 1) = X.epoch_sorted(i).latency/2;
       
    end


    start_idx = -500*BEGIN + 1;%konwersja jednostek , baseline 
    imagesc(new_to_image(:,start_idx:end))
    hold on
    plot( rts,1:length(rts),'black','LineWidth',2)
    title (ITEMS{I}{I_TITLE})
       caxis([-3,3]);
    colorbar
    xticklabels(sprintfc("%d", xticks*2));
    set(gcf,'color','w')
    %saveas(gcf, ['W:\Neurometric\Antisaccades\results\single_trial\three_electrodes_merged\' num2str(I)  '.png'] )
    hold on
    
end
%%
figure;
for I=1:length(ITEMS)
    subplot(2,2,I);
    %new method starts here
    X = MergeSets(ITEMS{I}{I_DATA}{:});%tutaj z tych: PODAWANEJEST EVENTS

    X = SortByLatencyOf(X, ITEMS{I}{I_TRIGGER})
    X.data_sorted;
    image_sc_data= squeeze(mean(X.data_sorted(ELECTRODE,:,:)))';
    save([ITEMS{I}{I_TITLE},'_data_cue.mat'],'X','-v7.3')
    colormap('jet');
    caxis([-25,25]);
    
    
    %adding rhe que time plot

    epoch_begin = -500 * BEGIN;
    new_to_image  = [];
    rts = [];
    for i = 26: size(image_sc_data,1) - 25
       limit = epoch_begin + X.epoch_sorted(i).latency/2; %to mi robi limit x, 
       %fakt ze chce zeby na koncu bylo bariera reaction time wszutskech
       %triali  w jednym miejscu
       new_to_image(end + 1, :) =  mean(image_sc_data(i-25:i+25,limit-400:limit),1); %lecimy 800 ms do tylu
       rts(end+1) = epoch_begin - limit + 400;%dupa z nazwa, teraz to jest po que, end+1 - budowanie
       %ciagu, dodawanie do konca, bo dla kazdego triala chce miec jedna
       %liczbe
    end


%adding tge cue plot
imagesc(new_to_image(:,:))
    hold on
    plot( rts,1:length(rts),'black','LineWidth',2)
    title (ITEMS{I}{I_TITLE})
       caxis([-3,3]);
    colorbar
    xticklabels(sprintfc("%d", xticks*2));
    set(gcf,'color','w')
    hold on
    
end
    %% Func
%here 2 functions - merging+sorting
function [X] = MergeSets(X, varargin)
    varmat = cell2mat(varargin);
    X.event = horzcat(X.event, varmat.event);
    X.trials = sum(vertcat(X.trials, varmat.trials));
    X.xmin = min(vertcat(X.xmin, varmat.xmin));
    X.xmax = max(vertcat(X.xmax, varmat.xmax));
    X.data = cat(3, X.data, varmat.data); %CONCATENATING ARRAYS , 3 = DIM
    X.epoch = horzcat(X.epoch, varmat.epoch); %Concatenate arrays horizontally
end


%here sorting
function [dataset_out] = SortByLatencyOf(dataset, events)
    dataset_out = EpochLatencyOf(dataset, events);
    [B,I] = sort([dataset_out.epoch.latency]);
    dataset_out.epoch_sorted = dataset_out.epoch(I);
    dataset_out.data_sorted = dataset_out.data(:,:,I);
end

%TERAZ SPRAWDZAMY JAKIE BYLO LATENCY DANEGO EVENTU
function [dataset] = EpochLatencyOf(dataset, events)
    for i=1:length(dataset.epoch)
        for e=1:length(events)
            event = events(e);
            sacc_idx = strcmp(dataset.epoch(i).eventtype, event);%WARUNEK ZEBY SPRAWDZIC INDEXMEVENTU
            %PORIWNUJE EVENT KTOREGO SZUKAM  Z EVENTEM KTORY WYSTAPIL
            if sum(sacc_idx) == 1 
                %JAK NAJDE TO NA INDEKSIE JEST ZAPALONY BIT, - TRUE - JESLI
                %TYLKO  JEDNEN  BIT JEST ZAPALONY, TO ZNACZY ZE ZNALAZL SIE
                %EVENT. 
                break;
            end
        end
        if sum(sacc_idx) == 1
            dataset.epoch(i).latency = dataset.epoch(i).eventlatency{sacc_idx};%JESLI ZNALAZLEM INDEX
            %TO BIORE LATRENCY TEGO EVENTU
        else
            dataset.epoch(i).latency = -1; %JESLI NIE, TO MINUS 1
        end
    end
end




