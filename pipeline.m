

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
     

        
        
        %% epoching
        
        BEGIN = -0.3;
        END = 0.5;
        try
        %teraz laczymy 10 i 11
        data_pro_mpe{i} = my_pop_epoch(FullEEG, '10  ', 'L_saccade_10');
        % mpev = EpochLatencyOf(mpev, 'L_saccade_10')
        data_pro_mpe{i} = clip_tb_cell(data_pro_mpe{i}, -500, 100);
        clear tmpdata 
        tmpdata = my_pop_epoch(FullEEG, '11  ', 'L_saccade_11');
        tmpdata = clip_tb_cell(tmpdata, -100, 500);
        
        data_pro_mpe{i} = [data_pro_mpe{i} tmpdata]
        suffix_pro_mpe{i} = SuffixData(data_pro_mpe{i})

        
        %teraz laczymy 12 i 13
        data_anti_mpe{i} = my_pop_epoch(FullEEG, '12  ', 'L_saccade_12');
        data_anti_mpe{i} = clip_tb_cell(data_anti_mpe{i}, -100, 500);
        clear tmpdata 
        tmpdata = my_pop_epoch(FullEEG, '13  ', 'L_saccade_13'); 
        tmpdata = clip_tb_cell(tmpdata, -500, 100);
                
        data_anti_mpe{i} = [data_anti_mpe{i} tmpdata]
        suffix_anti_mpe{i} = SuffixData(data_anti_mpe{i})

    

        catch
            invalid(i) = true;
            fprintf('my_pop_epoch() failed for %s\n', id);
            continue;    
        end
        invalid(i) = false;
        
        
      %dividing - old/young
        eval(['data_pro_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_pro_mpe{i}']);
        eval(['data_anti_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = data_anti_mpe{i}']);
        
        
        eval(['suffix_pro_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = suffix_pro_mpe{i}']);
        eval(['suffix_anti_mpe_' OLD_OR_YOUNG{young+1} '{end+1} = suffix_anti_mpe{i}']);         
        
    end
        
end

save('data_mpe_all_subjects.mat','data_anti_mpe','data_anti_mpe_old','data_anti_mpe_yng','data_pro_mpe_old','data_pro_mpe_yng','data_pro_mpe','all_subjects','young_subjects','old_subjects','-v7.3')

save('suffix_mpe_all_subjects.mat','suffix_anti_mpe','suffix_anti_mpe_old','suffix_anti_mpe_yng','suffix_pro_mpe_old','suffix_pro_mpe_yng','suffix_pro_mpe','all_subjects','young_subjects','old_subjects','-v7.3')


%% Compute Grand Average


  ELECTRODE = [19 9 104];

VAR_NAMES = {'pro_yng',  'anti_yng', 'pro_old', 'anti_old'};
         
         
for I=1:length(VAR_NAMES)
  mpe_name = [VAR_NAMES{I}(1:end-4) '_mpe_' VAR_NAMES{I}(end-2:end)];
  eval(['curdata = suffix_' mpe_name ';'])
  if isempty(curdata)
      continue;
  end
  for i=1:length(curdata)
    if isempty(curdata{i})
        continue;
    end
    erp(i,:,:) = mean(curdata{i}(:,end-50+1:end,:),3);
  end
  
  groupaverage_erp = squeeze(mean(mean(erp(:,ELECTRODE,:),2),1));
  eval(['erp_' VAR_NAMES{I} ' = erp;']);
  eval(['groupaverage_erp_' VAR_NAMES{I} ' = groupaverage_erp;']);
  
  clear erp groupaverage_erp curdata
end

%%
%PLOTS
 
nsubj_old=size(data_anti_mpe_old,2);
nsubj_yng=size(data_anti_mpe_yng,2);
 
figure
set(gcf,'color','w')
subplot(1,2,1)
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
plot(5:37, ones(1, 37-5+1) * -0.65, 'LineWidth', 5)%type
plot(1:50, ones(1, 50) * -0.50, 'LineWidth', 5)%age
plot(21:50, ones(1, 50-21+1) * -0.35, 'LineWidth', 5)%inter
legend({'Old', 'Young'})  
[hleg] = legend('Location','southoutside','NumColumns',3,'FontSize',7)
title ('PRO, Old Group vs Young Group')
ylim([-0.7 1])
xticklabels(sprintfc("%d", xticks*2-100));
hold off
 
subplot(1,2,2)
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
plot(5:37, ones(1, 37-5+1) * -0.65, 'LineWidth', 5)%addidgs pval bars for ohbm abstract
plot(1:50, ones(1, 50) * -0.50, 'LineWidth', 5)
plot(21:50, ones(1, 50-21+1) * -0.35, 'LineWidth', 5)

legend({'Old', 'Yng'})
[hleg] = legend('Location','southoutside','NumColumns',3,'FontSize',7)

title ('ANTI, Old Group vs Young Group')
ylim([-0.7 1])
xticklabels(sprintfc("%d", xticks*2-100));
hold off

%sgtitle('Old Group vs Young Group')

figure
set(gcf,'color','w')    
subplot(1,2,1)
shadedErrorBar(1:size(erp_pro_yng,3), groupaverage_erp_pro_yng,std(mean(erp_pro_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_yng,3), groupaverage_erp_anti_yng,std(mean(erp_anti_yng(:,[19,9,104],:),2))/sqrt(nsubj_yng) ,'lineprops', 'b');
legend({'Pro ', 'Anti '})
[hleg] = legend('Location','southoutside','NumColumns',3,'FontSize',7)

title ('PRO  vs ANTI  (Young Group)')

ylim([-0.5 1])
xticklabels(sprintfc("%d", xticks*2-100));
hold off
 
subplot(1,2,2)
shadedErrorBar(1:size(erp_pro_old,3), groupaverage_erp_pro_old,std(mean(erp_pro_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'r');
hold on
shadedErrorBar(1:size(erp_anti_old,3), groupaverage_erp_anti_old,std(mean(erp_anti_old(:,[19,9,104],:),2))/sqrt(nsubj_old) ,'lineprops', 'b');
legend({'Pro ', 'Anti'})
[hleg] = legend('Location','southoutside','NumColumns',3,'FontSize',7)

title ('PRO vs ANTI (Old Group)')
ylim([-0.5 1])
xticklabels(sprintfc("%d", xticks*2-100));
hold off
 
%% STATISTIC PREPARATION _TABLES
ELECTRODE = [104 9 19];
 
 
c_tables={}; %cell array of tables for each timepoint
for i=1:50
    id=0;
    type=0;
    age=0;
    trialdat=0;
    tbl=table(id,type,age,trialdat,'VariableNames',{'id','type','age','trialdat'});
    c_tables{i}=tbl;
end
 
%create matrix
for sj=1:length(suffix_pro_mpe_old)
    m=suffix_pro_mpe_old{sj};
    m_suffix_proold=squeeze(mean(m(ELECTRODE,:,:)));
    for tp=1:50 %size(m_suffix_proleftold, 1)
        currtbl=c_tables{tp};
        id=sj;
        type=0;%0 mean pro, 1 measn anti
        age=0;%1=old? check coding
        for trial=1:size(m_suffix_proold,2);
            trialdat=m_suffix_proold(end-51+tp,trial);
            tmptbl=table(id,type,age,trialdat,'VariableNames',{'id','type','age','trialdat'});
            c_tables{tp}=[c_tables{tp}; tmptbl];
        end
    end
    
    
   
    m=suffix_anti_mpe_old{sj};
    m_suffix_antiold=squeeze(mean(m(ELECTRODE,:,:)));
     for tp=1:50 %size(m_suffix_antirightold, 1)
        currtbl=c_tables{tp};
        id=sj;
        type=1;%0 mean pro, 1 means anti
        age=0;%1=old? check coding
        for trial=1:size(m_suffix_antiold,2);
            trialdat=m_suffix_antiold(end-51+tp,trial);
            tmptbl=table(id,type,age,trialdat,'VariableNames',{'id','type','age','trialdat'});
            c_tables{tp}=[c_tables{tp}; tmptbl];
        end
     end
end

for sj=1:length(suffix_pro_mpe_yng)
    %new . YNG
 
    m=suffix_anti_mpe_yng{sj};
    m_suffix_antiyng=squeeze(mean(m(ELECTRODE,:,:)));
     for tp=1: 50 %size(m_suffix_antileftyng, 1)
        currtbl=c_tables{tp};
        id=sj;
        type=1;%1 to anti
        age=1;%0=old? check coding
        for trial=1:size(m_suffix_antiyng,2);
            trialdat=m_suffix_antiyng(end-51+tp,trial);
            tmptbl=table(id,type,age,trialdat,'VariableNames',{'id','type','age','trialdat'});
            c_tables{tp}=[c_tables{tp}; tmptbl];
        end
     end
    
    
      %new . YNG
 
    
    m=suffix_pro_mpe_yng{sj};
    m_suffix_proyng=squeeze(mean(m(ELECTRODE,:,:)));
     for tp=1: 50 %size(m_suffix_proleftyng, 1)
        currtbl=c_tables{tp};
        id=sj;
        type=0;%1 to anti
        age=1;%1=old? check coding
        for trial=1:size(m_suffix_proyng,2);
            trialdat=m_suffix_proyng(end-51+tp,trial);
            tmptbl=table(id,type,age,trialdat,'VariableNames',{'id','type','age','trialdat'});
            c_tables{tp}=[c_tables{tp}; tmptbl];
        end
    end
    
 
 
end
 
 
 
%% STATISTICS, MIXED MODELS
 
%tutaj problem - must be of data type double.
for i = 1:50
    current_tbl = c_tables{i};
    current_tbl.age = categorical(current_tbl.age);
    current_tbl.type = categorical(current_tbl.type);
    current_tbl.trialdat = double(current_tbl.trialdat);
    glme{i} = fitglme( current_tbl,'trialdat~ 1 +type*age + (1|id)');
end


glme{50}


 clear age_ts age_ps type_ps type_ts direction_ps 
 for i=1:50
     currModel=glme{i};
    age_ts(i)=currModel.Coefficients{3,4};
    type_ts(i)=currModel.Coefficients{2,4};
    type_age_ts(i)= currModel.Coefficients{4,4}
    age_ps(i)=currModel.Coefficients{3,6};
    type_ps(i)=currModel.Coefficients{2,6};
   type_age_ps(i)= currModel.Coefficients{4,6};
 end
 
figure;
plot(age_ts)
hold on
plot(type_ts)
plot(type_age_ts)
legend({'age','type', 'inter'});title('T values of models on each timepoint')
 
 
figure;
plot(age_ps)
hold on
plot(type_ps)
plot(type_age_ps)
legend({'age','type', 'inter'});title('P values of models on each timepoint')
set(gcf,'color','w')


 
 
 