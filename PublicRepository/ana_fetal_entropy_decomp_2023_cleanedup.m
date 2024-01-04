% Script for running the fetal entropy decomposition 

%%% The script begins:
%     Clear existing variables in the workspace.
%     Load data from 'Data_traces.mat'.
%     Set the random number generator seed.
%     Define parameters: the number of surrogates (nsur) and sampling rate (fs).
%     Calculate quantiles to separate the data into bottom and top thirds based on the gestational age 'GA' column of the 'Datatable'.

clearvars
load Data_traces.mat
rng(2345894)
nsur = 250; % number of surrogates
fs = 610.3516; % sampling rate
as1 = quantile(Datatable.GA,1/3); % bottom 33% of ages
as2 = quantile(Datatable.GA,2/3); % top 33% of ages

% Suppresses warnings (otherwise we get lots of unnecessary warnings) 
warning('off', 'all')

% Identify fetal subjects with both young and old data by finding the intersection of 'young' and 'old' subjects.
young = Datatable(Datatable.GA < as1,:).ID;
old = Datatable(Datatable.GA > as2,:).ID;
use = intersect(young,old); % subjects with young and old data

SIDs = unique(Datatable.ID);
tar_hdecomp = nan(length(SIDs),1);
ssss = nan(length(SIDs),1);
sssS = nan(length(SIDs),1);
sssd = nan(length(SIDs),1);
sssD = nan(length(SIDs),1);

% Create new variables and organizes the data in 'DT2'.
recordings = unique([Datatable.ID Datatable.GA],'rows');
DT2 = Datatable(:,1:3);
DT2.recording = nan(size(DT2,1),1);
DT2.target = nan(size(DT2,1),1);


for irow = 1:size(Datatable,1)
    tmp = [Datatable.ID(irow) Datatable.GA(irow)];
    idx = find(ismember(recordings(:,1),tmp(1)) & ismember(recordings(:,2),tmp(2))); % row that contain the same recording
    DT2.recording(irow) = idx;
    if ismember(Datatable.ID(irow),use)
        ages = DT2.GA(DT2.ID==Datatable.ID(irow));
        if DT2.GA(irow) == max(ages) || DT2.GA(irow) == min(ages) % we want to the oldest or youngest recording
            DT2.target(irow) = true;
        else
            DT2.target(irow) = false;
        end
    else
        DT2.target(irow) = false;
    end
end

[~,NDX] = unique([DT2.ID DT2.GA DT2.recording DT2.target],'rows');
DT3 = DT2(NDX,:);
DT3.BlockA = nan(size(DT3,1),1);
DT3.BlockB = nan(size(DT3,1),1);
DT3.main = ones(size(DT3,1),1);

for irow = 1:size(DT3,1)
    ri = DT3.recording(irow);
    blocks = DT2(DT2.recording==ri,:).Blocklabel;
    if any(ismember(blocks,'A (ssss)','rows'))
        DT3.BlockA(irow) = true;
    else
        DT3.BlockA(irow) = false;
    end
    if any(ismember(blocks,'B (sssd)','rows'))
        DT3.BlockB(irow) = true;
    else
        DT3.BlockB(irow) = false;
    end
end

DT3.Blocklabel = [];
DT3 = DT3(:,[3 1 2 7 4 5 6]);


%% Get ready for the entropy decomp -- necessary preparation

%     Data Formatting:
%         Convert logical values in specific fields of the structure DT3 to string representations.
% 
%     Data Export:
%         Write the contents of the modified DT3 structure to a CSV file named 'SupplementalTableOfFetalSubjects.csv'.
%         Load a MATLAB file named 'newbornIDs.mat'.
% 
%     Data Modification:
%         Create a new logical column named 'NewbornData' in DT3 based on whether the recording is associated with newborns.
% 
%     Data Export (Latex):
%         Convert the resulting table to a LaTeX file named 'SupplementalTableOfFetalSubjects'.
% 
%     Data Filtering:
%         Separate data into two tables (Ty and To) based on gestational age criteria (as1 and as2).
% 
%     Iterative Processing:
%         Iterate through a list of subjects (use) to find the earliest (idxy) and latest (idxo) gestational ages in the tables.
% 
%     Result Storage:
%         Initialize cell arrays (pe_decomps, se_decomps, mse_decomps, lzc_decomps, and ctw_decomps) to store decomposition results for each subject.
% 
%     Execution Time Tracking:
%         Begin a timer (tic) to track the execution time of the code.

DT3.main = logical2str(DT3.main);


DT3.target = logical2str(DT3.target);
DT3.BlockA = logical2str(DT3.BlockA);
DT3.BlockB = logical2str(DT3.BlockB);
writetable(DT3,'SupplementalTableOfFetalSubjects.csv')
load newbornIDs.mat
DT3.NewbornData = logical2str(ismember(DT3.recording,newborns));
table2latex(DT3,'SupplementalTableOfFetalSubjects')

Ty = Datatable(Datatable.GA < as1,:);
To = Datatable(Datatable.GA > as2,:);

idxy = nan(1,length(use));
idxo = nan(1,length(use));

for isub = 1:length(use)
    SID = use(isub);

    rows = find(Ty.ID == SID);
    idxy(isub) = min(Ty.GA(rows)); % take the earliest gestational age for that subject

    rows = find(To.ID == SID);
    idxo(isub) = max(To.GA(rows));
end

header = {'ssss','sssS','sssd','sssD'};
pe_decomps = cell(4,length(use));
se_decomps = cell(4,length(use));
mse_decomps = cell(4,length(use));
lzc_decomps = cell(4,length(use));
ctw_decomps = cell(4,length(use));

tic


%% modified MSE (mMSE)

% (MSE) Decomposition:

%     Iterative Loop:
%         Iterate through each subject (use) and each condition (header).
%         Different cases based on the condition, e.g., 'ssss', 'sssS', 'sssd', 'sssD'.
%         Find corresponding rows in Ty and To for each condition, subject, and age.
%         If data is available for the condition, perform MSE decomposition using mMSEDecomposition_fMEG.
%         Store the output in the mse_decomps cell array.
% 
%     Print Progress:
%         Display a progress message for each subject.
%         Output the elapsed time.
% 
%     Save Results:
%         Save the decomposition results and headers in a MAT file named 'MSEDecompfMEG2023_Nsur=%i'.


for isub = 1:length(use)
    for icond = 1:length(header)
        switch header{icond}
            case 'ssss'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mMSEDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                mse_decomps{strcmp(header,'ssss'),isub} = out;
            case 'sssS'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mMSEDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                mse_decomps{strcmp(header,'sssS'),isub} = out;
            case 'sssd'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mMSEDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                mse_decomps{strcmp(header,'sssd'),isub} = out;

            case 'sssD'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mMSEDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                mse_decomps{strcmp(header,'sssD'),isub} = out;
            otherwise
                error('can''t read the condition')
        end
    end

    fprintf('Finished subject #%i (ID = %i)\n',isub,use(isub))
    toc
end
save(sprintf('MSEDecompfMEG2023_Nsur=%i',nsur),'mse_decomps','header')


%%% Now do the same below for mSampEn

% modified Sample Entropy
for isub = 1:length(use)
    for icond = 1:length(header)
        switch header{icond}
            case 'ssss'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mSampEnDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                se_decomps{strcmp(header,'ssss'),isub} = out;
            case 'sssS'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mSampEnDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                se_decomps{strcmp(header,'sssS'),isub} = out
            case 'sssd'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mSampEnDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                se_decomps{strcmp(header,'sssd'),isub} = out

            case 'sssD'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = mSampEnDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                se_decomps{strcmp(header,'sssD'),isub} = out;
            otherwise
                error('can''t read the condition')
        end
    end

    fprintf('Finished subject #%i (ID = %i)\n',isub,use(isub))
    toc
end
save(sprintf('mSampEnDecompfMEG2023_Nsur=%i',nsur),'se_decomps','header')

%%% Now do the same below for PermEn

% Permutation Entropy
for isub = 1:length(use)
    for icond = 1:length(header)
        switch header{icond}
            case 'ssss'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = PermEnDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                pe_decomps{strcmp(header,'ssss'),isub} = out;
            case 'sssS'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = PermEnDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                pe_decomps{strcmp(header,'sssS'),isub} = out;
            case 'sssd'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = PermEnDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                pe_decomps{strcmp(header,'sssd'),isub} = out;

            case 'sssD'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = PermEnDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                pe_decomps{strcmp(header,'sssD'),isub} = out;
            otherwise
                error('can''t read the condition')
        end
    end

    fprintf('Finished subject #%i (ID = %i)\n',isub,use(isub))
    toc
end
save(sprintf('PermEnDecompfMEG2023_Nsur=%i',nsur),'pe_decomps','header')


%%% Now do the same below for LZC

% Lempel-Ziv
for isub = 1:length(use)
    for icond = 1:length(header)
        switch header{icond}
            case 'ssss'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = LZDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                lzc_decomps{strcmp(header,'ssss'),isub} = out;
            case 'sssS'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = LZDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                lzc_decomps{strcmp(header,'sssS'),isub} = out;
            case 'sssd'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = LZDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                lzc_decomps{strcmp(header,'sssd'),isub} = out;

            case 'sssD'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = LZDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                lzc_decomps{strcmp(header,'sssD'),isub} = out;
            otherwise
                error('can''t read the condition')
        end
        fprintf('/')
    end

    fprintf('Finished subject #%i (ID = %i)\n',isub,use(isub))
    toc
end
save(sprintf('LZCDecompfMEG2023_Nsur=%i',nsur),'lzc_decomps','header')


%%% Now do the same below for CTW

% Context-tree weighting
for isub = 1:length(use)
    for icond = 1:length(header)
        switch header{icond}
            case 'ssss'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = CTWDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                ctw_decomps{strcmp(header,'ssss'),isub} = out;
            case 'sssS'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = CTWDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                ctw_decomps{strcmp(header,'sssS'),isub} = out;
            case 'sssd'
                % ssss block label for that subject and age
                rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='B (sssd)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='B (sssd)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = CTWDecomposition_fMEG(Ty.global_standards_normalized(rowy), use(isub), ...
                    To.global_standards_normalized(rowo), use(isub), nsur, fs);

                ctw_decomps{strcmp(header,'sssd'),isub} = out;

            case 'sssD'
                % ssss block label for that subject and age
                 rowy = find(Ty.ID == use(isub) & Ty.GA == idxy(isub) & all(Ty.Blocklabel=='A (ssss)',2));
                rowo = find(To.ID == use(isub) & To.GA == idxo(isub) & all(To.Blocklabel=='A (ssss)',2));
                if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                % young (high entropy) minus old (low entropy)
                out = CTWDecomposition_fMEG(Ty.global_deviants_normalized(rowy), use(isub), ...
                    To.global_deviants_normalized(rowo), use(isub), nsur, fs);

                ctw_decomps{strcmp(header,'sssD'),isub} = out;
            otherwise
                error('can''t read the condition')
                
        end
        fprintf('/')
    end

    fprintf('Finished subject #%i (ID = %i)\n',isub,use(isub))
    toc
end
save(sprintf('CTWDecompfMEG2023_Nsur=%i',nsur),'ctw_decomps','header')


%% Starting from this line? load the data here (if you already ran the decomp earlier)

% If the decomposition results are not already present, load them from corresponding MAT files.

if ~exist('ctw_decomps','var')
    load('CTWDecompfMEG2023_Nsur=250.mat')
    load('LZCDecompfMEG2023_Nsur=250.mat')
    load('MSEDecompfMEG2023_Nsur=250.mat')
    load('mSampEnDecompfMEG2023_Nsur=250.mat')
    load('PermEnDecompfMEG2023_Nsur=250.mat')
end




%%%  Analyze the results %%%

%% Extract output, mSampEn

%     Extract mSampEn components (phase, spec, interact) for each condition and subject.
%     Combine the data into a single vector.
%     Fit linear mixed-effects models (fitlme) for both the model with and without the component variable.
%     Perform a log-likelihood ratio test (compare) to compare the two models.
%     Save the results in variables (LRse and Pse).
%     This code repeats below for each entropy measure.

seP = nan(length(se_decomps),length(header))';
seS = nan(length(se_decomps),length(header))';
seX = nan(length(se_decomps),length(header))';

for icond = 1:length(header)
    for isub = 1:length(se_decomps)
        if ~isempty(se_decomps{icond,isub}) %&& ~strcmp(subs1{isub},bl)
            seP(icond,isub) = se_decomps{icond,isub}.phase;
            seS(icond,isub) = se_decomps{icond,isub}.spec;
            seX(icond,isub) = se_decomps{icond,isub}.interact;
        end
    end
end

sePX = seP + seX; % add phase + interaction

% Combine the data into a single vector
data = [seS; sePX];
cond = repmat([1:length(header)]',2,length(se_decomps));
comp = [zeros(length(header), length(se_decomps)); ones(length(header), length(se_decomps))];
ID   = repmat(1:12,8,1);
data = data(:);
cond = cond(:);
comp = comp(:);
ID = categorical(ID(:));

X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRse = out.LRStat(out.Model=='mdl2');
Pse = out.pValue(out.Model=='mdl2');


%% Extract output, mMSE

mseP = nan(length(mse_decomps),length(header))';
mseS = nan(length(mse_decomps),length(header))';
mseX = nan(length(mse_decomps),length(header))';

for icond = 1:length(header)
    for isub = 1:length(mse_decomps)
        if ~isempty(mse_decomps{icond,isub}) %&& ~strcmp(subs1{isub},bl)
            mseP(icond,isub) = mse_decomps{icond,isub}.phase;
            mseS(icond,isub) = mse_decomps{icond,isub}.spec;
            mseX(icond,isub) = mse_decomps{icond,isub}.interact;
        end
    end
end

msePX = mseP + mseX; % add phamse + interaction

% Combine the data into a single vector
data = [mseS; msePX];
cond = repmat([1:length(header)]',2,length(mse_decomps));
comp = [zeros(length(header), length(mse_decomps)); ones(length(header), length(mse_decomps))];
ID   = repmat(1:12,8,1);
data = data(:);
cond = cond(:);
comp = comp(:);
ID = categorical(ID(:));

X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRmse = out.LRStat(out.Model=='mdl2');
Pmse = out.pValue(out.Model=='mdl2');

%% Extract output, LZC

lzcP = nan(length(lzc_decomps),length(header))';
lzcS = nan(length(lzc_decomps),length(header))';
lzcX = nan(length(lzc_decomps),length(header))';

for icond = 1:length(header)
    for isub = 1:length(lzc_decomps)
        if ~isempty(lzc_decomps{icond,isub}) %&& ~strcmp(subs1{isub},bl)
            lzcP(icond,isub) = lzc_decomps{icond,isub}.phase;
            lzcS(icond,isub) = lzc_decomps{icond,isub}.spec;
            lzcX(icond,isub) = lzc_decomps{icond,isub}.interact;
        end
    end
end

lzcPX = lzcP + lzcX; % add phalzc + interaction

% Combine the data into a single vector
data = [lzcS; lzcPX];
cond = repmat([1:length(header)]',2,length(lzc_decomps));
comp = [zeros(length(header), length(lzc_decomps)); ones(length(header), length(lzc_decomps))];
ID   = repmat(1:12,8,1);
data = data(:);
cond = cond(:);
comp = comp(:);
ID = categorical(ID(:));

X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRlzc = out.LRStat(out.Model=='mdl2');
Plzc = out.pValue(out.Model=='mdl2');


%% Extract output, CTW

ctwP = nan(length(ctw_decomps),length(header))';
ctwS = nan(length(ctw_decomps),length(header))';
ctwX = nan(length(ctw_decomps),length(header))';

for icond = 1:length(header)
    for isub = 1:length(ctw_decomps)
        if ~isempty(ctw_decomps{icond,isub}) %&& ~strcmp(subs1{isub},bl)
            ctwP(icond,isub) = ctw_decomps{icond,isub}.phase;
            ctwS(icond,isub) = ctw_decomps{icond,isub}.spec;
            ctwX(icond,isub) = ctw_decomps{icond,isub}.interact;
        end
    end
end

ctwPX = ctwP + ctwX; % add phactw + interaction

% Combine the data into a single vector
data = [ctwS; ctwPX];
cond = repmat([1:length(header)]',2,length(ctw_decomps));
comp = [zeros(length(header), length(ctw_decomps)); ones(length(header), length(ctw_decomps))];
ID   = repmat(1:12,8,1);
data = data(:);
cond = cond(:);
comp = comp(:);
ID = categorical(ID(:));

X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRctw = out.LRStat(out.Model=='mdl2');
Pctw = out.pValue(out.Model=='mdl2');



%% Extract output, PermEn

pe32P = nan(length(pe_decomps),length(header))';
pe32S = nan(length(pe_decomps),length(header))';
pe32X = nan(length(pe_decomps),length(header))';

pe64P = nan(length(pe_decomps),length(header))';
pe64S = nan(length(pe_decomps),length(header))';
pe64X = nan(length(pe_decomps),length(header))';

for icond = 1:length(header)
    for isub = 1:length(pe_decomps)
        if ~isempty(pe_decomps{icond,isub}) %&& ~strcmp(subs1{isub},bl)
            idx32= find(pe_decomps{icond,isub}.lags==32);
            idx64= find(pe_decomps{icond,isub}.lags==64);
            pe32P(icond,isub) = pe_decomps{icond,isub}.phase(idx32);
            pe32S(icond,isub) = pe_decomps{icond,isub}.spec(idx32);
            pe32X(icond,isub) = pe_decomps{icond,isub}.interact(idx32);
            pe64P(icond,isub) = pe_decomps{icond,isub}.phase(idx64);
            pe64S(icond,isub) = pe_decomps{icond,isub}.spec(idx64);
            pe64X(icond,isub) = pe_decomps{icond,isub}.interact(idx64);
        end
    end
end

pe32PX = pe32P + pe32X; % add phape + interaction
pe64PX = pe64P + pe64X; % add phape + interaction

% Combine the data into a single vector (PermEn32)
data = [pe32S; pe32PX];
cond = repmat([1:length(header)]',2,length(pe_decomps));
comp = [zeros(length(header), length(pe_decomps)); ones(length(header), length(pe_decomps))];
ID   = repmat(1:12,8,1);
data = data(:);
cond = cond(:);
comp = comp(:);
ID = categorical(ID(:));

X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRpe32 = out.LRStat(out.Model=='mdl2');
Ppe32 = out.pValue(out.Model=='mdl2');


% Now do again for PermEn64
data = [pe64S; pe64PX];
data = data(:);


X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRpe64 = out.LRStat(out.Model=='mdl2');
Ppe64 = out.pValue(out.Model=='mdl2');

%% Write as a table

%     Create a table (T) containing the measures, p-values, and log-likelihood ratios.
%     Write the table to a CSV file named 'fMEGEntropyDecomp.csv'.

T = table({'CTW','LZC','PermEn32','PermEn64','mMSE','mSampEn'}',...
    round([Pctw Plzc Ppe32 Ppe64 Pmse Pse],3,'significant')',...
    round([LRctw LRlzc LRpe32 LRpe64 LRmse LRse],3,'significant')',...
    'VariableNames',{'Measures','P-values','LR-stat'})

writetable(T,'fMEGEntropyDecomp.csv')



%% Generate figures

% This block of code generates violin plots to visually represent the 
% change in different entropy measures' amplitude and non-amplitude 
% components across young/old conditions. 

%%% CTW %%%
%     Use myfigure2 to create a new figure.
%     Generate a violin plot using myviolin to visualize the change in CTW components (amplitude and non-amplitude).
%     Annotate the figure with Cohen's d effect size.
%     Save the figure in SVG format.
%     This code repeats below for each entropy measure

myfigure2
myviolin(-1.*[nanmean(ctwS); nanmean(ctwPX)]','bw',0.01,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
ylabel('Change in CTW')
title(sprintf('Cohen''s d = %1.2f',cohens_d(nanmean(ctwS), nanmean(ctwPX))),'fontsize',20)
al = max(abs(ylim))+0.01;
ylim([-al al])
plot(linspace(0,3,10),zeros(1,10),'k--')
xlim([0.5 2.5])
makefigpretty
tx = gca().YTick;
ylim([min(tx) max(tx)])
yticks(tx)
print('-dsvg','./Figures/CTW_decomp')

%% LZC
myfigure2
myviolin(-1.*[nanmean(lzcS); nanmean(lzcPX)]','bw',1.5,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
ylabel('Change in LZC')
title(sprintf('Cohen''s d = %1.2f',cohens_d(nanmean(lzcS), nanmean(lzcPX))),'fontsize',20)
al = max(abs(ylim))+0.01;
ylim([-al al])
plot(linspace(0,3,10),zeros(1,10),'k--')
xlim([0.5 2.5])
makefigpretty
tx = gca().YTick;
ylim([min(tx) max(tx)])
yticks(tx)
print('-dsvg','./Figures/LZC_decomp')


%% PE32
myfigure2
myviolin(-1.*[nanmean(pe32S(:,:,idx32)); nanmean(pe32PX(:,:,idx32))]','bw',0.025,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
ylabel('Change in PE32')
title(sprintf('Cohen''s d = %1.2f',cohens_d(nanmean(pe32S), nanmean(pe32PX))),'fontsize',20)
al = max(abs(ylim))+0.01;
ylim([-al al])
yticks([-0.5 -0.25 0 0.25 0.5])
plot(linspace(0,3,10),zeros(1,10),'k--')
xlim([0.5 2.5])
makefigpretty
tx = gca().YTick;
ylim([min(tx) max(tx)])
yticks(tx)
print('-dsvg','./Figures/PE32_decomp')

%% PE64

myfigure2
myviolin(-1.*[nanmean(pe64S(:,:,idx32)); nanmean(pe64PX(:,:,idx32))]','bw',0.025,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
ylabel('Change in PE64')
title(sprintf('Cohen''s d = %1.2f',cohens_d(nanmean(pe64S), nanmean(pe64PX))),'fontsize',20)
al = max(abs(ylim))+0.01;
ylim([-al al])
yticks([-0.5 -0.25 0 0.25 0.5])
plot(linspace(0,3,10),zeros(1,10),'k--')
xlim([0.5 2.5])
makefigpretty
tx = gca().YTick;
ylim([min(tx) max(tx)])
yticks(tx)
print('-dsvg','./Figures/PE64_decomp')

%% mSampEn
myfigure2
myviolin(-1.*[nanmean(seS); nanmean(sePX)]','bw',0.005,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
ylabel('Change in mSampEn')
title(sprintf('Cohen''s d = %1.2f',cohens_d(nanmean(seS), nanmean(sePX))),'fontsize',20)
al = max(abs(ylim))+0.01;
ylim([-al al])
plot(linspace(0,3,10),zeros(1,10),'k--')
xlim([0.5 2.5])
makefigpretty
tx = gca().YTick;
ylim([min(tx) max(tx)])
yticks(tx)
print('-dsvg','./Figures/mSampEn_decomp')

%% mMSE

myfigure2
myviolin(-1.*[nanmean(mseS); nanmean(msePX)]','bw',0.033,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
ylabel('Change in mMSE')
title(sprintf('Cohen''s d = %1.2f',cohens_d(nanmean(mseS), nanmean(msePX))),'fontsize',20)
al = max(abs(ylim))+0.01;
ylim([-al al])
plot(linspace(0,3,10),zeros(1,10),'k--')
xlim([0.5 2.5])
makefigpretty
tx = gca().YTick;
ylim([min(tx) max(tx)])
yticks(tx)
print('-dsvg','./Figures/mMSE_decomp')


%% HELPER FUNCTION

% This helper function takes an input array of logical values (x) and
% converts each element to a corresponding string value. The output is a 
% cell array of strings (str), where each element represents the logical 
% value in string form.

% Input:
% 
%     x: Array of logical values.
% 
% Output:
% 
%     str: Cell array of strings representing the logical values in x.
% 
% Processing:
% 
%     The function initializes a cell array str with the same size as the input x.
%     It then iterates through each element of x.
%     For each element, if it is true, the corresponding string in str is set to 'TRUE'. If it is false, the string is set to 'FALSE'.
%     The final cell array str contains string representations of the logical values in x.

function[str] = logical2str(x)
str = cell(size(x));
for i = 1:numel(x)
    if x(i)
        str{i} = 'TRUE';
    else
        str{i} = 'FALSE';
    end
end
end


