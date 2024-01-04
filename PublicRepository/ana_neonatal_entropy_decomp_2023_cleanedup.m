% Joel Frohlich
% University of Tuebingen
%
% This is the entropy decomp which shows how amplitude and non-amplitude
% components of the signal change with age.
% 
% Similar script as ana_fetal_entropy_decomp_2023.m, but applied to
% neonatal data -- refer to the other script as well, which may contain
% applicable comments

% Summary 

%     Loading Data:
%         Load MEG data from a file 'time_traces_data_manuscript.mat'.
%         Remove the last row of the data table, which doesn't contain
%         subject data. 
%         Set the random number generator seed.
%         Define the number of surrogates (nsur) and the sampling rate (fs).
% 
%     Age Information Addition:
%         Read a table ('MeasurementQualityTable.csv') containing age information.
%         Add age in days to the data table based on subject IDs.
%         Correct a typo in subject IDs.
% 
%     Group Division:
%         Divide the data into young and old groups based on the median age.
%         Select specific trial conditions from the data.
% 
%     Entropy Decomposition:
%         Perform entropy decomposition for all entropy measures across
%         young/old conditions.
%         The decomposition is done for each subject and age group.
% 
%     Statistical Analysis:
%         Fit linear mixed-effects models to the entropy data for different conditions.
%         Perform log-likelihood ratio tests to compare models with and without a component related to age group.
%         Extract and save the results for statistical analysis.
% 
%     Data Extraction and Analysis:
%         Extract relevant values from the entropy decompositions.
%         Combine the data into a table for each entropy measure (mSampEn, mMSE, LZC, CTW, PermEn).
%         Fit linear mixed-effects models and perform log-likelihood ratio tests for each entropy measure.
% 
%     Output Saving:
%         Save the results of entropy decomposition and statistical analysis into separate files.
% 
%     Loading Data (Optional):
%         Check if the output data is already loaded; if not, load it from saved files.
% 
%     Results Extraction and Analysis (Optional):
%         Extract relevant values from entropy decompositions for each measure.
%         Combine data and perform statistical analysis to compare models with and without a component related to age group.

%% NEONATAL ENTROPY DECOMP


clearvars
load time_traces_data_manuscript
data_tab(end,:) = []; % remove last row, doesn't contain a subject
rng(2345894)
nsur = 250; % number of surrogates
fs = 610.3516; % sampling rate

warning('off', 'all')

% Add age in days

T = readtable('MeasurementQualityTable.csv');
data_tab.age = nan(size(data_tab,1),1);
IDs = unique(data_tab.ID);
data_tab.SID = nan(size(data_tab.ID));

for iid = 1:length(IDs)
    idstr = sprintf('%s',IDs(iid,:));
    idstr = idstr(1:4) % somehow this is necessary for contains to work
    haso = idstr=='o';
    if any(haso)
        idstr(haso) = '0'; % switch o to 0 (correct typo)
    end
    findme = find(contains(T.filename,idstr));
    AGE = unique(T.age_days_(findme));
    assert(length(AGE)==1,'More than one age found')
    data_tab.age(data_tab.ID == IDs(iid)) = AGE;
    data_tab.SID(data_tab.ID == IDs(iid)) = str2num(replace(idstr,'C',''));
end

data_tab.ID = [];

% Take only sssS trials
as = median(data_tab.age); % median split
Ty = data_tab(data_tab.age <= as,:);
To = data_tab(data_tab.age > as,:);

fprintf('Dividing infant data into young (%i days or younger) and old groups\n',as)

assert(length(unique(Ty.SID))==length(unique(To.SID)),'Young and old infants are unbalanced');
N = length(unique(Ty.SID));

pe_decomps = cell(4,N,N);
se_decomps = cell(4,N,N);
mse_decomps = cell(4,N,N);
lzc_decomps = cell(4,N,N);
ctw_decomps = cell(4,N,N);
header = {'ssss','sssS','sssd','sssD'};
young = unique(Ty.SID);
old = unique(To.SID);

tic

%%% Each entropy measure below is decomposed into phase, amplitude, and
%%% interaction terms ... 
%     Iterative Loop:
%         Iterate through each subject (use) and each condition (header).
%         Different cases based on the condition, e.g., 'ssss', 'sssS', 'sssd', 'sssD'.
%         Find corresponding rows in Ty and To for each condition, subject, and age.
%         If data is available for the condition, perform the decomposition
%         using the appropriate function. 
%         Store the output in the mse_decomps cell array.
% 
%     Print Progress:
%         Display a progress message for each subject.
%         Output the elapsed time.
% 
%     Save Results:
%         Save the decomposition results and headers in a MAT file.


%% modified MSE
fprintf('Starting entropy decomp for mMSE\n')
for isub = 1:length(young)
    for jsub = 1:length(old)
        for icond = 1:length(header)
            switch header{icond}
                case 'ssss'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='standard',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='standard',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mMSEDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    mse_decomps{strcmp(header,'ssss'),isub,jsub} = out;
                case 'sssS'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mMSEDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    mse_decomps{strcmp(header,'sssS'),isub,jsub} = out;
                case 'sssd'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mMSEDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    mse_decomps{strcmp(header,'sssd'),isub,jsub} = out;

                case 'sssD'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mMSEDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    mse_decomps{strcmp(header,'sssD'),isub,jsub} = out;
                otherwise
                    error('can''t read the condition')
            end
        end
    end

    fprintf('      Finished subject #%i (ID = %i), %1.1%% complete\n',isub,young(isub),isub/N*100)
    toc
end
save(sprintf('MSEDecompneonate2023_Nsur=%i',nsur),'mse_decomps','header')


%%

% modified Sample Entropy
fprintf('Starting entropy decomp for mSampEn\n')
for isub = 1:length(young)
    for jsub = 1:length(old)
        for icond = 1:length(header)
            switch header{icond}
                case 'ssss'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='standard',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='standard',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mSampEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    se_decomps{strcmp(header,'ssss'),isub,jsub} = out;
                case 'sssS'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mSampEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    se_decomps{strcmp(header,'sssS'),isub,jsub} = out;
                case 'sssd'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mSampEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    se_decomps{strcmp(header,'sssd'),isub,jsub} = out;

                case 'sssD'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = mSampEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    se_decomps{strcmp(header,'sssD'),isub,jsub} = out;
                otherwise
                    error('can''t read the condition')
            end
        end
    end

   fprintf('      Finished subject #%i (ID = %i), %1.1%% complete\n',isub,young(isub),isub/N*100)
    toc
end
save(sprintf('mSampEnDecompneonate2023_Nsur=%i',nsur),'se_decomps','header')



%% Permutation Entropy
fprintf('Starting entropy decomp for mPermEn\n')
for isub = 1:length(young)
    for jsub = 1:length(old)
        for icond = 1:length(header)
            switch header{icond}
                case 'ssss'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='standard',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='standard',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = PermEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    pe_decomps{strcmp(header,'ssss'),isub,jsub} = out;
                case 'sssS'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = PermEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    pe_decomps{strcmp(header,'sssS'),isub,jsub} = out;
                case 'sssd'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = PermEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    pe_decomps{strcmp(header,'sssd'),isub,jsub} = out;

                case 'sssD'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = PermEnDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    pe_decomps{strcmp(header,'sssD'),isub,jsub} = out;
                otherwise
                    error('can''t read the condition')
            end
        end
    end
    fprintf('      Finished subject #%i (ID = %i), %1.1%% complete\n',isub,young(isub),isub/N*100)
    toc
end
save(sprintf('PermEnDecompneonate2023_Nsur=%i',nsur),'pe_decomps','header')


%% Lempel-Ziv
fprintf('Starting entropy decomp for LZC\n')
for isub = 1:length(young)
    for jsub = 1:length(old)
        for icond = 1:length(header)
            switch header{icond}
                case 'ssss'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='standard',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='standard',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = LZDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    lzc_decomps{strcmp(header,'ssss'),isub,jsub} = out;
                case 'sssS'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = LZDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    lzc_decomps{strcmp(header,'sssS'),isub,jsub} = out;
                case 'sssd'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = LZDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    lzc_decomps{strcmp(header,'sssd'),isub,jsub} = out;

                case 'sssD'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = LZDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    lzc_decomps{strcmp(header,'sssD'),isub,jsub} = out;
                otherwise
                    error('can''t read the condition')
            end
            fprintf('.')
        end
    end
    fprintf('      Finished subject #%i (ID = %i), %1.1%% complete\n',isub,young(isub),isub/N*100)
    toc
end
save(sprintf('LZCDecompneonate2023_Nsur=%i',nsur),'lzc_decomps','header')

%% CTW
fprintf('Starting entropy decomp for CTW\n')
for isub = 1:length(young)
    for jsub = 1:length(old)
        for icond = 1:length(header)
            switch header{icond}
                case 'ssss'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='standard',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='standard',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = CTWDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    ctw_decomps{strcmp(header,'ssss'),isub,jsub} = out;
                case 'sssS'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = CTWDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    ctw_decomps{strcmp(header,'sssS'),isub,jsub} = out;
                case 'sssd'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition=='local_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition=='local_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = CTWDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    ctw_decomps{strcmp(header,'sssd'),isub,jsub} = out;

                case 'sssD'
                    % ssss block label for that subject and age
                    rowy = find(Ty.SID == young(isub)  & all(Ty.condition =='global_dev',2));
                    rowo = find(To.SID == old(jsub) & all(To.condition =='global_dev',2));
                    if isempty(rowy) || isempty(rowo), continue, end % continue on if no data for this Block
                    % young (high entropy) minus old (low entropy)
                    out = CTWDecomposition_fMEG({Ty.data(rowy,:)}, young(isub), ...
                        {To.data(rowo,:)}, old(jsub), nsur, fs,'newborns');

                    ctw_decomps{strcmp(header,'sssD'),isub,jsub} = out;
                otherwise
                    error('can''t read the condition')
            end
        end
    end

    fprintf('      Finished subject #%i (ID = %i), %1.1%% complete\n',isub,young(isub),isub/N*100)
    toc
end
save(sprintf('CTWDecompneonate2023_Nsur=%i',nsur),'ctw_decomps','header')


%% Starting from this line? load the data here (if you already ran the decomp earlier)

% If the decomposition results are not already present, load them from corresponding MAT files.


if ~exist('ctw_decomps','var')
    load('CTWDecompneonate2023_Nsur=250.mat')
    load('LZCDecompneonate2023_Nsur=250.mat')
    load('MSEDecompneonate2023_Nsur=250.mat')
    load('mSampEnDecompneonate2023_Nsur=250.mat')
    load('PermEnDecompneonate2023_Nsur=250.mat')
end

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

        phase_values = nan(1, length(se_decomps));
        spec_values = nan(1, length(se_decomps));
        interact_values = nan(1, length(se_decomps));
        for jsub = 1:length(se_decomps)
            if ~isempty(se_decomps{icond, isub,jsub})
                phase_values(jsub) = se_decomps{icond, isub, jsub}.phase;
                spec_values(jsub) = se_decomps{icond, isub, jsub}.spec;
                interact_values(jsub) = se_decomps{icond, isub, jsub}.interact;
            end
        end
        seP(icond, isub) = nanmean(phase_values);
        seS(icond, isub) = nanmean(spec_values);
        seX(icond, isub) = nanmean(interact_values);
    end
end



sePX = seP + seX; % add phase + interaction

% Combine the data into a single vector
data = [seS; sePX];
cond = repmat([1:length(header)]',2,length(se_decomps));
comp = [zeros(length(header), length(se_decomps)); ones(length(header), length(se_decomps))];
ID   = repmat(1:length(se_decomps),length(header)*2,1);
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

        phase_values = nan(1, length(mse_decomps));
        spec_values = nan(1, length(mse_decomps));
        interact_values = nan(1, length(mse_decomps));
        for jsub = 1:length(mse_decomps)
            if ~isempty(mse_decomps{icond, isub,jsub})
                phase_values(jsub) = mse_decomps{icond, isub, jsub}.phase;
                spec_values(jsub) = mse_decomps{icond, isub, jsub}.spec;
                interact_values(jsub) = mse_decomps{icond, isub, jsub}.interact;
            end
        end
        mseP(icond, isub) = nanmean(phase_values);
        mseS(icond, isub) = nanmean(spec_values);
        mseX(icond, isub) = nanmean(interact_values);
    end
end



msePX = mseP + mseX; % add phase + interaction

% Combine the data into a single vector
data = [mseS; msePX];
cond = repmat([1:length(header)]',2,length(mse_decomps));
comp = [zeros(length(header), length(mse_decomps)); ones(length(header), length(mse_decomps))];
ID   = repmat(1:length(mse_decomps),length(header)*2,1);
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

        phase_values = nan(1, length(lzc_decomps));
        spec_values = nan(1, length(lzc_decomps));
        interact_values = nan(1, length(lzc_decomps));
        for jsub = 1:length(lzc_decomps)
            if ~isempty(lzc_decomps{icond, isub,jsub})
                phase_values(jsub) = lzc_decomps{icond, isub, jsub}.phase;
                spec_values(jsub) = lzc_decomps{icond, isub, jsub}.spec;
                interact_values(jsub) = lzc_decomps{icond, isub, jsub}.interact;
            end
        end
        lzcP(icond, isub) = nanmean(phase_values);
        lzcS(icond, isub) = nanmean(spec_values);
        lzcX(icond, isub) = nanmean(interact_values);
    end
end



lzcPX = lzcP + lzcX; % add phase + interaction

% Combine the data into a single vector
data = [lzcS; lzcPX];
cond = repmat([1:length(header)]',2,length(lzc_decomps));
comp = [zeros(length(header), length(lzc_decomps)); ones(length(header), length(lzc_decomps))];
ID   = repmat(1:length(lzc_decomps),length(header)*2,1);
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

        phase_values = nan(1, length(ctw_decomps));
        spec_values = nan(1, length(ctw_decomps));
        interact_values = nan(1, length(ctw_decomps));
        for jsub = 1:length(ctw_decomps)
            if ~isempty(ctw_decomps{icond, isub,jsub})
                phase_values(jsub) = ctw_decomps{icond, isub, jsub}.phase;
                spec_values(jsub) = ctw_decomps{icond, isub, jsub}.spec;
                interact_values(jsub) = ctw_decomps{icond, isub, jsub}.interact;
            end
        end
        ctwP(icond, isub) = nanmean(phase_values);
        ctwS(icond, isub) = nanmean(spec_values);
        ctwX(icond, isub) = nanmean(interact_values);
    end
end



ctwPX = ctwP + ctwX; % add phase + interaction

% Combine the data into a single vector
data = [ctwS; ctwPX];
cond = repmat([1:length(header)]',2,length(ctw_decomps));
comp = [zeros(length(header), length(ctw_decomps)); ones(length(header), length(ctw_decomps))];
ID   = repmat(1:length(ctw_decomps),length(header)*2,1);
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

        phase_values32 = nan(1, length(pe_decomps));
        spec_values32 = nan(1, length(pe_decomps));
        interact_values32 = nan(1, length(pe_decomps));
        phase_values64 = nan(1, length(pe_decomps));
        spec_values64 = nan(1, length(pe_decomps));
        interact_values64 = nan(1, length(pe_decomps));

        for jsub = 1:length(pe_decomps)
            if ~isempty(pe_decomps{icond, isub,jsub})
                idx32= find(pe_decomps{icond,isub,jsub}.lags==32);
                idx64= find(pe_decomps{icond,isub,jsub}.lags==64);
                phase_values32(jsub) = pe_decomps{icond, isub, jsub}.phase(idx32);
                spec_values32(jsub) = pe_decomps{icond, isub, jsub}.spec(idx32);
                interact_values32(jsub) = pe_decomps{icond, isub, jsub}.interact(idx32);
                phase_values64(jsub) = pe_decomps{icond, isub, jsub}.phase(idx64);
                spec_values64(jsub) = pe_decomps{icond, isub, jsub}.spec(idx64);
                interact_values64(jsub) = pe_decomps{icond, isub, jsub}.interact(idx64);
            end
        end
        pe32P(icond, isub) = nanmean(phase_values32);
        pe32S(icond, isub) = nanmean(spec_values32);
        pe32X(icond, isub) = nanmean(interact_values32);
        pe64P(icond, isub) = nanmean(phase_values64);
        pe64S(icond, isub) = nanmean(spec_values64);
        pe64X(icond, isub) = nanmean(interact_values64);
    end
end



pe32PX = pe32P + pe32X; % add phase + interaction

% Combine the data into a single vector
data = [pe32S; pe32PX];
cond = repmat([1:length(header)]',2,length(pe_decomps));
comp = [zeros(length(header), length(pe_decomps)); ones(length(header), length(pe_decomps))];
ID   = repmat(1:length(pe_decomps),length(header)*2,1);
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


pe64PX = pe64P + pe64X; % add phase + interaction

% Combine the data into a single vector
data = [pe64S; pe64PX];
cond = repmat([1:length(header)]',2,length(pe_decomps));
comp = [zeros(length(header), length(pe_decomps)); ones(length(header), length(pe_decomps))];
ID   = repmat(1:length(pe_decomps),length(header)*2,1);
data = data(:);
cond = cond(:);
comp = comp(:);
ID = categorical(ID(:));

X = table(data,cond,comp,ID);
mdl1 = fitlme(X,'data ~ cond + (1|ID)');
mdl2 = fitlme(X,'data ~ cond + comp + (1|ID)');
out = compare(mdl1,mdl2); % log-likelihood ratio test
LRpe64 = out.LRStat(out.Model=='mdl2');
Ppe64 = out.pValue(out.Model=='mdl2');



%% Write as a table

T = table({'CTW','LZC','PermEn32','PermEn64','mMSE','mSampEn'}',...
    round([Pctw Plzc Ppe32 Ppe64 Pmse Pse],3,'significant')',...
    round([LRctw LRlzc LRpe32 LRpe64 LRmse LRse],3,'significant')',...
    'VariableNames',{'Measures','P-values','LR-stat'})

writetable(T,'neonate_decomp.csv')

%% Generate figure

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
%%xtickangle(45)
%xlabel('Component')
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
print('-dsvg','./Figures/CTW_decomp_neonate')

%% LZC
myfigure2
myviolin(-1.*[nanmean(lzcS); nanmean(lzcPX)]','bw',1.5,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
%xtickangle(45)
%xlabel('Component')
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
print('-dsvg','./Figures/LZC_decomp_neonate')


%% PE32
myfigure2
myviolin(-1.*[nanmean(pe32S(:,:,idx32)); nanmean(pe32PX(:,:,idx32))]','bw',0.025,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
%xtickangle(45)
%xlabel('Component')
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
print('-dsvg','./Figures/PE32_decomp_neonate')

%% PE64

myfigure2
myviolin(-1.*[nanmean(pe64S(:,:,idx32)); nanmean(pe64PX(:,:,idx32))]','bw',0.025,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
%xtickangle(45)
%xlabel('Component')
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
print('-dsvg','./Figures/PE64_decomp_neonate')

%% mSampEn
myfigure2
myviolin(-1.*[nanmean(seS); nanmean(sePX)]','bw',0.005,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
%xtickangle(45)
%xlabel('Component')
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
print('-dsvg','./Figures/mSampEn_decomp_neonate')

%% mMSE

myfigure2
myviolin(-1.*[nanmean(mseS); nanmean(msePX)]','bw',0.033,'facecolor',[0.5 1 1; 1 0.5 1],'medc',[]) %,'paired',true)
xticks([1 2])
xticklabels({'Amplitude','Non-amplitude'})
%xtickangle(45)
%xlabel('Component')
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
print('-dsvg','./Figures/mMSE_decomp_neonate')

