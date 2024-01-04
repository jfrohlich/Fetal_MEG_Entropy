% Joel Frohlich
% University of Tuebingen
% This is the script to get all entropy measures for neonatal data. When
% surrogate data are computed with any outside function, they are generated
% using the 'IAAFT1' method.

clearvars
dbstop if error

% NOTE: the template subtraction step was NOT done in the manuscript.
% Although this code will compute a version of each entropy measure with
% the template (i.e., the average ERF from older subjects) subtracted, this
% field of the output is unnecessary if you are just trying to replicate the
% paper. If you want to ignore it, you can just feed the functions zeros for the
% templates. 

warning ON VERBOSE
warning off MATLAB:colon:nonIntegerIndex

% These lines below fetch the ERF templates (again, the code uses this
% info, but you can ignore the corresponding output) 
load time_traces_data_manuscript
load Newborn_ERFs.mat % ERFs from newborns
% Take the average ERF across datasets to use as template
ERFdev = nestedmean(ERFdev,IDdev);
ERFstd = nestedmean(ERFstd,IDstd);

% Correct 0s written as letter o
data_tab.ID(find(data_tab.ID=="Co12-2-1A")) = "C012-2-1A";
nsur = 100; % number of surrogates to generate
rng(111) % seed added 10.03.22

% Each trial-averaged ERF traces is -200 to 3000 ms references to the onset
% of the first tone in the sequence. Each tone is 200 ms in duration with
% 400 ms in between tones (thus, delta-T betwee tone *onsets* is 600 ms). The
% data are sampled at 610.3516 Hz. To take the data from the fourth tone
% onward *only*, set portion to 'end' and the signal will only be analyzed
% from 1800 ms - 3000 ms. 

fs = 610.3516; % sampling rate (Hz)
portion = 'all'; % just analyze response to fourth tone? 
endstr = ''; % to append to end of file name

% Define the variable portion to determine which portion of the data to analyze.
switch portion
    case 'beginning' % just analyze before the 4th tone
        start = 1;
        stop = round(fs*2);
        endstr = strcat(endstr,'_T1-3only');
    case 'end'
        % if we want just the fourth tone onward, discard the first 2000 ms
        % i.e. -200 -- +1800 ms. 
        start = round(fs*2);
        stop = ceil(fs*3.2);
        endstr = strcat(endstr,'_T4only');
    case 'all' % use all the data from each trial (this is what we used)
        start = 1;
        stop = ceil(fs*3.2);
end

Aidx = [];
Bidx = [];
Cidx = [];
Didx = [];

for irow = 1:size(data_tab,1)
    checkme = data_tab.condition(irow,:);
    if ~strcmp(checkme,"time")
        if strcmp(checkme,"global_dev") % ssss (rule: sssd)neo
            Aidx = [Aidx irow];
        elseif strcmp(checkme,"local_global_dev") % sssd (rule: ssss)
            Bidx = [Bidx irow];
        elseif strcmp(checkme,"sandard") % ssss (rule: ssss)
            Cidx = [Cidx irow];
        elseif strcmp(checkme,"local_dev") % sssd (rule: sssd)
            Didx = [Didx irow];
        end
    end
end

keyboard

%%% Note - the code structure below closely mirrors the code in the corresponding
%%% fetal script and uses the same helper functions -- please refer to
%%% comments in that script if you need further clarification! 

%%
%%%%%%%%%%%%%%%% Lempel-Ziv %%%%%%%%%%%%%%%%%  

% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
oldsub = nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = data_tab.data(Aidx(i),start:end);
    try
        ID = sscanf(data_tab.ID{Aidx(i)},'C00%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    catch
        ID = sscanf(data_tab.ID{Aidx(i)},'C0%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    end
    oldsub(i) = sscanf(data_tab.ID{Aidx(i)},'C%i');
end
LZC_devA = LempelZiv(dat,nsur,fs,ERFdev);
save(sprintf('LZC_neonate_IAAFT1_deviantA%s_%s',endstr,date),'LZC_devA','sub','oldsub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
oldsub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = data_tab.data(Bidx(i),start:end);
    try
        ID = sscanf(data_tab.ID{Bidx(i)},'C00%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    catch
        ID = sscanf(data_tab.ID{Bidx(i)},'C0%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    end
    oldsub(i) = sscanf(data_tab.ID{Bidx(i)},'C%i');
end
LZC_devB = LempelZiv(dat,nsur,fs,ERFdev);
save(sprintf('LZC_neonate_IAAFT1_deviantB%s_%s',endstr,date),'LZC_devB','sub','oldsub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Cidx));
oldsub =nan(1,length(Cidx));
for i = 1:length(Cidx)
    dat{i} = data_tab.data(Cidx(i),start:end);
    try
        ID = sscanf(ddata_tab.ID{Cidx(i)},'C00%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    catch
        ID = sscanf(data_tab.ID{Cidx(i)},'C0%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    end
    oldsub(i) = sscanf(data_tab.ID{Cidx(i)},'C%i');
end
LZC_stdA = LempelZiv(dat,nsur,fs,ERFstd);
save(sprintf('LZC_neonate_IAAFT1_standardA%s_%s',endstr,date),'LZC_stdA','sub','oldsub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Didx));
oldsub =nan(1,length(Didx));
for i = 1:length(Didx)
    dat{i} = data_tab.data(Didx(i),start:end);
    try
        ID = sscanf(data_tab.ID{Didx(i)},'C00%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    catch
        ID = sscanf(data_tab.ID{Didx(i)},'C0%i');
        assert(~isempty(ID),'Failed to extract ID')
        sub(i) = ID;
    end
    oldsub(i) = sscanf(data_tab.ID{Didx(i)},'C%i');
end
LZC_stdB = LempelZiv(dat,nsur,fs,ERFstd);
save(sprintf('LZC_neonate_IAAFT1_standardB%s_%s',endstr,date),'LZC_stdB','sub','oldsub')


%%%%%%%%%%%%%%%% Permutation Entropy 32 ms %%%%%%%%%%%%%%%%%

fprintf('PermEn32\n')
t = 32; % timescale (ms)

% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = data_tab.data(Aidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Aidx(i)},'C%i');
end
PermEn32_devA = PermutationH(dat,nsur,fs,t,ERFdev);
save(sprintf('PermEn32_neonate_IAAFT1_deviantA%s_%s',endstr,date),'PermEn32_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = data_tab.data(Bidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Bidx(i)},'C%i');
end
PermEn32_devB = PermutationH(dat,nsur,fs,t,ERFdev);
save(sprintf('PermEn32_neonate_IAAFT1_deviantB%s_%s',endstr,date),'PermEn32_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Cidx));
for i = 1:length(Cidx)
    dat{i} = data_tab.data(Cidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Cidx(i)},'C%i');
end
PermEn32_stdA = PermutationH(dat,nsur,fs,t,ERFstd);
save(sprintf('PermEn32_neonate_IAAFT1_standardA%s_%s',endstr,date),'PermEn32_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Didx));
for i = 1:length(Didx)
    dat{i} = data_tab.data(Didx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Didx(i)},'C%i');
end
PermEn32_stdB = PermutationH(dat,nsur,fs,t,ERFstd);
save(sprintf('PermEn32_neonate_IAAFT1_standardB%s_%s',endstr,date),'PermEn32_stdB','sub')

%%%%%%%%%%%%%%%% Permutation Entropy 64 ms %%%%%%%%%%%%%%%%%

fprintf('PermEn64\n')
t = 64; % timescale (ms)

% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = data_tab.data(Aidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Aidx(i)},'C%i');
end
PermEn64_devA = PermutationH(dat,nsur,fs,t,ERFdev);
save(sprintf('PermEn64_neonate_IAAFT1_deviantA%s_%s',endstr,date),'PermEn64_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = data_tab.data(Bidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Bidx(i)},'C%i');
end
PermEn64_devB = PermutationH(dat,nsur,fs,t,ERFdev);
save(sprintf('PermEn64_neonate_IAAFT1_deviantB%s_%s',endstr,date),'PermEn64_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Cidx));
for i = 1:length(Cidx)
    dat{i} = data_tab.data(Cidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Cidx(i)},'C%i');
end
PermEn64_stdA = PermutationH(dat,nsur,fs,t,ERFstd);
save(sprintf('PermEn64_neonate_IAAFT1_standardA%s_%s',endstr,date),'PermEn64_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Didx));
for i = 1:length(Didx)
    dat{i} = data_tab.data(Didx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Didx(i)},'C%i');
end
PermEn64_stdB = PermutationH(dat,nsur,fs,t,ERFstd);
save(sprintf('PermEn64_neonate_IAAFT1_standardB%s_%s',endstr,date),'PermEn64_stdB','sub')


%%%%%%%%%%%%% CTW %%%%%%%%%%%%%%%

fprintf('CTW\n')
% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = data_tab.data(Aidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Aidx(i)},'C%i');
end
CTW_devA = CTW(dat,nsur,fs,ERFdev);
save(sprintf('CTW_neonate_IAAFT1_deviantA%s_%s',endstr,date),'CTW_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
sub = [];
for i = 1:length(Bidx)
    dat{i} = data_tab.data(Bidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Bidx(i)},'C%i');
end
CTW_devB = CTW(dat,nsur,fs,ERFdev);
save(sprintf('CTW_neonate_IAAFT1_deviantB%s_%s',endstr,date),'CTW_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Cidx));
for i = 1:length(Cidx)
    dat{i} = data_tab.data(Cidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Cidx(i)},'C%i');
end
CTW_stdA = CTW(dat,nsur,fs,ERFstd);
save(sprintf('CTW_neonate_IAAFT1_standardA%s_%s',endstr,date),'CTW_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Didx));
for i = 1:length(Didx)
    dat{i} = data_tab.data(Didx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Didx(i)},'C%i');
end
CTW_stdB = CTW(dat,nsur,fs,ERFstd);
save(sprintf('CTW_neonate_IAAFT1_standardB%s_%s',endstr,date),'CTW_stdB','sub')


%%%%%%%%%%%%% MSE %%%%%%%%%%%%%%%
fprintf('MSE\n')
% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = data_tab.data(Aidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Aidx(i)},'C%i');
end
MSE_devA = MSEntropy(dat,nsur,fs,ERFdev);
save(sprintf('MSE_neonate_IAAFT1_deviantA%s_%s',endstr,date),'MSE_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = data_tab.data(Bidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Bidx(i)},'C%i');
end
MSE_devB = MSEntropy(dat,nsur,fs,ERFdev);
save(sprintf('MSE_neonate_IAAFT1_deviantB%s_%s',endstr,date),'MSE_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Cidx));
for i = 1:length(Cidx)
    dat{i} = data_tab.data(Cidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Cidx(i)},'C%i');
end
MSE_stdA = MSEntropy(dat,nsur,fs,ERFstd);
save(sprintf('MSE_neonate_IAAFT1_standardA%s_%s',endstr,date),'MSE_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Didx));
for i = 1:length(Didx)
    dat{i} = data_tab.data(Didx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Didx(i)},'C%i');
end
MSE_stdB = MSEntropy(dat,nsur,fs,ERFstd);
save(sprintf('MSE_neonate_IAAFT1_standardB%s_%s',endstr,date),'MSE_stdB','sub')

%%%%%%%%%%%%% Stochasticity %%%%%%%%%%%%%%%
%%
fprintf('Stochasticity\n')
% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = data_tab.data(Aidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Aidx(i)},'C%i');
end

Sto_devA = StoTest(dat,fs,ERFdev);
save(sprintf('Sto_neonate_deviantA%s_IAAFT2_%s',endstr,date),'Sto_devA','sub')


% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = data_tab.data(Bidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Bidx(i)},'C%i');
end

Sto_devB = StoTest(dat,fs,ERFdev);
save(sprintf('Sto_neonate_deviantB%s_IAAFT2_%s',endstr,date),'Sto_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Cidx));
for i = 1:length(Cidx)
    dat{i} = data_tab.data(Cidx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Cidx(i)},'C%i');
end

Sto_stdA = StoTest(dat,fs,ERFstd);
save(sprintf('Sto_neonate_standardA%s_IAAFT2_%s',endstr,date),'Sto_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Didx));
for i = 1:length(Didx)
    dat{i} = data_tab.data(Didx(i),start:end);
    sub(i) =  sscanf(data_tab.ID{Didx(i)},'C%i');
end

Sto_stdB = StoTest(dat,fs,ERFstd);
save(sprintf('Sto_neonate_standardB%s_IAAFT2_%s',endstr,date),'Sto_stdB','sub')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions

%% Stochaticity 

function[output] = StoTest(alldata,fs,meanERF)
seed = 123;

result1 = cell(1,size(alldata,2)); % subjects
result2 = cell(1,size(alldata,2)); % subjects

for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))

    datr = alldata{isub};

    % lowpass filter at 10 Hz to match fetal data
    b = fir1(round(fs),10/fs*2,'low');
    dat = filtfilt(b,1,datr')'; 

   tmpres = stochastic_test_JF(dat,fs,seed);
    switch tmpres
        case 1
            result1{isub} = 'stochastic';
        case 0
            result1{isub} = 'deterministic';
    end
    clear tmpres % just to be safe
    tmpres = stochastic_test_JF(dat-meanERF,fs,seed);
    switch tmpres
        case 1
            result2{isub} = 'stochastic';
        case 0
            result2{isub} = 'deterministic';
    end

end

output.full = result1;
output.subtr = result2;

end


%% LZC


function[LZC] = LempelZiv(alldata,nsur,fs,meanERF)
tic
seed = 123;

LZC.full  = nan(1,size(alldata,2)); % subjects
LZC.subtr = nan(1,size(alldata,2)); % subjects
LZC.trunc = nan(1,size(alldata,2)); % subjects
rsurr     = nan(size(alldata,2),nsur); % subjects x surrogates


for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))

    dat = alldata{isub};
    dat2 = dat - meanERF;

    LZC.full(isub) = LZ76(dat>=median(dat));
    LZC.subtr(isub) = LZ76(dat2>=median(dat2));

    % Genreate surrogates
    [dat_surr,params] = surrogate_JF(dat,nsur,'IAAFT1',true,fs,seed); % preprocessing on
    dat_surr = dat_surr';
    cutsignal = params.cutsig; % Use truncated signal to compare with surrogates
    LZC.trunc(isub) = LZ76(cutsignal>=median(cutsignal));
    
    parfor isur = 1:nsur
        javaaddpath('/home/joel/Documents/MATLAB/MATLAB-20220201T161226Z-001/MATLAB/Universal_scripts/EntRate-master/private/vmm/vmm.jar')
        rng(isur,'twister') % we must include this line for rng to be set inside the parfor loop
        if mod(isur,nsur/5)==0
            fprintf('          Surrogate %i out %i\n',isur,nsur)
        end
        sr = dat_surr(:,isur);
        assert(length(sr)==length(cutsignal),'Surrogate has different length than truncated signal')
        % Surrogate
        rsurr(isub,isur) = LZ76(sr>=median(sr));
    end

end

LZC.rsurr = rsurr;

end

%% PermEn

function[PE] = PermutationH(alldata,nsur,fs,t,meanERF)
seed = 123;
m = 3;
tau = round((t/1000)/(1/fs));

PE.full  = nan(1,size(alldata,2)); % subjects
PE.subtr = nan(1,size(alldata,2)); % subjects
PE.trunc = nan(1,size(alldata,2)); % subjects
rsurr    = nan(size(alldata,2),nsur); % subjects x surrogates


for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))
    dat = alldata{isub};

    PE.full(isub) = PermEn(dat,m,tau);
    PE.subtr(isub) = PermEn(dat-meanERF,m,tau);

    % Genreate surrogates
    [dat_surr,params] = surrogate_JF(dat,nsur,'IAAFT1',true,fs,seed); %
    dat_surr = dat_surr';
    cutsignal = params.cutsig; % use truncated signal;
    PE.trunc(isub) = PermEn(cutsignal,m,tau);

    parfor isur = 1:nsur
        javaaddpath('/home/joel/Documents/MATLAB/MATLAB-20220201T161226Z-001/MATLAB/Universal_scripts/EntRate-master/private/vmm/vmm.jar')
        rng(isur,'twister') % we must include this line for rng to be set inside the parfor loop
        if mod(isur,nsur/5)==0
            fprintf('          Surrogate %i out %i\n',isur,nsur)
        end
        sr = dat_surr(:,isur);
        assert(length(sr)==length(cutsignal),'Surrogate has different length than truncated signal')
        % Surrogate
        rsurr(isub,isur) =  PermEn(sr',m,tau);
    end
end

PE.rsurr = rsurr;
PE.m = m;
PE.tau = tau;
PE.t = t;
PE.fs = fs;

end


%% CTW

function[CTW] = CTW(alldata,nsur,fs,meanERF)
seed = 123;
nsym = 2; % how many symbols to use
vmo = 100; % Max model order

CTW.full= nan(1,size(alldata,2)); % subjects
CTW.trunc = nan(1,size(alldata,2)); % subjects
rsurr = nan(size(alldata,2),nsur); % subjects x surrogates
N = length(alldata{1});

wsize = 100; % window for time resolved entropy
overlap = 0.9; % overlap for time resolved entropy

for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))
    dat = alldata{isub};
    tmp = timeresentropy_preissl(dat,wsize,overlap);
    if ~isfield(CTW,'trh')
        CTW.trh = nan(size(alldata,2),length(tmp));
    end
    if ~isfield(CTW,'erf')
        CTW.erf = nan(size(alldata,2),N);
    end
    CTW.trh(isub,:) = tmp; % store time-resolved entropy
    CTW.erf(isub,:) = dat; % store ERF signal
    
    CTW.full(isub)  = CTWEntropyRate(dat,nsym,vmo);
    CTW.subtr(isub) = CTWEntropyRate(dat-meanERF,nsym,vmo);

    % Genreate surrogates
    [dat_surr,params] = surrogate_JF(dat,nsur,'IAAFT1',true,fs,seed); %
    dat_surr = dat_surr';
    cutsignal = params.cutsig; % use truncated signal
    CTW.trunc(isub) = CTWEntropyRate(cutsignal,nsym,vmo);
    
    parfor isur = 1:nsur
        javaaddpath('/home/joel/Documents/MATLAB/MATLAB-20220201T161226Z-001/MATLAB/Universal_scripts/EntRate-master/private/vmm/vmm.jar')
        rng(isur,'twister') % we must include this line for rng to be set inside the parfor loop
        if mod(isur,nsur/5)==0
            fprintf('          Surrogate %i out %i\n',isur,nsur)
        end
        sr = dat_surr(:,isur)';
        assert(length(sr)==length(cutsignal),'Surrogate has different length than truncated signal')
        % Surrogate
        rsurr(isub,isur) = CTWEntropyRate(dat_surr(:,isur)',nsym,vmo);
    end
        
end
CTW.rsurr = rsurr;

end


%% MSE

% This gets us both mSampEn (the first timescale of mMSE) and mMSE

function[MSE] = MSEntropy(alldata,nsur,fs,meanERF)
seed = 123;
nscale = 20;

MSE.full   = nan(nscale,size(alldata,2)); % scales x subjects
MSE.subtr  = nan(nscale,size(alldata,2)); % scales x subjects
MSE.trunc  = nan(nscale,size(alldata,2)); % scales x subjects
rsurr      = nan(nscale,size(alldata,2),nsur); % scales x subjects x surrogates

% Configuration for MSE
cfg = [];
cfg.fs = fs;
cfg.scales = 1:nscale;
cfg.type = 'Xie';
cfg.r = 0.15;
cfg.dynr = true;
cfg.gpu = true;


for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))
    dat = alldata{isub};
    MSE.full(:,isub) = ro_mse(zscore(dat),cfg);
    MSE.subtr(:,isub) = ro_mse(zscore(dat-meanERF),cfg);

    [dat_surr,params] = surrogate_JF(dat,nsur,'IAAFT1',true,fs,seed);
    dat_surr = dat_surr';
    cutsignal = params.cutsig; % IMPORTANT: Use the truncated signal for comparison with surrogates
    MSE.trunc(:,isub) = ro_mse(zscore(cutsignal),cfg);
    for isur = 1:nsur
    % parfor doesn't work here because of sliced variable - we have to use
    % a nomral for loop
        if mod(isur,nsur/5)==0
            fprintf('          Surrogate %i out %i\n',isur,nsur)
        end
        sr = dat_surr(:,isur)';
        assert(length(sr)==length(cutsignal),'Surrogate has different length than truncated signal')
        % Surrogate
        rsurr(:,isub,isur) = ro_mse(zscore(sr),cfg);
    end
end
MSE.rsurr = rsurr;

end



