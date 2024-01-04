% Joel Frohlich
% University of Tuebingen
% This is the script to get all entropy measures for fetal data. When
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

% These lines below fetch the ERF templates (again, the code uses this
% info, but you can ignore the corresponding output) 
load grand_average_ERFs ERFdev ERFstd ID % load ERFs from older fetuses 
% Take the average of ERF templates, averging first within subjects
ERFdev = nestedmean(ERFdev,ID);
ERFstd = nestedmean(ERFstd,ID);

%%% Now we get started ... 
%
%     Load Data:
%         Load data from the file "Data_traces.mat."
% 
%     Set Parameters:
%         nsur: Number of surrogates to generate (set to 100).
%         Set the random number generator seed using rng(111).

load Data_traces.mat
nsur = 100; % number of surrogates to generate 
rng(111) % random number generator seed

% Each trial-averaged ERF traces is -200 to 3000 ms references to the onset
% of the first tone in the sequence. Each tone is 200 ms in duration with
% 400 ms in between tones (thus, delta-T betwee tone *onsets* is 600 ms). The
% data are sampled at 610.3516 Hz. To take the data from the fourth tone
% onward *only*, set portion to 'end' and the signal will only be analyzed
% from 1800 ms - 3000 ms. 

fs = 610.3516; % sampling rate
portion = 'all'; % which part of the trial should we analyze? 
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

%     Label Blocks:
%         Identify A and B blocks based on the block labels in the Datatable.
%         A block corresponds to the global rule 'ssss,' and B block corresponds to 'sssd.'
%         Populate arrays Aidx and Bidx with indices corresponding to A and B blocks.

Aidx = [];
Bidx = [];

for irow = 1:size(Datatable,1)
    if strcmp(Datatable.Blocklabel(irow,:),'A (ssss)')
        Aidx = [Aidx irow];
    elseif strcmp(Datatable.Blocklabel(irow,:),'B (sssd)')
        Bidx = [Bidx irow];
    else
        error('condition not recognized')
    end
end

%% Stochaticity test 
% This block of code conducts a stochasticity test on the global deviants 
% and global standards for two conditions ('A' and 'B').

fprintf('Stochaticity test\n')

% Stochasticity Test 
% 
%     Extract data (dat) 
%     Extract corresponding subject IDs (sub) for condition.
%     Apply the stochasticity test (StoTest) to the data using the function StoTest(dat,fs,ERFdev).
%     NOTE: The function takes the ERF template as the third input argument
%     and does the computation both with and without template subtraction,
%     you can just IGNORE the output with the template subtraction! 
%     Save the results to a file with a formatted name.

% Global deviant, ssss/sssd
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_deviants_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
Sto_devA = StoTest(dat,fs,ERFdev);
save(sprintf('Sto_deviantA%s_IAAFT2_%s',endstr,date),'Sto_devA','sub')

% Global deviant, sssd/ssss
dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_deviants_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
Sto_devB = StoTest(dat,fs,ERFdev);
save(sprintf('Sto_deviantB%s_IAAFT2_%s',endstr,date),'Sto_devB','sub')

% Global standard, ssss/ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_standards_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
Sto_stdA = StoTest(dat,fs,ERFstd);
save(sprintf('Sto_standardA%s_IAAFT2_%s',endstr,date),'Sto_stdA','sub')

% Global standard, sssd/sssd
dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_standards_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
Sto_stdB = StoTest(dat,fs,ERFstd);
save(sprintf('Sto_standardB%s_IAAFT2_%s',endstr,date),'Sto_stdB','sub')

%% Lempel-Ziv 

% As in the block of code before, we are extracting data for each subject
% in each condition; the ERF template is needed as input, but the
% corresponding output can be ignored. 
        
fprintf('LZC\n')
% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_deviants_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
LZC_devA = LempelZiv(dat,nsur,fs,ERFdev);

save(sprintf('LZC_deviantA%s_IAAFT1_%s',endstr,date),'LZC_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_deviants_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
LZC_devB = LempelZiv(dat,nsur,fs,ERFdev);

save(sprintf('LZC_deviantB%s_IAAFT1_%s',endstr,date),'LZC_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_standards_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
LZC_stdA = LempelZiv(dat,nsur,fs,ERFstd);


save(sprintf('LZC_standardA%s_IAAFT1_%s',endstr,date),'LZC_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_standards_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
LZC_stdB = LempelZiv(dat,nsur,fs,ERFstd);

save(sprintf('LZC_standardB%s_IAAFT1_%s',endstr,date),'LZC_stdB','sub')
 

%% Permutation Entropy 32 ms 

% Now we do everything again for PermEn32, same caveat about template
% subtraction applies here. 

fprintf('PermEn32\n')
t = 32; % timescale (ms)

% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_deviants_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
PermEn32_devA = PermutationH(dat,nsur,fs,t,ERFdev);

save(sprintf('PermEn32_deviantA%s_IAAFT1_%s',endstr,date),'PermEn32_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_deviants_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
PermEn32_devB = PermutationH(dat,nsur,fs,t,ERFdev);

save(sprintf('PermEn32_deviantB%s_IAAFT1_%s',endstr,date),'PermEn32_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_standards_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
PermEn32_stdA = PermutationH(dat,nsur,fs,t,ERFstd);

save(sprintf('PermEn32_standardA%s_IAAFT1_%s',endstr,date),'PermEn32_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_standards_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
PermEn32_stdB = PermutationH(dat,nsur,fs,t,ERFstd);

save(sprintf('PermEn32_standardB%s_IAAFT1_%s',endstr,date),'PermEn32_stdB','sub')

%% Permutation Entropy 64 ms 

% You get the idea ... 

fprintf('PermEn64\n')
t = 64; % timescale (ms)

% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_deviants_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
PermEn64_devA = PermutationH(dat,nsur,fs,t,ERFdev);

save(sprintf('PermEn64_deviantA%s_IAAFT1_%s',endstr,date),'PermEn64_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_deviants_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
PermEn64_devB = PermutationH(dat,nsur,fs,t,ERFdev);

save(sprintf('PermEn64_deviantB%s_IAAFT1_%s',endstr,date),'PermEn64_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_standards_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
PermEn64_stdA = PermutationH(dat,nsur,fs,t,ERFstd);

save(sprintf('PermEn64_standardA%s_IAAFT1_%s',endstr,date),'PermEn64_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_standards_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
PermEn64_stdB = PermutationH(dat,nsur,fs,t,ERFstd);

save(sprintf('PermEn64_standardB%s_IAAFT1_%s',endstr,date),'PermEn64_stdB','sub')

%% CTW

fprintf('CTW\n')
% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_deviants_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
CTW_devA = CTW(dat,nsur,fs,ERFdev);

save(sprintf('CTW_deviantA%s_IAAFT1_%s',endstr,date),'CTW_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_deviants_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
CTW_devB = CTW(dat,nsur,fs,ERFdev);

save(sprintf('CTW_deviantB%s_IAAFT1_%s',endstr,date),'CTW_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_standards_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
CTW_stdA = CTW(dat,nsur,fs,ERFstd);

save(sprintf('CTW_standardA%s_IAAFT1_%s',endstr,date),'CTW_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_standards_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
CTW_stdB = CTW(dat,nsur,fs,ERFstd);

save(sprintf('CTW_standardB%s_IAAFT1_%s',endstr,date),'CTW_stdB','sub')

%% mMSE (modified multiscale entropy, also encompasses mSampEn as the first time scale) 

% Note - I sometimes abbreviate mMSE just as MSE in scripts 

fprintf('MSE\n')
% Global deviant, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_deviants_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
MSE_devA = MSEntropy(dat,nsur,fs,ERFdev);

save(sprintf('MSE_deviantA%s_IAAFT1_%s',endstr,date),'MSE_devA','sub')

% Global deviant, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_deviants_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
MSE_devB = MSEntropy(dat,nsur,fs,ERFdev);

save(sprintf('MSE_deviantB%s_IAAFT1_%s',endstr,date),'MSE_devB','sub')

% Global standard, ssss
dat =cell(1,1);
sub =nan(1,length(Aidx));
for i = 1:length(Aidx)
    dat{i} = Datatable.global_standards_normalized{Aidx(i)}(start:stop);
    sub(i) = Datatable.ID(Aidx(i));
end
MSE_stdA = MSEntropy(dat,nsur,fs,ERFstd);

save(sprintf('MSE_standardA%s_IAAFT1_%s',endstr,date),'MSE_stdA','sub')

% Global standard, sssd

dat =cell(1,1);
sub =nan(1,length(Bidx));
for i = 1:length(Bidx)
    dat{i} = Datatable.global_standards_normalized{Bidx(i)}(start:stop);
    sub(i) = Datatable.ID(Bidx(i));
end
MSE_stdB = MSEntropy(dat,nsur,fs,ERFstd);

save(sprintf('MSE_standardB%s_IAAFT1_%s',endstr,date),'MSE_stdB','sub')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% helper functions

%% Stochaticity 

% Input Parameters: 
%   alldata (cell array of data for multiple subjects)
%   fs (sampling rate)
%   meanERF (mean Evoked Response Field, it gets used here for the template subtraction, even though the template subtraction idea was ulimtately abandoned)

function[output] = StoTest(alldata,fs,meanERF)

% SUMMARY 
% 
% Initialization:
% 
%     Set a fixed seed for the random number generator, to be fed to the
%     surrogate data function. Initialize two cell arrays (result1 and
%     result2) to store the results of the stochasticity test for each
%     subject.
% 
% Loop Over Subjects:
% 
%     Iterate through each subject's data in the input alldata.
% 
% Subject-specific Processing:
% 
%     Print a message indicating the current subject being processed:
%     fprintf('Subject %i out of %i\n',isub,size(alldata,2)). Extract the
%     data for the current subject: dat = alldata{isub}.
% 
% Stochasticity Test 1:
% 
%     Apply the stochasticity test (stochastic_test_JF) to the raw data
%     (dat) using a specified seed (seed). Depending on the result
%     (tmpres), categorize the subject as either 'stochastic' or
%     'deterministic' in result1{isub}.
% 
% Stochasticity Test 2:
% 
%     Apply the same stochasticity test to the data after subtracting the
%     mean Evoked Response Field (meanERF). Categorize the subject based on
%     the result in result2{isub}.
% 
% Output:
% 
%     Store the results in a structure (output) with fields full (for the
%     test on the raw data) and subtr (for the test on data after
%     subtracting meanERF).
% 

seed = 123;

result1 = cell(1,size(alldata,2)); % subjects
result2 = cell(1,size(alldata,2)); % subjects

for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))

    dat = alldata{isub};


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

% Same inputs as before with stochasticity test, save caveat regarding data
% subtraction ... 

function[LZC] = LempelZiv(alldata,nsur,fs,meanERF)

% Initialization:
% 
%     Set a fixed seed for the random number generator: seed = 123.
%     Initialize variables for storing Lempel-Ziv complexity results:
%     LZC.full, LZC.subtr, LZC.trunc, and rsurr.
% 
% Loop Over Subjects:
% 
%     Iterate through each subject's data in the input alldata.
% 
% Subject-specific Processing:
% 
%     Print a message indicating the current subject being processed:
%     fprintf('Subject %i out of %i\n',isub,size(alldata,2)). Extract the
%     raw data and data after subtracting meanERF: dat = alldata{isub} and
%     dat2 = dat - meanERF.
% 
% Lempel-Ziv Complexity (Full Data):
% 
%     Calculate Lempel-Ziv complexity for the full data using the LZ76
%     function. Store the result in LZC.full(isub).
% 
% Lempel-Ziv Complexity (Subtracted MeanERF):
% 
%     Calculate Lempel-Ziv complexity for the data after subtracting
%     meanERF. Store the result in LZC.subtr(isub).
% 
% Generate Surrogates:
% 
%     Generate surrogates using the surrogate_JF function with the 'IAAFT1'
%     method. Calculate Lempel-Ziv complexity for the truncated signal
%     after surrogates generation. Store the result in LZC.trunc(isub).
% 
% Lempel-Ziv Complexity (Surrogates):
% 
%     Use parallel processing (parfor) to iterate over surrogates and
%     calculate Lempel-Ziv complexity. Store the results in the rsurr
%     matrix.
% 
% Output:
% 
%     Store the Lempel-Ziv complexity results and surrogate results in the
%     LZC structure.

seed = 123;
tic

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
        rng(isur,'twister') % we must include this line for rng to be set inside the parfor loop
        javaaddpath('./EntRate-master/private/vmm/vmm.jar')
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

% Similar to the above helper function, but for permutation entropy rather
% than Lempel-Ziv. 
%
% Noteable differences: 
% This function has additional parameters: t (time lag parameter in ms) and m (embedding dimension).

function[PE] = PermutationH(alldata,nsur,fs,t,meanERF)
seed = 123;
m = 3;
tau = round((t/1000)/(1/fs)); % convert from ms to samples

PE.full  = nan(1,size(alldata,2)); % subjects
PE.subtr = nan(1,size(alldata,2)); % subjects
PE.trunc = nan(1,size(alldata,2)); % subjects
rsurr    = nan(size(alldata,2),nsur); % subjects x surrogates


for isub = 1:size(alldata,2)
    fprintf('Subject %i out of %i\n',isub,size(alldata,2))
    dat = alldata{isub};

    PE.full(isub) = PermEn(dat,m,tau);
    PE.subtr(isub) = PermEn(dat-meanERF,m,tau);

    % scramble amplitudes (NOTE: as best I recall, this analysis was abandoned
    % and we did NOT use PE.randamp in the manuscript) 
    phi = angle(fft(dat));
    amp = 10.^randn(size(dat,1),size(dat,2)); % log-normal distribution of random amplitudes 
    dat_ra = ifft(amp.*exp(1j.*phi), 'symmetric'); % random amplitdues
    PE.randamp(isub) = PermEn(dat_ra,m,tau);

    % Genreate surrogates
    [dat_surr,params] = surrogate_JF(dat,nsur,'IAAFT1',true,fs,seed); 
    dat_surr = dat_surr';
    cutsignal = params.cutsig; % use truncated signal;
    PE.trunc(isub) = PermEn(cutsignal,m,tau);

    parfor isur = 1:nsur
        rng(isur,'twister') % we must include this line for rng to be set inside the parfor loop
        javaaddpath('./EntRate-master/private/vmm/vmm.jar')
        if mod(isur,nsur/5)==0
            fprintf('          Surrogate %i out %i\n',isur,nsur)
        end
        sr = dat_surr(:,isur);
        assert(length(sr)==length(cutsignal),'Surrogate has different length than truncated signal')
        % Surrogate
        rsurr(isub,isur) =  PermEn(sr',m,tau);
    end
end

% add output to data structure, including parametric information such as PE.m, PE.tau, PE.t, PE.fs (embedding dimension, delay, time window, and sampling rate, respectively).
PE.rsurr = rsurr;
PE.m = m;
PE.tau = tau;
PE.t = t;
PE.fs = fs;

end


%% CTW

% Does like above, but for context-tree weighted entropy. 
% This function has unique parameters: nsym (number of symbols to use), vmo (max model order), wsize (window size for time-resolved entropy), and overlap (overlap for time-resolved entropy).

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
        rng(isur,'twister') % we must include this line for rng to be set inside the parfor loop
        javaaddpath('./EntRate-master/private/vmm/vmm.jar')
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

% The helper function for mMSE/mSampEn

function[MSE] = MSEntropy(alldata,nsur,fs,meanERF)
seed = 123;
nscale = 20; % number of timescales to compute 
rng(909023920) % just to be safe, and also since we don't have rng being 
% called from each line of a parfor loop (parfor doesn't work here, sliced variable)


MSE.full   = nan(nscale,size(alldata,2)); % scales x subjects
MSE.subtr  = nan(nscale,size(alldata,2)); % scales x subjects
MSE.trunc  = nan(nscale,size(alldata,2)); % scales x subjects
rsurr      = nan(nscale,size(alldata,2),nsur); % scales x subjects x surrogates

% MSE Configuration:
% 
%     Configuration structure (cfg) for MSE:
%         cfg.fs: Sampling rate.
%         cfg.scales: Timescales for MSE calculation.
%         cfg.type: Type of MSE algorithm ('Xie').
%         cfg.r: state space radius parameter for MSE, in standard deviations of data
%         cfg.dynr: Dynamic R parameter (true = set new state space radius for each timescale).
%         cfg.gpu: Flag for GPU usage (true = use GPU for calculations).

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
    % Note: this loop doesn't work as a parfor loop because there is a sliced variable
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


