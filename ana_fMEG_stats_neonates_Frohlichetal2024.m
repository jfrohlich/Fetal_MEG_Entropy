% This is the script that plots the figures and does the stats
% Joel Frohlich
% University of Tuebingen

% General note - this script closely mirrors ana_fMEG_stats_fetal_2023.m,
% the comments of which may also be useful for understanding this script. 

clearvars
close all

%%% How to interpret the stimulus/rule labels:
% stim/rule
% DevA = ssss/sssd (global deviant)
% DevB = sssd/ssss (local and global deviant)
% StdA = ssss/ssss (standard)
% StdB = sssd/sssd (local deviant)

% The script begins by loading in the preprocessed MEG data from different
% stimulus/rule conditions that were presented in an auditory oddball paradigm.

T = readtable('HypothesisTesting.csv');
% Only uncomment below to correct errors in reading subject IDs
%D = readtable('NeonateSubjectDecoder.csv');

portion = 'all';
fs = 610.3516; % sampling rate
k = 3; % "spatial" (i.e., time x freq) interpolation factor for time-freq representations
warning off MATLAB:colon:nonIntegerIndex

% portion denotes which part of signal should be analysed. To replicate
% results from the manuscript, use 'all' (-200 - 3000 ms)
switch portion
    case 'beginning'
        tstr = '_T1-3only';
    case 'end'
        tstr = '_T4only';
    case 'all'
        tstr = '';
end

%%% LOAD THE DATA %%%
% Note: datestr must be set to the date of files you want to load

% LZC
datestr = '_23-Jun-2022';
load(sprintf('LZC_neonate_IAAFT1_deviantA%s%s',tstr,datestr))
mslzc_devA = median(LZC_devA.rsurr,2); % median across all surrogates
sub1 = sub;
oldsub1 = oldsub;

load(sprintf('LZC_neonate_IAAFT1_deviantB%s%s',tstr,datestr))
mslzc_devB = median(LZC_devB.rsurr,2); % median across all surrogates
sub2 = sub;
oldsub2 = oldsub;

load(sprintf('LZC_neonate_IAAFT1_standardA%s%s',tstr,datestr))
mslzc_stdA = median(LZC_stdA.rsurr,2); % median across all surrogates
sub3 = sub;
oldsub3 = oldsub;

load(sprintf('LZC_neonate_IAAFT1_standardB%s%s',tstr,datestr))
mslzc_stdB = median(LZC_stdB.rsurr,2); % median across all surrogates
sub4 = sub;
oldsub4 = oldsub;


load(sprintf('CTW_neonate_IAAFT1_deviantA%s%s',tstr,datestr))
msctw_devA = median(CTW_devA.rsurr,2); % median across all surrogates
assert(all(sub==oldsub1),'Subjects loaded in different order')

load(sprintf('CTW_neonate_IAAFT1_deviantB%s%s',tstr,datestr))
msctw_devB = median(CTW_devB.rsurr,2); % median across all surrogates
assert(all(sub==oldsub2),'Subjects loaded in different order')

load(sprintf('CTW_neonate_IAAFT1_standardA%s%s',tstr,datestr))
msctw_stdA = median(CTW_stdA.rsurr,2); % median across all surrogates
assert(all(sub==oldsub3),'Subjects loaded in different order')

load(sprintf('CTW_neonate_IAAFT1_standardB%s%s',tstr,datestr))
msctw_stdB = median(CTW_stdB.rsurr,2); % median across all surrogates
assert(all(sub==oldsub4),'Subjects loaded in different order')


% PermEn32

load(sprintf('PermEn32_neonate_IAAFT1_deviantA%s%s',tstr,datestr))
mspe32_devA = median(PermEn32_devA.rsurr,2); % median across all surrogates
assert(all(sub==oldsub1),'Subjects loaded in different order')

load(sprintf('PermEn32_neonate_IAAFT1_deviantB%s%s',tstr,datestr))
mspe32_devB = median(PermEn32_devB.rsurr,2); % median across all surrogates
assert(all(sub==oldsub2),'Subjects loaded in different order')

load(sprintf('PermEn32_neonate_IAAFT1_standardA%s%s',tstr,datestr))
mspe32_stdA = median(PermEn32_stdA.rsurr,2); % median across all surrogates
assert(all(sub==oldsub3),'Subjects loaded in different order')

load(sprintf('PermEn32_neonate_IAAFT1_standardB%s%s',tstr,datestr))
mspe32_stdB = median(PermEn32_stdB.rsurr,2); % median across all surrogates
assert(all(sub==oldsub4),'Subjects loaded in different order')

% PermEn64

load(sprintf('PermEn64_neonate_IAAFT1_deviantA%s%s',tstr,datestr))
mspe64_devA = median(PermEn64_devA.rsurr,2); % median across all surrogates
assert(all(sub==oldsub1),'Subjects loaded in different order')

load(sprintf('PermEn64_neonate_IAAFT1_deviantB%s%s',tstr,datestr))
mspe64_devB = median(PermEn64_devB.rsurr,2); % median across all surrogates
assert(all(sub==oldsub2),'Subjects loaded in different order')

load(sprintf('PermEn64_neonate_IAAFT1_standardA%s%s',tstr,datestr))
mspe64_stdA = median(PermEn64_stdA.rsurr,2); % median across all surrogates
assert(all(sub==oldsub3),'Subjects loaded in different order')

load(sprintf('PermEn64_neonate_IAAFT1_standardB%s%s',tstr,datestr))
mspe64_stdB = median(PermEn64_stdB.rsurr,2); % median across all surrogates
assert(all(sub==oldsub4),'Subjects loaded in different order')

% MSE

load(sprintf('MSE_neonate_IAAFT1_deviantA%s%s',tstr,datestr))
msmse_devA = squeeze(median(mean(MSE_devA.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_devA  = squeeze(median(MSE_devA.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==oldsub1),'Subjects loaded in different order')

load(sprintf('MSE_neonate_IAAFT1_deviantB%s%s',tstr,datestr))
msmse_devB = squeeze(median(mean(MSE_devB.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_devB  = squeeze(median(MSE_devB.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==oldsub2),'Subjects loaded in different order')

load(sprintf('MSE_neonate_IAAFT1_standardA%s%s',tstr,datestr))
msmse_stdA = squeeze(median(mean(MSE_stdA.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_stdA  = squeeze(median(MSE_stdA.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==oldsub3),'Subjects loaded in different order')

load(sprintf('MSE_neonate_IAAFT1_standardB%s%s',tstr,datestr))
msmse_stdB = squeeze(median(mean(MSE_stdB.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_stdB  = squeeze(median(MSE_stdB.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==oldsub4),'Subjects loaded in different order')

% Add stochaticity test
datestr = '_08-Apr-2023';
load(sprintf('Sto_neonate_10HzLP_deviantA%s_IAAFT2%s',tstr,datestr))
assert(all(sub==oldsub1),'Subjects loaded in a different order')

load(sprintf('Sto_neonate_10HzLP_deviantB%s_IAAFT2%s',tstr,datestr))
assert(all(sub==oldsub2),'Subjects loaded in a different order')

load(sprintf('Sto_neonate_10HzLP_standardA%s_IAAFT2%s',tstr,datestr))
assert(all(sub==oldsub3),'Subjects loaded in a different order')

load(sprintf('Sto_neonate_10HzLP_standardB%s_IAAFT2%s',tstr,datestr))
assert(all(sub==oldsub4),'Subjects loaded in a different order')


%% Time-freq transform

% This code performs time-frequency analysis on magnetoencephalography (MEG) data and generates corresponding plots. Here's a summary of each major part:
% 
%     Time-Frequency Analysis:
%         Sets up parameters for the time-frequency analysis, including the sampling frequency (fs), frequency of interest (foi_start to foi_end), and window shift (window_shift).
%         Calls the function ro_freq_wavelet_TFT to perform the time-frequency transformation on the data (CTW_devA.erf).
% 
%     Plotting for Condition Stimulus (sssd Rule: ssss):
%         Creates a figure with two subplots.
%         Plots the condition/stimulus signal for 'sssd' in the first subplot.
%         Plots the time-frequency representation in the second subplot, including contour lines and tone markers.
%         Calculates correlation between power and entropy (CTW_devA.trh) for each frequency.
% 
%     Correlation coefficients and p-values are calculated for each frequency and stored in variables R1, Pr1, R2, Pr2, R3, Pr3, R4, Pr4.
% 
%     Repeats the Above Steps for Other Conditions:
%         Repeats the time-frequency analysis, plotting, and correlation calculation for conditions 'ssss, rule: sssd' (CTW_devB.erf), 'ssss, rule: ssss' (CTW_stdA.erf), and 'sssd, rule: sssd' (CTW_stdB.erf).
% 
%     Saving Plots:
%         Saves the generated figures in PNG and SVG formats.


% parameters
cfg = [];
cfg.fsample = fs;
cfg.foi_start = 1;
cfg.foi_end = 15;
cfg.window_shift=0.1;
[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_devA.erf),cfg);



% plot the condition/stimulus signal
myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_devA.erf)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (T)')
title('stimulus: sssd, rule: ssss','fontsize',18)
makefighandsome

% plot the condition/stimulus time-freq representation
subplot(2,1,2)
hold on
TFT1 = interp2(squeeze(mean(log10(pow)))',k);
imagesc(TFT1)
contour(TFT1,10,'LineColor','k');
xax = 1:size(TFT1,2);
tmp=linspace(-200,3000,size(TFT1,2));
xidx = find(mod(round(tmp,-1),600)==0);
tmp2 = 1:size(TFT1,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')
yax = linspace(0,size(TFT1,1),length(foi));
yidx = find(mod(foi,1)==0);
yticks(yax(yidx))
yticklabels(round(foi(yidx)))
ylabel('Frequency (Hz)')
title('stimulus: sssd, rule: ssss','fontsize',18)
% Mark the tones
plot(ones(1,size(TFT1,1)).*xax(xidx(1)),1:size(TFT1,1),'w','linewidth',2)
plot(ones(1,size(TFT1,1)).*xax(xidx(2)),1:size(TFT1,1),'w','linewidth',2)
plot(ones(1,size(TFT1,1)).*xax(xidx(3)),1:size(TFT1,1),'w','linewidth',2)
plot(ones(1,size(TFT1,1)).*xax(xidx(4)),1:size(TFT1,1),'w','linewidth',2)
% mycolorbar
caxis([-30 -27])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_neonatal_devA.png')
print('-dsvg','./Figures/TFT_neonatal_devA.svg')

% correlation between power and entropy
R1 = nan(1,length(foi));
Pr1 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(mean(log10(pow(:,:,ifrq))));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_devA.trh,2)),'spline');
    [r,p] = corr(pw2',mean(squeeze(CTW_devA.trh))');
    R1(ifrq) = r;
    Pr1(ifrq) = p;
end

[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_devB.erf),cfg);


myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_devB.erf)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (T)')
title('stimulus: ssss, rule: sssd','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT2 = interp2(squeeze(mean(log10(pow)))',k);
imagesc(TFT2)
contour(TFT2,10,'LineColor','k');
xax = 1:size(TFT2,2);
tmp=linspace(-200,3000,size(TFT2,2));
xidx = find(mod(round(tmp,-1),600)==0);
tmp2 = 1:size(TFT2,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')
yax = linspace(0,size(TFT2,1),length(foi));
yidx = find(mod(foi,1)==0);
yticks(yax(yidx))
yticklabels(round(foi(yidx)))
ylabel('Frequency (Hz)')
title('stimulus: ssss, rule: sssd','fontsize',18)
% Mark the tones
plot(ones(1,size(TFT2,1)).*xax(xidx(1)),1:size(TFT2,1),'w','linewidth',2)
plot(ones(1,size(TFT2,1)).*xax(xidx(2)),1:size(TFT2,1),'w','linewidth',2)
plot(ones(1,size(TFT2,1)).*xax(xidx(3)),1:size(TFT2,1),'w','linewidth',2)
plot(ones(1,size(TFT2,1)).*xax(xidx(4)),1:size(TFT2,1),'w','linewidth',2)
%mycolorbar
caxis([-30 -27])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_neonatal_devB.png')
print('-dsvg','./Figures/TFT_neonatal_devB.svg')

% correlation between power and entropy
R2 = nan(1,length(foi));
Pr2 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(mean(log10(pow(:,:,ifrq))));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_devB.trh,2)),'spline');
    [r,p] = corr(pw2',mean(squeeze(CTW_devB.trh))');
    R2(ifrq) = r;
    Pr2(ifrq) = p;
end

[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_stdA.erf),cfg);


myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_stdA.erf)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (T)')
title('stimulus: ssss, rule: ssss','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT3 = interp2(squeeze(mean(log10(pow)))',k);
imagesc(TFT3)
contour(TFT3,10,'LineColor','k');
xax = 1:size(TFT3,2);
tmp=linspace(-200,3000,size(TFT3,2));
xidx = find(mod(round(tmp,-1),600)==0);
tmp2 = 1:size(TFT3,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')
yax = linspace(0,size(TFT3,1),length(foi));
yidx = find(mod(foi,1)==0);
yticks(yax(yidx))
yticklabels(round(foi(yidx)))
ylabel('Frequency (Hz)')
title('stimulus: ssss, rule: ssss','fontsize',18)
% Mark the tones
plot(ones(1,size(TFT3,1)).*xax(xidx(1)),1:size(TFT3,1),'w','linewidth',2)
plot(ones(1,size(TFT3,1)).*xax(xidx(2)),1:size(TFT3,1),'w','linewidth',2)
plot(ones(1,size(TFT3,1)).*xax(xidx(3)),1:size(TFT3,1),'w','linewidth',2)
plot(ones(1,size(TFT3,1)).*xax(xidx(4)),1:size(TFT3,1),'w','linewidth',2)
%mycolorbar
caxis([-30 -27])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_neonatal_stdA.png')
print('-dsvg','./Figures/TFT_neonatal_stdA.svg')

% correlation between power and entropy
R3 = nan(1,length(foi));
Pr3 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(mean(log10(pow(:,:,ifrq))));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_stdA.trh,2)),'spline');
    [r,p] = corr(pw2',mean(squeeze(CTW_stdA.trh))');
    R3(ifrq) = r;
    Pr3(ifrq) = p;
end

[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_stdB.erf),cfg);


myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_stdB.erf)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (T)')
title('stimulus: sssd, rule: sssd','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT4 = interp2(squeeze(mean(log10(pow)))',k);
imagesc(TFT4)
contour(TFT4,10,'LineColor','k');
xax = 1:size(TFT4,2);
tmp=linspace(-200,3000,size(TFT4,2));
xidx = find(mod(round(tmp,-1),600)==0);
tmp2 = 1:size(TFT4,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')
yax = linspace(0,size(TFT4,1),length(foi));
yidx = find(mod(foi,1)==0);
yticks(yax(yidx))
yticklabels(round(foi(yidx)))
ylabel('Frequency (Hz)')
title('stimulus: sssd, rule: sssd','fontsize',18)
% Mark the tones
plot(ones(1,size(TFT4,1)).*xax(xidx(1)),1:size(TFT4,1),'w','linewidth',2)
plot(ones(1,size(TFT4,1)).*xax(xidx(2)),1:size(TFT4,1),'w','linewidth',2)
plot(ones(1,size(TFT4,1)).*xax(xidx(3)),1:size(TFT4,1),'w','linewidth',2)
plot(ones(1,size(TFT4,1)).*xax(xidx(4)),1:size(TFT4,1),'w','linewidth',2)
%mycolorbar
caxis([-30 -27])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_neonatal_stdB.png')
print('-dsvg','./Figures/TFT_neonatal_stdB.svg')

% correlation between power and entropy
R4 = nan(1,length(foi));
Pr4 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(mean(log10(pow(:,:,ifrq))));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_stdB.trh,2)),'spline');
    [r,p] = corr(pw2',mean(squeeze(CTW_stdB.trh))');
    R4(ifrq) = r;
    Pr4(ifrq) = p;
end

%% Correlation between power and entropy
% Correlation coefficients
foi_hd = 2.^linspace(0,log2(max(foi)),200);
R1hd = interp1(foi,R1,foi_hd,'spline');
R2hd = interp1(foi,R2,foi_hd,'spline');
R3hd = interp1(foi,R3,foi_hd,'spline');
R4hd = interp1(foi,R4,foi_hd,'spline');

%%% A friendly reminder about the labels ;-) %%%
% DevA = ssss/sssd (global deviant)
% DevB = sssd/ssss (local and global deviant)
% StdA = ssss/ssss (standard)
% StdB = sssd/sssd (local deviant)

R = [R1hd' R2hd' R3hd' R4hd']';
myfigure
plot(log2(foi_hd),R,'linewidth',2)
plot(log2(foi_hd),mean(R),'k','linewidth',4)
xticks(0:4)
xticklabels(2.^[0:4])
xlabel('Frequency (Hz)')
ylabel('Correlation (r)')
legend({'sssd/ssss','ssss/sssd','ssss/ssss','sssd/sssd','Mean'},'fontsize',...
    16,'location','northeastoutside','autoupdate','off')
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k:','linewidth',3)
title('Neonatal log10(power) vs CTW correlation','fontsize',20)
ylim([-0.6 0.6])
yticks([-0.6:0.2:0.6])
xlim([0 log2(max(foi))])
legend boxoff
makefighandsome
print('-dpng','./Figures/CTW_neonatal_Power_R.png')
print('-dsvg','./Figures/CTW_neonatal_Power_R.svg')


%% Build table

% Assign values to variables based on experimental conditions, sort the 
% table, split it into smaller tables, and reassemble it.

Datatable = T(strcmp(T.time_window,'800-1000ms'),:); % the variables we are interested are static within a trial
Datatable.Stimulus = nan(size(Datatable,1),1);
Datatable.Violation = nan(size(Datatable,1),1);

% the stimulus is ssss
idx1 = strcmp(Datatable.condition,'a_standard') | strcmp(Datatable.condition,'global_deviant');
Datatable.Stimulus(idx1) = 0; % ssss

% the stimulus is sssd
idx2 = strcmp(Datatable.condition,'local_deviant') | strcmp(Datatable.condition,'local_and_global_deviant');
Datatable.Stimulus(idx2) = 1; % sssd

idx1 = contains(Datatable.condition,'global'); % global deviant
Datatable.Violation(idx1) = 1; % rule violation
idx2 = ~idx1;
Datatable.Violation(idx2) = 0; % congruent with rule

% First sort by violation, then by stimulus
Datatable = sortrows(Datatable,{'Violation','Stimulus'});

% split into four smaller tables

TstdA = Datatable(strcmp(Datatable.condition,'a_standard'),:);
TstdB = Datatable(strcmp(Datatable.condition,'local_deviant'),:);
TdevA = Datatable(strcmp(Datatable.condition,'global_deviant'),:);
TdevB = Datatable(strcmp(Datatable.condition,'local_and_global_deviant'),:);

% sort each smaller table by subject ID

sub1(sub1==0) = 9; % somehow 9 gets misread as 0
sub2(sub2==0) = 9; % somehow 9 gets misread as 0
sub3(sub3==0) = 9; % somehow 9 gets misread as 0
sub4(sub4==0) = 9; % somehow 9 gets misread as 0


idx1 = sortref(TdevA.ID,sub1');
TdevA = TdevA(idx1,:);

idx2 = sortref(TdevB.ID,sub2');
TdevB = TdevB(idx2,:);

idx3 = sortref(TstdA.ID,sub3');
TstdA = TstdA(idx3,:);

idx4 = sortref(TstdB.ID,sub4');
TstdB = TstdB(idx4,:);


% Put the big table back together again
Datatable = [TstdA; TstdB; TdevA; TdevB];
SIDs_in = [sub3'; sub4'; sub1'; sub2'];
assert(all([Datatable.ID] == SIDs_in),'Subject IDs don''t match')

Datatable.CTW = [CTW_stdA.full'; CTW_stdB.full'; CTW_devA.full'; CTW_devB.full'];
Datatable.LZC = [LZC_stdA.full'; LZC_stdB.full'; LZC_devA.full'; LZC_devB.full'];
Datatable.PermEn32 = [PermEn32_stdA.full'; PermEn32_stdB.full'; PermEn32_devA.full'; PermEn32_devB.full'];
Datatable.PermEn64 = [PermEn64_stdA.full'; PermEn64_stdB.full'; PermEn64_devA.full'; PermEn64_devB.full'];
Datatable.MSE = [mean(squeeze(MSE_stdA.full))'; mean(squeeze(MSE_stdB.full))'; ...
    mean(squeeze(MSE_devA.full))'; mean(squeeze(MSE_devB.full))'];
Datatable.SE  = [squeeze(MSE_stdA.full(1,:))'; squeeze(MSE_stdB.full(1,:))'; ...
    squeeze(MSE_devA.full(1,:))'; squeeze(MSE_devB.full(1,:))'];
Datatable.Sto  = [Sto_stdA.full'; Sto_stdB.full'; Sto_devA.full'; Sto_devB.full'];

Datatable.CTWtrn= [CTW_stdA.trunc'; CTW_stdB.trunc'; CTW_devA.trunc'; CTW_devB.trunc'];
Datatable.LZCtrn = [LZC_stdA.trunc'; LZC_stdB.trunc'; LZC_devA.trunc'; LZC_devB.trunc'];
Datatable.PermEn32trn = [PermEn32_stdA.trunc'; PermEn32_stdB.trunc'; PermEn32_devA.trunc'; PermEn32_devB.trunc'];
Datatable.PermEn64trn = [PermEn64_stdA.trunc'; PermEn64_stdB.trunc'; PermEn64_devA.trunc'; PermEn64_devB.trunc'];
Datatable.MSEtrn = [mean(squeeze(MSE_stdA.trunc))'; mean(squeeze(MSE_stdB.trunc))'; ...
    mean(squeeze(MSE_devA.trunc))'; mean(squeeze(MSE_devB.trunc))'];
Datatable.SEtrn  = [squeeze(MSE_stdA.trunc(1,:))'; squeeze(MSE_stdB.trunc(1,:))'; ...
    squeeze(MSE_devA.trunc(1,:))'; squeeze(MSE_devB.trunc(1,:))'];

Datatable.CTWsbt= [CTW_stdA.subtr'; CTW_stdB.subtr'; CTW_devA.subtr'; CTW_devB.subtr'];
Datatable.LZCsbt = [LZC_stdA.subtr'; LZC_stdB.subtr'; LZC_devA.subtr'; LZC_devB.subtr'];
Datatable.PermEn32sbt = [PermEn32_stdA.subtr'; PermEn32_stdB.subtr'; PermEn32_devA.subtr'; PermEn32_devB.subtr'];
Datatable.PermEn64sbt = [PermEn64_stdA.subtr'; PermEn64_stdB.subtr'; PermEn64_devA.subtr'; PermEn64_devB.subtr'];
Datatable.MSEsbt = [mean(squeeze(MSE_stdA.subtr))'; mean(squeeze(MSE_stdB.subtr))'; ...
    mean(squeeze(MSE_devA.subtr))'; mean(squeeze(MSE_devB.subtr))'];
Datatable.SEsbt  = [squeeze(MSE_stdA.subtr(1,:))'; squeeze(MSE_stdB.subtr(1,:))'; ...
    squeeze(MSE_devA.subtr(1,:))'; squeeze(MSE_devB.subtr(1,:))'];
Datatable.Stosbt  = [Sto_stdA.subtr'; Sto_stdB.subtr'; Sto_devA.subtr'; Sto_devB.subtr'];

% In the manuscript we log-scaled HRV
Datatable.SDNN = log10(Datatable.hrv_vec); %log-10 scale hrv
Datatable.HR = log10(Datatable.hr_vec); %log-10 scale heart rate


% Add sex data
Datatable.Sex = cell(size(Datatable,1),1);
Datatable.XChrom = nan(size(Datatable,1),1);
T2 = readtable('./DataManuscript/Exploratoryanalysis.csv');
BMI = readtable('./BMI.xls');

for irow = 1:size(Datatable,1)
    for jrow = 1:size(T2,1)
        if T2.ID(jrow) == Datatable.ID(irow)
            switch T2.gender{jrow}
                case 'm'
                    Datatable.Sex{irow} = 'M';
                    Datatable.XChrom(irow) = 1;
                case 'f'
                    Datatable.Sex{irow} = 'F';
                    Datatable.XChrom(irow) = 2;
                otherwise
                    error('Can''t read participant sex data')
            end
        end
    end

    % Also add maternal BMI as a covariate
    for jrow = 1:size(BMI,1)
        tmpid = str2double(BMI.Pat_ID{jrow}(2:end));
        if tmpid == Datatable.ID(irow)
            Datatable.BMI(irow) = BMI.BMIVorSS(jrow);
        end
    end
end



%% Relationships between measures

    % Calculate and visualize the correlation matrix of entropy measures.
    % Generate a heatmap of the correlation matrix and save it as PNG and SVG files.

% Correlation matrix of entropy measures
r = corr([Datatable.LZC Datatable.CTW Datatable.SE Datatable.MSE Datatable.PermEn32 Datatable.PermEn64]);
myfigure2
mypcolor(rot90(r))
xticks([1.5 2.5 3.5 4.5 5.5 6.5])
xticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
yticks([1.5 2.5 3.5 4.5 5.5 6.5])
yticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
makefigpretty
chrome = customcolormap_preset('red-white-blue');
colormap(chrome)
caxis([-1 1])
mycolorbar
print('-dpng','./Figures/neonate_EntropyMeasuresRMatrix.png')
print('-dsvg','./Figures/neonate_EntropyMeasuresRMatrix.svg')

%% Compute between subjects correlation between entropy measures and power

%     Perform wavelet analysis on MEG data for different experimental conditions.
%     Compute correlations between entropy measures and log-transformed power.
%     Create a new table for power and saves it to a CSV file.
%     Visualize and save the correlation matrix between complexity measures and log-transformed power.

cfg.oct_bw = 1;

for isub = 1:size(CTW_devA.erf,1)

    if isub == 1
        [~,~,~,~,~,~,foi] = ro_freq_wavelet(squeeze(CTW_devA.erf(1,isub,:))',cfg);  % to make sure we get foi is defined
        powDevA = nan(size(CTW_devA.erf,1),length(foi));
        powDevB = nan(size(CTW_devB.erf,1),length(foi));
        powStdA = nan(size(CTW_stdA.erf,1),length(foi));
        powStdB = nan(size(CTW_stdB.erf,1),length(foi));
    end
    powDevA(isub,:) = ro_freq_wavelet(squeeze(CTW_devA.erf(isub,:)),cfg);
end

for isub = 1:size(CTW_devB.erf,1)
    powDevB(isub,:) = ro_freq_wavelet(squeeze(CTW_devB.erf(isub,:)),cfg);
end

for isub = 1:size(CTW_stdA.erf,1)
    powStdA(isub,:) = ro_freq_wavelet(squeeze(CTW_stdA.erf(isub,:)),cfg);
end

for isub = 1:size(CTW_stdB.erf,1)
    powStdB(isub,:) = ro_freq_wavelet(squeeze(CTW_stdB.erf(isub,:)),cfg);
end

Powertable = Datatable; % copy table

Powertable.power = log10([powStdA; powStdB; powDevA; powDevB]);
writetable(Powertable, 'NeonetalEntropyPower.csv')


P = readtable('NeonetalEntropyPower.csv'); % trick to get columns seprated

cmpx = {'LZC','CTW','SE','MSE','PermEn32','PermEn64'};
R = nan(length(cmpx),length(foi));
for i = 1:length(cmpx)
    for j = 1:length(foi)
        eval(sprintf('R(i,j) = corr(P.%s,P.power_%i);',cmpx{i},j))
    end
end

%
myfigure2
mypcolor(flipud(R))
title('Neonatal multicolinearity','fontsize',18)
xticks([1:16]+0.5)
xlim([1 17])
set(gca,'xticklabel',num2str(foi','%.2f'))
%xticklabels(round(foi,2))
xlabel('log_{10}(power) by frequency (Hz)')
yticks([1:length(cmpx)]+0.5)
yticklabels(cmpx)
ylabel('Complexity measure')
chrome = customcolormap_preset('red-white-blue');
colormap(chrome)
mycolorbar
caxis([-1 1])
makefigpretty
print('-dpng','./Figures/Neonatal_multicolinearity.png')
print('-dsvg','./Figures/Neonatal_multicolinearity.svg')


%% Add age in days

%     Read an external table containing age information.
%     Add age data to the main data table based on subject IDs.

T2 = readtable('MeasurementQualityTable.csv');
Datatable.age = nan(size(Datatable,1),1);
IDs = unique(Datatable.ID);

for iid = 1:length(IDs)
    idstr = sprintf('C%03.0f\n',IDs(iid));
    idstr = idstr(1:4); % somehow this is necessary for contains to work
    findme = find(contains(T2.filename,idstr));
    AGE = unique(T2.age_days_(findme));
    assert(length(AGE)==1,'More than one age found')
    Datatable.age(Datatable.ID == IDs(iid)) = AGE;
end


%% Supplemental table

   % Create a supplemental table with information about age, experimental blocks, and hierarchical decomposition.
   % Save the supplemental table as a CSV file.

[~,~,simpid] = unique(T.ID);

T.simpid = simpid;

Tsup = table();
Tsup.SID = [min(simpid):max(simpid)]';
Tsup.age = nan(size(Tsup,1),1);
Tsup.BlockA = nan(size(Tsup,1),1);
Tsup.BlockB = nan(size(Tsup,1),1);
Tsup.HDecomp = nan(size(Tsup,1),1);

for irow = 1:size(Tsup,1)
    id = unique(T.ID(T.simpid==Tsup.SID(irow)));
    Tsup.age(irow) = unique(Datatable.age(Datatable.ID==id));
    if any(contains(Datatable.condition(Datatable.ID==id),'a_condition')) |...
            any(strcmp(Datatable.condition(Datatable.ID==id),'local_and_global_deviant'))
        Tsup.BlockA(irow) = true;
    else
        Tsup.BlockA(irow)= false;
    end
    if any(strcmp(Datatable.condition(Datatable.ID==id),'local_deviant')) |...
            any(strcmp(Datatable.condition(Datatable.ID==id),'global_deviant'))
        Tsup.BlockB(irow) = true;
        strcmp(Datatable.condition(Datatable.ID==id),'global_deviant')
        if any(strcmp(Datatable.condition(Datatable.ID==id),'global_deviant'))
            Tsup.HDecomp(irow) = true;
        else
            Tsup.HDecomp(irow) = false;
        end
    else
        Tsup.BlockB(irow)= false;
        Tsup.HDecomp(irow) = false;
    end
end

Tsup

writetable(Tsup,'./SupplementalNewbornTable.csv')


%% Linear mixed model

   % Compute the rule variable based on stimulus and violation information.
   % Initialize variables for statistical outputs, including p-values, t-stats, and confidence intervals.

fprintf('Running LMMs ...\n\n\n\n')
% Get rule from stim and violation info
Datatable.Rule = xor(Datatable.Stimulus, Datatable.Violation);

measures = {'CTW','LZC','PermEn32','PermEn64','MSE','SE'};
PvalsAGE = nan(1,length(measures));
PvalsSDNN = nan(1,length(measures));
PvalsSXR = nan(1,length(measures));
tStatsAGE = nan(1,length(measures));
tStatsSDNN = nan(1,length(measures));
tStatsSXR = nan(1,length(measures));
tStatsSur = nan(1,length(measures));
PvalsSur = nan(1,length(measures));
QvalsSur = nan(1,length(measures));

mm = cell(6,1);

%% Update table with demographic data
    % Read demographic data from Excel file.
    % Add relevant information such as birth age, mom's age, and birth weight to the main data table.
Tdemo = readtable('./DataPaper/Overview_Participants.xlsx','Sheet',2);

% add some info (birth age, mom's age)
for irow = 1:size(Datatable,1)
    Tidx = find(Tdemo.Var1 == Datatable.ID(irow));
    if isempty(Tidx), keyboard,end
    try
        tmp = sscanf(Tdemo.Var10{Tidx},'%i+%i');
        GAAB = tmp(1)*7 + tmp(2)'; % gestational age at birth (in days)
    catch 
        GAAB = nan;
    end

    Datatable.BirthAge(irow) = GAAB;
    Datatable.MomsAge(irow) = Tdemo.Var8(Tidx);
    Datatable.BirthWeight(irow) = Tdemo.Var11(Tidx);
end

%% stepwise modeling

   % Construct an overall model formula based on predictors used for at 
   % least half of the entropy measures.
   % Run the final linear mixed model and extracts relevant statistics.

for i = 1:length(measures)
    % Stepwise model without random effects
    mdl = stepwiselm(Datatable,sprintf('%s ~ 1',measures{i}),'ResponseVar',...
        measures{i},'PredictorVars',{'age','SDNN','Sex','Rule',...
        'Stimulus','BirthAge','MomsAge','BirthWeight','BMI'},...
        'CategoricalVar',{'Sex','Rule','Stimulus'});
    %sprintf('Using this formula: %s\n',mdl.Formula)
    formula1 = char(mdl.Formula);
    start = strfind(formula1,'~');
    mm{i,1} = char(formula1(start+1:end));
end

allf = strcat(mm{1},mm{2},mm{3},mm{4},mm{5},mm{6});
C = strsplit(allf,' + ');
tok = unique(C);
tok = strip(tok,'1');
tok = strip(tok);
rmv = [];
for i = 1:length(tok)
    if isempty(tok{i})
        rmv = [rmv i];
    end
end
tok(rmv) = [];
tokens = setdiff(tok,{'(1|ID)','(1|ID)~ 1','~ 1'});
Ntok = nan(1,length(tokens));
for itok = 1:length(tokens)
    Ntok(itok) = sum(contains(C,tokens{itok}));
end

mdlstr = '';
cnt = 0;
for itok = 1:length(tokens)
    if Ntok(itok) >= length(measures)/2 % if this predictor is present for at least half the entropy features
        cnt = cnt + 1;
        if cnt == 1
            mdlstr = sprintf('%s ~ %s',mdlstr,tokens{itok});
        else
            mdlstr = sprintf('%s + %s',mdlstr,tokens{itok});
        end
    end
end
mdlstr = sprintf('%s + (1|ID)',mdlstr)

PvalsAGE = nan(6,1);
tStatAGE = nan(6,1);
betaAGE = nan(6,3);
CRIT = 0.024; % The critical P-value determined elsewhere by an FDR correction

% run the modal model
for i = 1:length(measures)
    lme_formula = sprintf('%s %s',measures{i},mdlstr);
    out = fitlme(Datatable,lme_formula);
    PvalsAGE(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'age'));
    tStatAGE(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'age'));
    betaAGE(i,1) = out.Coefficients.Estimate(strcmp(out.Coefficients.Name,'age'));
    betaAGE(i,2) = out.Coefficients.Lower(strcmp(out.Coefficients.Name,'age'));
    betaAGE(i,3) = out.Coefficients.Upper(strcmp(out.Coefficients.Name,'age'));
end

%% Build table from model results

    % Format and export the results of the linear mixed models, including p-values, t-stats, and confidence intervals.
    % Convert the output table to LaTeX format for better presentation.

Tout = table(measures',tStatAGE,PvalsAGE,PvalsAGE < CRIT, betaAGE(:,1),betaAGE(:,2),betaAGE(:,3))
writetable(Tout, 'NeonatalModel.csv')

tmp = table2cell(Tout);
Tltx = cell2table( table2cell(Tout)' );
Tltx.Measure = {'Measure','Age t-stat','Age P-value','Age survive FDR','Age \beta est.','Age \beta LB','Age \beta UB'}';
Tltx = Tltx(:,[end 1:end-1])
Tltx.Properties.VariableNames = Tltx{1,:};
Tltx(1,:) = []
for irow = 1:size(Tltx,1)
    for icol = 2:size(Tltx,2)
        try
            Tltx{irow,icol}{1} = round(Tltx{irow,icol}{1},3,'significant');
        catch
           if islogical(Tltx{irow,icol}{1})
               if  Tltx{irow,icol}{1}
                    Tltx{irow,icol}{1} = 'TRUE';
               else
                    Tltx{irow,icol}{1} = 'FALSE';
               end
           end
        end
    end
end


table2latex(Tltx, 'NeonatalModel')

%% Test whether entropy is noise using surrogate data

   % Run models to test effect of surrogacy on entropy
   % Construct a table (Toutsur) containing statistics for surrogate data.
   % Write the results to a CSV file (NeonatalSurrogate.csv).

test = 'mixedmodels';

Datatable.msctw = [msctw_stdA; msctw_stdB; msctw_devA; msctw_devB];
Datatable.mslzc = [mslzc_stdA; mslzc_stdB; mslzc_devA; mslzc_devB];
Datatable.mspe32 = [mspe32_stdA; mspe32_stdB; mspe32_devA; mspe32_devB];
Datatable.mspe64 = [mspe64_stdA; mspe64_stdB; mspe64_devA; mspe64_devB];
Datatable.msmse = [msmse_stdA'; msmse_stdB'; msmse_devA'; msmse_devB'];
Datatable.msse = [msse_stdA'; msse_stdB'; msse_devA'; msse_devB'];


Biggertable = [Datatable; Datatable]; % stack two big tables
Biggertable.LZCtrn = [Datatable.LZCtrn; Datatable.mslzc];
Biggertable.CTWtrn = [Datatable.CTWtrn; Datatable.msctw];
Biggertable.PermEN32trn = [Datatable.PermEn32trn; Datatable.mspe32];
Biggertable.PermEN64trn = [Datatable.PermEn64trn; Datatable.mspe64];
Biggertable.MSEtrn = [Datatable.MSEtrn; Datatable.msmse];
Biggertable.SEtrn  = [Datatable.SEtrn; Datatable.msse];
Biggertable.Surrogate = logical([zeros(size(Datatable,1),1); ones(size(Datatable,1),1)]);

for i = 1:length(measures)
    eval(sprintf('model = fitlme(Biggertable,''%strn ~ Surrogate + (1|ID)'')',measures{i}))
    SURidx = find(strcmp(model.Coefficients.Name,'Surrogate_1'));
    PvalsSur(i)   = model.Coefficients.pValue(SURidx);
    tStatsSur(i)   = model.Coefficients.tStat(SURidx);
end

QvalsSur = mafdr(PvalsSur,'bhfdr',1);
Toutsur = table(measures',tStatsSur',PvalsSur',...
    'variablenames',{'Measure','Sur_tstat','Sur_Pvalue'})
writetable(Toutsur, 'NeonatalSurrogate.csv')


%% Test relationship between dynamics and age

    % Perform logistic regression to test the relationship between stochastic behavior (IsSto) and age.
    % Create a table (Tout) containing statistics for age in relation to stochastic behavior.
    % Write the results to a CSV file (NeonatalModelDyanmics.csv).

Pvalsage = nan;
PvalsSDNN = nan;
tStatsage = nan;
tStatsSDNN = nan;

Pvalsage2 = nan;
PvalsSDNN2 = nan;
tStatsage2 = nan;
tStatsSDNN2 = nan;

Datatable.IsSto = strcmp(Datatable.Sto,'stochastic');


model = fitglme(Datatable,'IsSto ~ age  + (1|ID)','distribution','binomial','link','logit')
ageidx = find(strcmp(model.Coefficients.Name,'age'));

Pvalsage    = model.Coefficients.pValue(ageidx);
tStatsage   = model.Coefficients.tStat(ageidx);


Tout = table([tStatsage], [Pvalsage],'variablenames',{'age_tstat','age_Pvalue'})
writetable(Tout, 'NeonatalModelDyanmics.csv')

%%  Make histograms of deterministic vs stochastic datasets

   % Generate histograms of age distribution for different conditions 
   % (e.g., stochastic vs. deterministic) and violation types.
   % Save the histograms as SVG files in the Figures directory.

myfigure2
age_sto_std = Datatable.age(strcmp(Datatable.Sto,'stochastic') & ~Datatable.Violation);
age_dtm_std = Datatable.age(strcmp(Datatable.Sto,'deterministic') & ~Datatable.Violation);
histogram(age_sto_std,[12:4:60],'facecolor','b')
histogram(age_dtm_std,[12:4:60],'facecolor','r')

N = length(age_sto_std) + length(age_dtm_std);
stoN = length(age_sto_std);
dtmN = length(age_dtm_std);

title(sprintf('Global standards (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
legend({'Stochastic','Deterministic'},'fontsize',16)
legend box off

title(sprintf('Global standards (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
xlabel('Age (days)')
ylabel('Count (datasets)')
makefighandsome
print('-dsvg','./Figures/NeonatalDynamicsHistogramStdFull')

myfigure2
age_sto_dev = Datatable.age(strcmp(Datatable.Sto,'stochastic') & Datatable.Violation);
age_dtm_dev = Datatable.age(strcmp(Datatable.Sto,'deterministic') & Datatable.Violation);
histogram(age_sto_dev,12:4:60,'facecolor','b')
histogram(age_dtm_dev,12:4:60,'facecolor','r')

N = length(age_sto_dev) + length(age_dtm_dev);
stoN = length(age_sto_dev);
dtmN = length(age_dtm_dev);

title(sprintf('Global deviants (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
legend({'Stochastic','Deterministic'},'fontsize',16)
legend box off

xlabel('Age (days)')
ylabel('Count (datasets)')
makefighandsome
print('-dsvg','./Figures/NeonatalDynamicsHistogramDevFull')

myfigure2
age_sto_std = Datatable.age(strcmp(Datatable.Stosbt,'stochastic') & ~Datatable.Violation);
age_dtm_std = Datatable.age(strcmp(Datatable.Stosbt,'deterministic') & ~Datatable.Violation);
histogram(age_sto_std,12:4:60,'facecolor','b')
histogram(age_dtm_std,12:4:60,'facecolor','r')

N = length(age_sto_std) + length(age_dtm_std);
stoN = length(age_sto_std);
dtmN = length(age_dtm_std);

title(sprintf('Global standards w/subtrac. (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
legend({'Stochastic','Deterministic'},'fontsize',16)
legend box off

xlabel('Age (days)')
ylabel('Count (datasets)')
makefighandsome
print('-dsvg','./Figures/NeonatalDynamicsHistogramStdSubtrFull')

myfigure2
age_sto_dev = Datatable.age(strcmp(Datatable.Stosbt,'stochastic') & Datatable.Violation);
age_dtm_dev = Datatable.age(strcmp(Datatable.Stosbt,'deterministic') & Datatable.Violation);
histogram(age_sto_dev,12:4:60,'facecolor','b')
histogram(age_dtm_dev,12:4:60,'facecolor','r')

N = length(age_sto_dev) + length(age_dtm_dev);
stoN = length(age_sto_dev);
dtmN = length(age_dtm_dev);

title(sprintf('Global deviants w/subtrac. (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
legend({'Stochastic','Deterministic'},'fontsize',16)
legend box off

xlabel('Age (days)')
ylabel('Count (datasets)')
makefighandsome
print('-dsvg','./Figures/NeonatalDynamicsHistogramDevSubtrFull')

% Save the data table

writetable(Datatable,'./NeonatalTable.csv')


%% Time-frequency results %%

% This code performs threshold-free cluster enhancement (TFCE) on the results of linear mixed-effects models applied to MEG data.

%     Setup:
%         Set random number generator seed (rng(111)).
%         Define parameters such as the number of permutations (Nperm), minimum value for lasso regression (lambda_min), and whether to allow interactions (interact).
%         Set up frequency transforms using the ro_freq_wavelet_TFT function.
% 
%     Data Preprocessing:
%         Apply frequency transforms to obtain time-frequency representations (powDevA, powStdA, powDevB, powStdB).
%         Ensure the time-frequency representations have matching sizes.
% 
%     Create a Table for Predictors:
%         Create a table (Tm) with predictors such as age, standard deviation of NN intervals (SDNN), sex, rule, stimulus, BMI, birth age, and birth weight.
% 
%     Lasso Regression:
%         Perform lasso regression using lassoglm with cross-validation to select the optimal regularization parameter.
%         Extract the selected independent variables for the model.
% 
%     Linear Mixed-Effects Model:
%         Fit a linear mixed-effects model for each time-frequency point.
%         Extract p-values and t-statistics for specific predictors (e.g., MomsAge, BirthAge, etc.).
% 
%     Threshold-Free Cluster Enhancement (TFCE):
%         Apply TFCE to the p-values and t-statistics for each predictor.
%         Separate the results into positive and negative clusters.
% 
%     Permutations:
%         Perform permutations to assess the significance of the clusters.
%         Save the results of the permutations.
% 
%     Post-Permutation Analysis:
%         Sum the number of positive and negative clusters for each predictor across permutations.


rng(111) % seed the random number generator again for the cross-validation (in case we load old data and start from here)

Nperm = 200; % how many permutations to run
lambda_min = 0.033; 
interact = false; % allow interaction terms? 

% frequency transforms
cfg = [];
cfg.fsample = fs;
cfg.foi_start = 1;
cfg.foi_end = 15.5;
cfg.window_shift=0.05;

% ssss
[powDevA,~,~,~,~,~,foi] = ro_freq_wavelet_TFT(squeeze(CTW_devA.erf),cfg);
powStdA = ro_freq_wavelet_TFT(squeeze(CTW_stdA.erf),cfg);
assert(all( size(powDevA) == size(powStdA)),'Time-frequency representations don''t match in size')

% sssd
powDevB = ro_freq_wavelet_TFT(squeeze(CTW_devB.erf),cfg);
powStdB = ro_freq_wavelet_TFT(squeeze(CTW_stdB.erf),cfg);
assert(all( size(powDevB) == size(powStdB)),'Time-frequency representations don''t match in size')

powDevA = log10(powDevA);
powDevB = log10(powDevB);
powStdA = log10(powStdA);
powStdB = log10(powStdB);


Tm = table(Datatable.ID,Datatable.age,Datatable.SDNN,Datatable.Stimulus,...
    Datatable.Rule,Datatable.XChrom,Datatable.MomsAge,Datatable.BMI, ...
    Datatable.BirthAge,Datatable.BirthWeight,...
   'VariableNames',{'ID','Age','SDNN','Stimulus','Rule','Sex','MomsAge','BMI','BirthAge','BirthWeight'});

kfold = 10;
POW = [powStdA; powStdB; powDevA; powDevB];
P = POW(:);
Tm2 = repmat(Tm,size(POW,2)*size(POW,3),1);
predictors = [Tm2.Age Tm2.SDNN Tm2.Sex Tm2.Rule Tm2.Stimulus Tm2.BirthAge Tm2.BirthWeight Tm2.MomsAge Tm2.BMI];
terms = {'Age','SDNN','Sex','Rule','Stimulus','BirthAge','BirthWeight','MomsAge','BMI'};

switch interact
    case true
        % compute interaction terms
        X = nan(size(predictors,1), size(predictors,2)^2);
        newterms = cell(1,size(predictors,2)^2);
        counter = 1;
        for i = 1:size(predictors,2)
            for j = 1:size(predictors,2)
                if i > j
                    X(:, counter) = predictors(:, i) .* predictors(:, j);
                    newterms{counter} = sprintf('%s:%s',terms{i},terms{j});
                elseif i == j
                    X(:,counter) = predictors(:,i); 
                    newterms{counter} = terms{i};
                end
                counter = counter + 1;
            end
        end
        
        rmv = all(isnan(X),1);
        X(:,rmv) = [];
        newterms(rmv) = [];
        predictors = X;

    case false
        newterms = terms;
end

L = logspace(log10(lambda_min),0,25);
[B0,FitInfo] =lassoglm(predictors, P,'normal','CV',kfold,'lambda',L,'Alpha',1);
[~,useme] = min(FitInfo.Deviance); % find optimal regularization
lambda = FitInfo.Lambda(useme); % regularization parameter
[B, FitInfo] = lassoglm(predictors, P, 'normal', 'Alpha', 1, 'Lambda', lambda);

IVs = newterms(find(abs(B)>0)) % select independent variables for our model

formula = '1';
for ivar = 1:length(IVs)
    formula = sprintf('%s + %s',formula,IVs{ivar});
end

B

%% Run linear mixed models

% Initialization:
% 
%     Several matrices (PvalsBirthAge, PvalsAgeMA, ..., tStatBMI) are initialized with NaN values. These matrices will store p-values and t-statistics for different variables and model parameters.
% 
% Nested Loop Iteration:
% 
%     The code uses nested loops to iterate over the time and frequency dimensions of a 3D matrix POW.
%     For each combination of itime (time) and ifreq (frequency), a linear mixed-effects model (lme) is fitted using the fitlme function.
% 
% Model Formulation:
% 
%     The formula for the LME model is dynamically generated using the lme_formula variable. It includes a response variable (Power_%i), fixed effects (formula), and a random effect (1|ID).
% 
% Coefficient Analysis:
% 
%     The code checks for specific terms in the fitted model (e.g., 'MomsAge', 'BirthAge', etc.) using ismember.
%     If a term is present, the corresponding p-value and t-statistic are extracted and stored in the respective matrices (e.g., PvalsMomsAge, tStatMomsAge).
% 
% Progress Tracking:
% 
%     The code prints the completion percentage of the nested loop to track the progress of the analysis.
% 
% Threshold-Free Cluster Enhancement (TFCE):
% 
%     After obtaining p-values and t-statistics, the code applies the TFCE method separately for positive and negative values for each variable.
% 
% Saving Results:
% 
%     The results are saved in variables which contain the TFCE results for positive values, and similarly for negative values.
%     The entire process's execution time is measured using tic and toc.
% 
% Saving to File:
% 
%     The final results are saved to a file named based on the current date (fetal_meg_stats_progress_%s).

PvalsBirthAge = nan(size(POW,2),size(POW,3));
PvalsAgeMA = nan(size(POW,2),size(POW,3));
PvalsMABA = nan(size(POW,2),size(POW,3));
PvalsMomsAge = nan(size(POW,2),size(POW,3));
PvalsAge = nan(size(POW,2),size(POW,3));
PvalsSex = nan(size(POW,2),size(POW,3));
PvalsSDNN = nan(size(POW,2),size(POW,3));
PvalsRule = nan(size(POW,2),size(POW,3));
PvalsStimulus = nan(size(POW,2),size(POW,3));
PvalsBMI = nan(size(POW,2),size(POW,3));

tStatBirthAge = nan(size(POW,2),size(POW,3));
tStatAgeMA = nan(size(POW,2),size(POW,3));
tStatMABA = nan(size(POW,2),size(POW,3));
tStatMomsAge = nan(size(POW,2),size(POW,3));
tStatAge = nan(size(POW,2),size(POW,3));
tStatSex = nan(size(POW,2),size(POW,3));
tStatSDNN = nan(size(POW,2),size(POW,3));
tStatRule = nan(size(POW,2),size(POW,3));
tStatStimulus = nan(size(POW,2),size(POW,3));
tStatBMI = nan(size(POW,2),size(POW,3));

tic
ncnt = 0;
for itime = 1:size(POW,2)
    for ifreq = 1:size(POW,3) % do outside parfor loop
        eval(sprintf('Tm.Power_%i = [POW(:,itime,ifreq)];',ifreq));

        lme_formula = sprintf('Power_%i ~ %s + (1|ID)',ifreq,formula);
        lmeout = fitlme(Tm,lme_formula);

        if any(ismember(lmeout.Coefficients.Name,'Age:MomsAge'))
            PvalsAgeMA(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Age:MomsAge'));
            tStatAgeMA(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Age:MomsAge'));
        end

        if any(ismember(lmeout.Coefficients.Name,'MomsAge:BirthAge'))
            PvalsMABA(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'MomsAge:BirthAge'));
            tStatMABA(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'MomsAge:BirthAge'));
        end

        if any(ismember(lmeout.Coefficients.Name,'MomsAge'))
            PvalsMomsAge(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'MomsAge'));
            tStatMomsAge(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'MomsAge'));
        end

        if any(ismember(lmeout.Coefficients.Name,'BirthAge'))
            PvalsBirthAge(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'BirthAge'));
            tStatBirthAge(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'BirthAge'));
        end

        if any(ismember(lmeout.Coefficients.Name,'Age'))
            PvalsAge(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Age'));
            tStatAge(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Age'));
        end

        if any(ismember(lmeout.Coefficients.Name,'Sex'))
            PvalsSex(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Sex'));
            tStatSex(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Sex'));
        end

        if any(ismember(lmeout.Coefficients.Name,'SDNN'))
            PvalsSDNN(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'SDNN'));
            tStatSDNN(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'SDNN'));
        end

        if any(ismember(lmeout.Coefficients.Name,'Rule'))
            PvalsRule(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Rule'));
            tStatRule(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Rule'));
        end

        if any(ismember(lmeout.Coefficients.Name,'Stimulus'))
            PvalsStimulus(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Stimulus'));
            tStatStimulus(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Stimulus'));
        end

        if any(ismember(lmeout.Coefficients.Name,'BMI'))
            PvalsBMI(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'BMI'));
            tStatBMI(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'BMI'));
        end

        ncnt = ncnt + 1;
        fprintf('%1.2f%% complete\n',ncnt/( size(POW,2) * size(POW,3) )*100)
    end
end

[~,~,Age_statsum_pos] = TFCE(PvalsAge,tStatAge,'positive',0);
[~,~,AgeMA_statsum_pos] = TFCE(PvalsAgeMA,tStatAgeMA,'positive',0);
[~,~,MABA_statsum_pos] = TFCE(PvalsMABA,tStatMABA,'positive',0);
[~,~,MomsAge_statsum_pos] = TFCE(PvalsMomsAge,tStatMomsAge,'positive',0);
[~,~,BirthAge_statsum_pos] = TFCE(PvalsBirthAge,tStatBirthAge,'positive',0);
[~,~,Sex_statsum_pos] = TFCE(PvalsSex,tStatSex,'positive',0);
[~,~,SDNN_statsum_pos] = TFCE(PvalsSDNN,tStatSDNN,'positive',0);
[~,~,Rule_statsum_pos] = TFCE(PvalsRule,tStatRule,'positive',0);
[~,~,Stimulus_statsum_pos] = TFCE(PvalsStimulus,tStatStimulus,'positive',0);
[~,~,BMI_statsum_pos] = TFCE(PvalsBMI,tStatBMI,'positive',0);


[~,~,Age_statsum_neg] = TFCE(PvalsAge,tStatAge,'negative',0);
[~,~,AgeMA_statsum_neg] = TFCE(PvalsAgeMA,tStatAgeMA,'negative',0);
[~,~,MABA_statsum_neg] = TFCE(PvalsMABA,tStatMABA,'negative',0);
[~,~,MomsAge_statsum_neg] = TFCE(PvalsMomsAge,tStatMomsAge,'negative',0);
[~,~,BirthAge_statsum_neg] = TFCE(PvalsBirthAge,tStatBirthAge,'negative',0);
[~,~,Sex_statsum_neg] = TFCE(PvalsSex,tStatSex,'negative',0);
[~,~,SDNN_statsum_neg] = TFCE(PvalsSDNN,tStatSDNN,'negative',0);
[~,~,Rule_statsum_neg] = TFCE(PvalsRule,tStatRule,'negative',0);
[~,~,Stimulus_statsum_neg] = TFCE(PvalsStimulus,tStatStimulus,'negative',0);
[~,~,BMI_statsum_neg] = TFCE(PvalsBMI,tStatBMI,'negative',0);

toc

saveme = sprintf('neonatal_meg_stats_progress_%s',date);
save(saveme)

%% now run permutations

% This code block performs permutations for threshold-free cluster enhancement (TFCE) for statistical analysis

% Initialization:
% 
%     Matrices are initialized to store the TFCE results for negative and positive values for different variables and model parameters. These matrices represent the results of permutation tests.
% 
% Parallel Pool Setup:
% 
%     The code attempts to set up a new parallel pool using parpool to enable parallel computing. If unsuccessful, it prints a message indicating that the parallel pool is already open.
% 
% Data Preparation:
% 
%     For each combination of itime (time) and ifreq (frequency) in the outer loop, a temporary variable Tm.(sprintf('Power_%i_%i',itime,ifreq)) is created by concatenating different data matrices (powStdA, powStdB, powDevA, powDevB).
% 
% Permutations (parfor Loop):
% 
%     The code uses a parallel loop (parfor) to perform the permutation Nperm times.
%     For each permutation (iperm), a temporary copy Tm_tmp is created to avoid issues with parallel execution.
%     Inside the permutation loop, for each combination of itime and ifreq, the response variable is permuted using randperm and a new linear mixed-effects model (lme) is fitted.
%     P-values and t-statistics for various coefficients are stored in matrices.
% 
% TFCE Application:
% 
%     TFCE is applied separately for negative and positive values for each variable and permutation. This involves enhancing clusters of significant values in the permutation results.
% 
% Results Storage:
% 
%     The TFCE results for negative and positive values are stored in the
%     corresponding matrices.
% 
% Progress Tracking:
% 
%     The code prints the completion percentage of the permutation loop to track the progress of the analysis.
% 
% Results Saving:
% 
%     After completing all permutations, the results are saved to a file with a name based on the current date (fetal_perm_stats_TFCEs_%s).



tic

Age_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
AgeMA_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
MABA_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
MomsAge_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
BirthAge_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
Sex_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
SDNN_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
Rule_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
Stimulus_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
BMI_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));


Age_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
AgeMA_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
MABA_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
MomsAge_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
BirthAge_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
Sex_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
SDNN_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
Rule_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
Stimulus_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
BMI_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));

try
    parpool
catch
    fprintf('Looks like parpool is already open\n')
end

% Do this outside the parfor loop below
for itime = 1:size(POW,2)
    for ifreq = 1:size(POW,3)
        tmp = [powStdA(:,itime,ifreq); powStdB(:,itime,ifreq); powDevA(:,itime,ifreq); powDevB(:,itime,ifreq)];
        Tm.(sprintf('Power_%i_%i',itime,ifreq)) = tmp;
    end
end

parfor iperm = 1:Nperm
    %ncnt = 0;
    PvalsMABAPerm = nan(size(POW,2),size(POW,3));
    PvalsAgeMAPerm = nan(size(POW,2),size(POW,3));
    PvalsBirthAgePerm = nan(size(POW,2),size(POW,3));
    PvalsMomsAgePerm  = nan(size(POW,2),size(POW,3));
    PvalsAgePerm       = nan(size(POW,2),size(POW,3));
    PvalsSexPerm = nan(size(POW,2),size(POW,3));
    PvalsSDNNPerm = nan(size(POW,2),size(POW,3));
    PvalsRulePerm  = nan(size(POW,2),size(POW,3));
    PvalsStimulusPerm       = nan(size(POW,2),size(POW,3));
    PvalsBMIPerm       = nan(size(POW,2),size(POW,3));

    tStatMABAPerm = nan(size(POW,2),size(POW,3));
    tStatAgeMAPerm = nan(size(POW,2),size(POW,3));
    tStatBirthAgePerm = nan(size(POW,2),size(POW,3));
    tStatMomsAgePerm  = nan(size(POW,2),size(POW,3));
    tStatAgePerm       = nan(size(POW,2),size(POW,3));
    tStatSexPerm = nan(size(POW,2),size(POW,3));
    tStatSDNNPerm = nan(size(POW,2),size(POW,3));
    tStatRulePerm  = nan(size(POW,2),size(POW,3));
    tStatStimulusPerm       = nan(size(POW,2),size(POW,3));
    tStatBMIPerm       = nan(size(POW,2),size(POW,3));

    Tm_tmp = Tm; % Create a temporary copy of Tm for this iteration (only way to get it work with parfor)

    for itime = 1:size(POW,2)
        for ifreq = 1:size(POW,3)
            tmp = Tm_tmp.(sprintf('Power_%i_%i',itime,ifreq));
            NDX = randperm(size(POW,1));
            Tm_tmp.(sprintf('Power_%i_%i',itime,ifreq)) = tmp(NDX);
            lme_formula = sprintf('Power_%i_%i ~ %s + (1|ID)',itime,ifreq,formula);
            lmeout = fitlme(Tm_tmp,lme_formula);


            if any(ismember(lmeout.Coefficients.Name,'MomsAge'))
                PvalsMomsAgePerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'MomsAge'));
                tStatMomsAgePerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'MomsAge'));
            end

            if any(ismember(lmeout.Coefficients.Name,'BirthAge'))
                PvalsBirthAgePerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'BirthAge'));
                tStatBirthAgePerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'BirthAge'));
            end


            if any(ismember(lmeout.Coefficients.Name,'Age'))
                PvalsAgePerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Age'));
                tStatAgePerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Age'));
            end

            if any(ismember(lmeout.Coefficients.Name,'Sex'))
                PvalsSexPerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Sex'));
                tStatSexPerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Sex'));
            end

            if any(ismember(lmeout.Coefficients.Name,'SDNN'))
                PvalsSDNNPerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'SDNN'));
                tStatSDNNPerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'SDNN'));
            end

            if any(ismember(lmeout.Coefficients.Name,'Rule'))
                PvalsRulePerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Rule'));
                tStatRulePerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Rule'));
            end

            if any(ismember(lmeout.Coefficients.Name,'Stimulus'))
                PvalsStimulusPerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'Stimulus'));
                tStatStimulusPerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'Stimulus'));
            end

            if any(ismember(lmeout.Coefficients.Name,'BMI'))
                PvalsBMIPerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'BMI'));
                tStatBMIPerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'BMI'));
            end

        end

    end



    [~,~,Age_ssPerm_neg(iperm,:,:)] = TFCE(PvalsAgePerm,tStatAgePerm,'negative',0);
    [~,~,AgeMA_ssPerm_neg(iperm,:,:)] = TFCE(PvalsAgeMAPerm,tStatAgeMAPerm,'negative',0);
    [~,~,MABA_ssPerm_neg(iperm,:,:)] = TFCE(PvalsMABAPerm,tStatMABAPerm,'negative',0);
    [~,~,MomsAge_ssPerm_neg(iperm,:,:)] = TFCE(PvalsMomsAgePerm,tStatMomsAgePerm,'negative',0);
    [~,~,BirthAge_ssPerm_neg(iperm,:,:)] = TFCE(PvalsBirthAgePerm,tStatBirthAgePerm,'negative',0);
    [~,~,Sex_ssPerm_neg(iperm,:,:)] = TFCE(PvalsSexPerm,tStatSexPerm,'negative',0);
    [~,~,SDNN_ssPerm_neg(iperm,:,:)] = TFCE(PvalsSDNNPerm,tStatSDNNPerm,'negative',0);
    [~,~,Rule_ssPerm_neg(iperm,:,:)] = TFCE(PvalsRulePerm,tStatRulePerm,'negative',0);
    [~,~,Stimulus_ssPerm_neg(iperm,:,:)] = TFCE(PvalsStimulusPerm,tStatStimulusPerm,'negative',0);
    [~,~,BMI_ssPerm_neg(iperm,:,:)] = TFCE(PvalsBMIPerm,tStatBMIPerm,'negative',0)

    [~,~,Age_ssPerm_pos(iperm,:,:)] = TFCE(PvalsAgePerm,tStatAgePerm,'positive',0);
    [~,~,AgeMA_ssPerm_pos(iperm,:,:)] = TFCE(PvalsAgeMAPerm,tStatAgeMAPerm,'positive',0);
    [~,~,MABA_ssPerm_pos(iperm,:,:)] = TFCE(PvalsMABAPerm,tStatMABAPerm,'positive',0);
    [~,~,MomsAge_ssPerm_pos(iperm,:,:)] = TFCE(PvalsMomsAgePerm,tStatMomsAgePerm,'positive',0);
    [~,~,BirthAge_ssPerm_pos(iperm,:,:)] = TFCE(PvalsBirthAgePerm,tStatBirthAgePerm,'positive',0);
    [~,~,Sex_ssPerm_pos(iperm,:,:)] = TFCE(PvalsSexPerm,tStatSexPerm,'positive',0);
    [~,~,SDNN_ssPerm_pos(iperm,:,:)] = TFCE(PvalsSDNNPerm,tStatSDNNPerm,'positive',0);
    [~,~,Rule_ssPerm_pos(iperm,:,:)] = TFCE(PvalsRulePerm,tStatRulePerm,'positive',0);
    [~,~,Stimulus_ssPerm_pos(iperm,:,:)] = TFCE(PvalsStimulusPerm,tStatStimulusPerm,'positive',0);
    [~,~,BMI_ssPerm_pos(iperm,:,:)] = TFCE(PvalsBMIPerm,tStatBMIPerm,'positive',0);
    
    fprintf('%1.2f%% complete\n',iperm/Nperm*100)

end

toc

saveme = sprintf('neonatal_perm_stats_TFCEs_%s',date);
save(saveme)

%% Gather results

% Results Aggregation:
% 
%     Several matrices are initialized to store the count of significant clusters for each variable and permutation.
%     The matrices are initialized with zeros and will be incremented based on the significance of the TFCE clusters.
% 
% Counting Significant Clusters:
% 
%     Two separate loops are used to iterate over each permutation (iperm).
%     For each permutation, the code compares the TFCE results for each variable with the corresponding permuted results.
%     The count matrices are incremented based on whether a cluster in the original data is more extreme than the clusters in the permuted data.


AGEposN = zeros(size(POW,2),size(POW,3));
MAposN = zeros(size(POW,2),size(POW,3));
BAposN = zeros(size(POW,2),size(POW,3));
SEXposN = zeros(size(POW,2),size(POW,3));
SDNNposN = zeros(size(POW,2),size(POW,3));
STIMposN = zeros(size(POW,2),size(POW,3));
RuleposN = zeros(size(POW,2),size(POW,3));
MABAposN = zeros(size(POW,2),size(POW,3));
BMIposN = zeros(size(POW,2),size(POW,3));

AGEnegN = zeros(size(POW,2),size(POW,3));
MAnegN = zeros(size(POW,2),size(POW,3));
BAnegN = zeros(size(POW,2),size(POW,3));
SEXnegN = zeros(size(POW,2),size(POW,3));
SDNNnegN = zeros(size(POW,2),size(POW,3));
STIMnegN = zeros(size(POW,2),size(POW,3));
RulenegN = zeros(size(POW,2),size(POW,3));
MABAnegN = zeros(size(POW,2),size(POW,3));
BMInegN = zeros(size(POW,2),size(POW,3));

for iperm = 1:Nperm
    AGEposN  = AGEposN +  (Age_statsum_pos > squeeze(Age_ssPerm_pos(iperm,:,:)));
    MAposN   = MAposN +  (MomsAge_statsum_pos > squeeze(MomsAge_ssPerm_pos(iperm,:,:)));
    BAposN   = BAposN +  (BirthAge_statsum_pos > squeeze(BirthAge_ssPerm_pos(iperm,:,:)));
    SEXposN  = SEXposN +  (Sex_statsum_pos > squeeze(Sex_ssPerm_pos(iperm,:,:)));
    SDNNposN = SDNNposN +  (SDNN_statsum_pos > squeeze(SDNN_ssPerm_pos(iperm,:,:)));
    STIMposN = STIMposN +  (Stimulus_statsum_pos > squeeze(Stimulus_ssPerm_pos(iperm,:,:)));
    BMIposN  = BMIposN +  (BMI_statsum_pos > squeeze(BMI_ssPerm_pos(iperm,:,:)));
    RuleposN = RuleposN +  (Rule_statsum_pos > squeeze(Rule_ssPerm_pos(iperm,:,:)));
    MABAposN = MABAposN +  (MABA_statsum_pos > squeeze(MABA_ssPerm_pos(iperm,:,:)));
end


for iperm = 1:Nperm
    AGEnegN = AGEnegN +  (Age_statsum_neg > squeeze(Age_ssPerm_neg(iperm,:,:)));
    MAnegN = MAnegN +  (MomsAge_statsum_neg > squeeze(MomsAge_ssPerm_neg(iperm,:,:)));
    BAnegN = BAnegN +  (BirthAge_statsum_neg > squeeze(BirthAge_ssPerm_neg(iperm,:,:)));
    SEXnegN = SEXnegN +  (Sex_statsum_neg > squeeze(Sex_ssPerm_neg(iperm,:,:)));
    SDNNnegN = SDNNnegN +  (SDNN_statsum_neg > squeeze(SDNN_ssPerm_neg(iperm,:,:)));
    STIMnegN = STIMnegN +  (Stimulus_statsum_neg > squeeze(Stimulus_ssPerm_neg(iperm,:,:)));
    BMInegN  = BMInegN +  (BMI_statsum_neg > squeeze(BMI_ssPerm_neg(iperm,:,:)));
    RulenegN = RulenegN +  (Rule_statsum_neg > squeeze(Rule_ssPerm_neg(iperm,:,:)));
    MABAnegN = MABAnegN +  (MABA_statsum_neg > squeeze(MABA_ssPerm_neg(iperm,:,:)));
end


%% Visualization preparation

%     After counting the significant clusters, the code computes p-values by comparing the counts to the total number of permutations.
%     The p-values are transformed to the negative logarithm (-log10) for better visualization.
%     Several variables store the transformed p-values for positive and negative clusters.


tmp = (Nperm - BAposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
BAposPval = -log10(tmp);


tmp = (Nperm - SDNNposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
SDNNposPval = -log10(tmp);


tmp = (Nperm - AGEposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
AGEposPval = -log10(tmp);


tmp = (Nperm - MAposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
MAposPval = -log10(tmp);


tmp = (Nperm - SEXposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
SEXposPval = -log10(tmp);


tmp = (Nperm - STIMposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
STIMposPval = -log10(tmp);


tmp = (Nperm - RuleposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
RuleposPval = -log10(tmp);

tmp = (Nperm - MABAposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
MABAposPval = -log10(tmp);

tmp = (Nperm - BMIposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
BMIposPval = -log10(tmp);

%%%%% Negative

tmp = (Nperm - BAnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
BAnegPval = -log10(tmp);


tmp = (Nperm - SDNNnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);


tmp = (Nperm - AGEnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
AGEnegPval = -log10(tmp);


tmp = (Nperm - MAnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
MAnegPval = -log10(tmp);


tmp = (Nperm - SEXnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
SEXnegPval = -log10(tmp);


tmp = (Nperm - STIMnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
STIMnegPval = -log10(tmp);


tmp = (Nperm - RulenegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
RulenegPval = -log10(tmp);


tmp = (Nperm - MABAnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
MABAnegPval = -log10(tmp);

tmp = (Nperm - BMInegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
BMInegPval = -log10(tmp);



%% clusters

% Set up parameters such as interpolation style, significance threshold (Pcrit), and a custom red-white-blue colormap.

interpstyle = 'linear';
Pcrit = 0.05; % adjust this later as needed
RWB = customcolormap_preset('red-white-blue');

%%% For each fgures ...
%     Plot the relevant map of P-values in time-frequency space
%     Use a custom colormap, adjusts color axis, and adds contours based on statistical significance.
%     Generate and save the figure as an SVG file.

%%%% SDNN/HRV %%%%
myfigure2
plotme1 = flipud(rot90(interp2(SDNNposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(SDNNnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
makefighandsome
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0  log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])

xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
makefighandsome
title('SDNN','fontsize',20)
print('-dsvg','./Figures/SDNNTFCEnewborn')

%%
%%%% Maternal AGE %%%%
myfigure2
plotme1 = flipud(rot90(interp2(MAposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(MAnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])


xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
makefighandsome
title('Maternal age','fontsize',20)
print('-dsvg','./Figures/MATFCEnewborn')


%%
%%%% BMI %%%%
myfigure2
plotme1 = flipud(rot90(interp2(BMIposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(BMInegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])


xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
makefighandsome
title('Maternal BMI','fontsize',20)
print('-dsvg','./Figures/BMITFCEnewborn')

%%
%%% Birth AGE %%%%
myfigure2
plotme1 = flipud(rot90(interp2(BAposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(BAnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])

xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
set(h,'YMinorTick','on')
makefighandsome
title('Birth age','fontsize',20)
print('-dsvg','./Figures/BATFCEnewborn')

%%
%%% AGE %%%%
myfigure2
plotme1 = flipud(rot90(interp2(AGEposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(AGEnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])

xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
set(h,'YMinorTick','on')
makefighandsome
title('Age','fontsize',20)
print('-dsvg','./Figures/AGETFCEnewborn')

%%

%%% MOTHERS AGE %%%%
myfigure2
plotme1 = flipud(rot90(interp2(MAposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(MAnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])

xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
set(h,'YMinorTick','on')
makefighandsome
title('Mother''s age','fontsize',20)
print('-dsvg','./Figures/MATFCEnewborn')


%%
%%%% Sex %%%%
myfigure2
plotme1 = flipud(rot90(interp2(SEXposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(SEXnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])

xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
makefighandsome
title('Sex','fontsize',20)
print('-dsvg','./Figures/SEXTFCEnewborn')



%%%% Stimulus %%%%
myfigure2
plotme1 = flipud(rot90(interp2(STIMposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(STIMnegPval,4,interpstyle)));
plotme = plotme1 + plotme2.*-1;
imagesc(plotme)
caxis([log10(0.005) -log10(0.005)])
mycolorbar;
set(c,'Ticks',sort([log10([0.1 0.05 0.01 0.005]) 0 log10([0.1 0.05 0.01 0.005]).*-1]))
set(c,'TickLabels',[fliplr([0.1 0.05 0.01 0.005]) 0 [0.1 0.05 0.01 0.005]])

xax = 1:size(plotme,2);
tmp=linspace(-200,3000,size(plotme,2));
[~,I1] = min(abs(tmp));
[~,I2] = min(abs(tmp-600));
[~,I3] = min(abs(tmp-1200));
[~,I4] = min(abs(tmp-1800));
[~,I5] = min(abs(tmp-2400));
[~,I6] = min(abs(tmp-3000));
xidx = [I1 I2 I3 I4 I5 I6];

tmp2 = 1:size(plotme,2);
xticks(tmp2(xidx))
xticklabels([0:600:3000])
xlabel('time (ms)')

foihd = interp1(linspace(1,size(plotme,1),length(foi)),foi,1:size(plotme,1));
yt = [];
ytm = [];
target_values = 1:2:13; % Creates a vector [1, 3, 5, 7, 9, 11, 13]
alt_values = 2:2:14; % Creates a vector [2, 4, 6, ..., 14]
rounded_foihd = round(foihd, 1);

for value = target_values
    yt = [yt, find(rounded_foihd == value, 1, 'first')];
end
for value = alt_values
    ytm = [ytm, find(rounded_foihd == value, 1, 'first')];
end

yticks(yt)
yticklabels([1:2:13])
ax = gca;
ax.YAxis.MinorTick = 'on';
ax.YAxis.MinorTickValues = ytm;
ylabel('Frequency (Hz)')

contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
makefighandsome
title('Stimulus','fontsize',20)
print('-dsvg','./Figures/STIMTFCEnewborn')

