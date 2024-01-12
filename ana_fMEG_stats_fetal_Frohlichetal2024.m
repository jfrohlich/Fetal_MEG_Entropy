% This is the script that plots the figures and does the stats
% Joel Frohlich
% University of Tuebingen

clearvars
close all
rng(675)

% Create a directory for the figures
if ~isfolder('Figures'), mkdir('Figures'), end

%%% How to interpret the stimulus/rule labels:
% stim/rule
% DevA = ssss/sssd (global deviant)
% DevB = sssd/ssss (local and global deviant)
% StdA = ssss/ssss (standard)
% StdB = sssd/sssd (local deviant)

% Other abbreviations
% GA = gestational age
% SDNN = standard deviation of the interpeak interval in the MCG (heart rate variability)

% To replicate the manuscript, set this to 'false'
neoLP = false; % load the Neonatal data that were lowpass filtered at 10 Hz to match the fetal filtering settings

load grand_average_ERFs ERFdev ERFstd ID % load ERFs from older fetuses

% Take the average of ERF templates, averging first within subjects
ERFdev = nestedmean(ERFdev,ID);
ERFstd = nestedmean(ERFstd,ID);

T = readtable('dataset_differencesT3T4.csv');
T2 = readtable('Overview_Participants.xlsx');
portion = 'all';
fs = 610.3516; % sampling rate
k = 3; % "spatial" (i.e., time x freq) interpolation factor for time-freq representations

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

datestr = '_24-Jun-2022';

% Add CTW

load(sprintf('CTW_deviantA%s_IAAFT1%s',tstr,datestr))
sub1 = sub;
msctw_devA = median(CTW_devA.rsurr,2); % median across all surrogates

load(sprintf('CTW_deviantB%s_IAAFT1%s',tstr,datestr))
sub2 = sub;
msctw_devB = median(CTW_devB.rsurr,2); % median across all surrogates

load(sprintf('CTW_standardA%s_IAAFT1%s',tstr,datestr))
sub3 = sub;
msctw_stdA = median(CTW_stdA.rsurr,2); % median across all surrogates

load(sprintf('CTW_standardB%s_IAAFT1%s',tstr,datestr))
sub4 = sub;
msctw_stdB = median(CTW_stdB.rsurr,2); % median across all surrogates

% add mean signal value

mval_devA = mean(CTW_devA.erf,2);
mval_devB = mean(CTW_devB.erf,2); 
mval_stdA = mean(CTW_stdA.erf,2); 
mval_stdB = mean(CTW_stdB.erf,2); 


% Now also add LZC

load(sprintf('LZC_deviantA%s_IAAFT1%s',tstr,datestr))
subdevA = sub;
mslzc_devA = median(LZC_devA.rsurr,2); % median across all surrogates
assert(all(sub==sub1),'Subjects loaded in a different order')

load(sprintf('LZC_deviantB%s_IAAFT1%s',tstr,datestr))
subdevB = sub;
mslzc_devB = median(LZC_devB.rsurr,2); % median across all surrogates
assert(all(sub==sub2),'Subjects loaded in a different order')

load(sprintf('LZC_standardA%s_IAAFT1%s',tstr,datestr))
substdA = sub;
mslzc_stdA = median(LZC_stdA.rsurr,2); % median across all surrogates
assert(all(sub==sub3),'Subjects loaded in a different order')

load(sprintf('LZC_standardB%s_IAAFT1%s',tstr,datestr))
substdB = sub;
mslzc_stdB = median(LZC_stdB.rsurr,2); % median across all surrogates
assert(all(sub==sub4),'Subjects loaded in a different order')


% Now also add MSE

load(sprintf('MSE_deviantA%s_IAAFT1%s',tstr,datestr))
msmse_devA = squeeze(median(mean(MSE_devA.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_devA  = squeeze(median(MSE_devA.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==sub1),'Subjects loaded in a different order')

load(sprintf('MSE_deviantB%s_IAAFT1%s',tstr,datestr))
msmse_devB = squeeze(median(mean(MSE_devB.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_devB  = squeeze(median(MSE_devB.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==sub2),'Subjects loaded in a different order')

load(sprintf('MSE_standardA%s_IAAFT1%s',tstr,datestr))
msmse_stdA = squeeze(median(mean(MSE_stdA.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_stdA  = squeeze(median(MSE_stdA.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==sub3),'Subjects loaded in a different order')

load(sprintf('MSE_standardB%s_IAAFT1%s',tstr,datestr))
msmse_stdB = squeeze(median(mean(MSE_stdB.rsurr,1),3)); % median across all surrogates (20 timescale average)
msse_stdB  = squeeze(median(MSE_stdB.rsurr(1,:,:),3)); % median across all surrogates (1st timescale only)
assert(all(sub==sub4),'Subjects loaded in a different order')

% Add stochaticity test

load(sprintf('Sto_deviantA%s_IAAFT2%s',tstr,datestr))
assert(all(sub==sub1),'Subjects loaded in a different order')

load(sprintf('Sto_deviantB%s_IAAFT2%s',tstr,datestr))
assert(all(sub==sub2),'Subjects loaded in a different order')

load(sprintf('Sto_standardA%s_IAAFT2%s',tstr,datestr))
assert(all(sub==sub3),'Subjects loaded in a different order')

load(sprintf('Sto_standardB%s_IAAFT2%s',tstr,datestr))
assert(all(sub==sub4),'Subjects loaded in a different order')

% Add permen32

load(sprintf('PermEn32_deviantA%s_IAAFT1%s',tstr,datestr))
mspe32_devA = median(PermEn32_devA.rsurr,2); % median across all surrogates
assert(all(sub==sub1),'Subjects loaded in a different order')

load(sprintf('PermEn32_deviantB%s_IAAFT1%s',tstr,datestr))
mspe32_devB = median(PermEn32_devB.rsurr,2); % median across all surrogates
assert(all(sub==sub2),'Subjects loaded in a different order')

load(sprintf('PermEn32_standardA%s_IAAFT1%s',tstr,datestr))
mspe32_stdA = median(PermEn32_stdA.rsurr,2); % median across all surrogates
assert(all(sub==sub3),'Subjects loaded in a different order')

load(sprintf('PermEn32_standardB%s_IAAFT1%s',tstr,datestr))
mspe32_stdB = median(PermEn32_stdB.rsurr,2); % median across all surrogates
assert(all(sub==sub4),'Subjects loaded in a different order')

% Add permen64

load(sprintf('PermEn64_deviantA%s_IAAFT1%s',tstr,datestr))
mspe64_devA = median(PermEn64_devA.rsurr,2); % median across all surrogates
assert(all(sub==sub1),'Subjects loaded in a different order')

load(sprintf('PermEn64_deviantB%s_IAAFT1%s',tstr,datestr))
mspe64_devB = median(PermEn64_devB.rsurr,2); % median across all surrogates
assert(all(sub==sub2),'Subjects loaded in a different order')

load(sprintf('PermEn64_standardA%s_IAAFT1%s',tstr,datestr))
mspe64_stdA = median(PermEn64_stdA.rsurr,2); % median across all surrogates
assert(all(sub==sub3),'Subjects loaded in a different order')

load(sprintf('PermEn64_standardB%s_IAAFT1%s',tstr,datestr))
mspe64_stdB = median(PermEn64_stdB.rsurr,2); % median across all surrogates
assert(all(sub==sub4),'Subjects loaded in a different order')

% These are the fetuses in ALL conditions
common = intersect( intersect(sub1,sub2), intersect(sub3,sub4) );
fprintf('%i unique subjects entered the analysis\n',length(unique([sub1 sub2 sub3 sub4])))
fprintf('%i subjects common to all conditions entered the analysis\n',length(common))


%% Time-freq transform

% This code performs time-frequency analysis on fetal magnetoencephalography (MEG) data and generates corresponding plots. Here's a summary of each major part:
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
cfg.foi_end = 10;
cfg.window_shift=0.1;
[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_devA.erf),cfg);


% plot the condition/stimulus signal
% sssD
myfigure2
subplot(2,1,1)
hold on
plot(nestedmean(squeeze(CTW_devA.erf),sub1),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: sssd, rule: ssss','fontsize',18)
makefighandsome

% plot the condition/stimulus time-freq representation
subplot(2,1,2)
hold on
TFT1 = interp2(squeeze(nestedmean(log10(pow),sub1))',k,'linear');
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
%mycolorbar
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_devA.png')
print('-dsvg','./Figures/TFT_devA.svg')

% correlation between power and entropy
R1 = nan(1,length(foi));
Pr1 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(nestedmean(log10(pow(:,:,ifrq)),sub1));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_devA.trh,2)),'linear');
    [r,p] = corr(pw2',nestedmean(squeeze(CTW_devA.trh),sub1)');
    R1(ifrq) = r;
    Pr1(ifrq) = p;
end

% sssS
[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_devB.erf),cfg);
myfigure2
subplot(2,1,1)
hold on
plot(nestedmean(squeeze(CTW_devB.erf),sub2),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: ssss, rule: sssd','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT2 = interp2(squeeze(nestedmean(log10(pow),sub2))',k,'linear');
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
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_devB.png')
print('-dsvg','./Figures/TFT_devB.svg')

% correlation between power and entropy
R2 = nan(1,length(foi));
Pr2 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(nestedmean(log10(pow(:,:,ifrq)),sub2));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_devB.trh,2)),'spline');
    [r,p] = corr(pw2',nestedmean(squeeze(CTW_devB.trh),sub2)');
    R2(ifrq) = r;
    Pr2(ifrq) = p;
end

% ssss
[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_stdA.erf),cfg);
myfigure2
subplot(2,1,1)
hold on
plot(nestedmean(squeeze(CTW_stdA.erf),sub3),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: ssss, rule: ssss','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT3 = interp2(squeeze(nestedmean(log10(pow),sub3))',k,'linear');
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
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_stdA.png')
print('-dsvg','./Figures/TFT_stdA.svg')

% correlation between power and entropy
R3 = nan(1,length(foi));
Pr3 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(nestedmean(log10(pow(:,:,ifrq)),sub3));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_stdA.trh,2)),'spline');
    [r,p] = corr(pw2',nestedmean(squeeze(CTW_stdA.trh),sub3)');
    R3(ifrq) = r;
    Pr3(ifrq) = p;
end

% sssd
[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_stdB.erf),cfg);
myfigure2
subplot(2,1,1)
hold on
plot(nestedmean(squeeze(CTW_stdB.erf),sub4),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: sssd, rule: sssd','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT4 = interp2(squeeze(nestedmean(log10(pow),sub4))',k,'linear');
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
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_stdB.png')
print('-dsvg','./Figures/TFT_stdB.svg')

% correlation between power and entropy
R4 = nan(1,length(foi));
Pr4 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(nestedmean(log10(pow(:,:,ifrq)),sub4));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_stdB.trh,2)),'spline');
    [r,p] = corr(pw2',nestedmean(squeeze(CTW_stdB.trh),sub4)');
    R4(ifrq) = r;
    Pr4(ifrq) = p;
end

%% Time-freq transform with ERF template subtraction
% Template subtraction means that an average ERF across subjects is
% subtracted from the signal
% NOTE: Template subtraction was NOT performed in the manuscript! If you
% are just trying to replicate our findings, this block of code is
% extraneous.

cfg = [];
cfg.fsample = fs;
cfg.foi_start = 1;
cfg.foi_end = 10;
cfg.window_shift=0.1;
[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_devA.erf-ERFdev),cfg);

myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_devA.erf)-ERFdev),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: sssd, rule: ssss; ERF subtracted','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT1 = interp2(squeeze(mean(log10(pow)))',k,'linear');
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
%mycolorbar
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_devA_subtr.png')
print('-dsvg','./Figures/TFT_devA_subtr.svg')

% correlation between power and entropy
R1 = nan(1,length(foi));
Pr1 = nan(1,length(foi));
for ifrq = 1:size(pow,3) % for each frequency
    pw = squeeze(mean(log10(pow(:,:,ifrq))));
    pw2 = interp1(linspace(-200,3000,size(pow,2)),pw,linspace(-200,3000,size(CTW_devA.trh,2)),'linear');
    [r,p] = corr(pw2',mean(squeeze(CTW_devA.trh))');
    R1(ifrq) = r;
    Pr1(ifrq) = p;
end

[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_devB.erf-ERFdev),cfg);
myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_devB.erf-ERFdev)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: ssss, rule: sssd, ERF subtracted','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT2 = interp2(squeeze(mean(log10(pow)))',k,'linear');
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
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_devB_subtr.png')
print('-dsvg','./Figures/TFT_devB_subtr.svg')

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

[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_stdA.erf-ERFstd),cfg);
myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_stdA.erf-CTW_stdA.erf-ERFstd)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: ssss, rule: ssss, ERF subtracted','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT3 = interp2(squeeze(mean(log10(pow)))',k,'linear');
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
title('stimulus: ssss, rule: ssss, ERF subtracted','fontsize',18)
% Mark the tones
plot(ones(1,size(TFT3,1)).*xax(xidx(1)),1:size(TFT3,1),'w','linewidth',2)
plot(ones(1,size(TFT3,1)).*xax(xidx(2)),1:size(TFT3,1),'w','linewidth',2)
plot(ones(1,size(TFT3,1)).*xax(xidx(3)),1:size(TFT3,1),'w','linewidth',2)
plot(ones(1,size(TFT3,1)).*xax(xidx(4)),1:size(TFT3,1),'w','linewidth',2)
%mycolorbar
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_stdA_subtr.png')
print('-dsvg','./Figures/TFT_stdA_subtr.svg')

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

[pow,pow_full,n,unit,foi_target,foi_delta,foi,DAT] = ro_freq_wavelet_TFT(squeeze(CTW_stdB.erf-ERFstd),cfg);
myfigure2
subplot(2,1,1)
hold on
plot(mean(squeeze(CTW_stdB.erf-ERFstd)),'b','linewidth',2)
xticks([])
xlabel('time (ms)')
ylabel('Field strength (AU)')
title('stimulus: sssd, rule: sssd, ERF subtracted','fontsize',18)
makefighandsome

subplot(2,1,2)
hold on
TFT4 = interp2(squeeze(mean(log10(pow)))',k,'linear');
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
caxis([1.5  3])
makefighandsome
colormap hot
print('-dpng','./Figures/TFT_stdB_subtr.png')
print('-dsvg','./Figures/TFT_stdB_subtr.svg')

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

%% Build table

load Data_traces

% Count how many fetuses had more than one visit
%     Identifies unique IDs in the dataset (Datatable.ID).
%     Counts the number of fetuses with more than one visit (nages > 1) and prints the count.

IDs = unique(Datatable.ID);
cnt = 0;
for isub = 1:length(IDs)
    idx = Datatable.ID == IDs(isub);
    nages = length(unique(Datatable.GA(idx)));
    if nages > 1, cnt = cnt + 1; end
end

fprintf('%i fetuses gave data at more than one visit\n',cnt)


% Next, we add columns to Datatable for conditions, SDNN (standard deviation of NN intervals), XChrom (sex), DD (difference between global deviant and standard), and BMI (maternal body mass index).

% add SDNN to the main table and create a simpler column with condition labels and add SDNN to the main table
Datatable.condition = nan(size(Datatable,1),1);
Datatable.SDNN = nan(size(Datatable,1),1);

Datatable.XChrom = nan(size(Datatable,1),1);
Datatable.DD = nan(size(Datatable,1),1);
Datatable.BMI = nan(size(Datatable,1),1);
BMI = readtable('./BMI.xls');

for irow = 1:size(Datatable,1)
    % Add condition label (simple to read)
    if strcmp(Datatable.Blocklabel(irow,:),'A (ssss)')
        Datatable.condition(irow) = 1;
    elseif strcmp(Datatable.Blocklabel(irow,:),'B (sssd)')
        Datatable.condition(irow) = 2;
    else
        error('condition not recognized')
    end

    % Add SDNN (note: SDNN is the same for all data from one block)
    for jrow = 1:size(T,1)
        if T.ID(jrow) == Datatable.ID(irow) && T.GA(jrow) == Datatable.GA(irow) ...
                && strcmp(T.Block{jrow},Datatable.Blocklabel(irow,1))
            Datatable.SDNN(irow) = T.SDNN(jrow);
            % Also add the difference (global dev vs. std) of differences (T4 - T3)
            NDX = find(T.ID == Datatable.ID(irow) & T.GA == Datatable.GA(irow) ...
                & strcmp(T.Block,Datatable.Blocklabel(irow,1)) & T.tonepart==2);
            GD = NDX(find(contains(T.condition(NDX),'global_deviant'))); % global deviant
            GS = NDX(find(~contains(T.condition(NDX),'global_deviant'))); % global standard

            Datatable.DD(irow) = T.value(GD) - T.value(GS);

        end
    end

    % Add sex as a covariate
    for jrow = 1:size(T2,1)
        if T2.ID(jrow) == Datatable.ID(irow)
            switch T2.FetalSex{jrow}
                case 'm'
                    Datatable.XChrom(irow) = 1;
                    Datatable.Sex{irow} = 'M';
                    break
                case 'f'
                    Datatable.XChrom(irow) = 2;
                    Datatable.Sex{irow} = 'F';
                    break
                otherwise
                    continue
            end
        end
    end

    % Add maternal BMI as a covariate
    for jrow = 1:size(BMI,1)
        tmpid = str2double(BMI.Pat_ID{jrow}(2:end));
        if tmpid == Datatable.ID(irow)
            Datatable.BMI(irow) = BMI.BMIVorSS(jrow);
        end
    end
end

% Concatenate corresponding values

%%% Mean signal value %%%
Datatable.mval_dev = [mval_devA; mval_devB];
Datatable.mval_std = [mval_stdA; mval_stdB];

%%% CTW %%%
% add CTW
Datatable.CTW_dev = [CTW_devA.full CTW_devB.full]';
Datatable.CTW_std = [CTW_stdA.full CTW_stdB.full]';

% Median CTW surrogate
Datatable.msctw_dev = [msctw_devA; msctw_devB];
Datatable.msctw_std = [msctw_stdA; msctw_stdB];

%%% LZC %%%
% add LZC
Datatable.LZC_dev = [LZC_devA.full LZC_devB.full]';
Datatable.LZC_std = [LZC_stdA.full LZC_stdB.full]';

% median surrogate LZC
Datatable.mslzc_dev = [mslzc_devA; mslzc_devB];
Datatable.mslzc_std = [mslzc_stdA; mslzc_stdB];

% add PermEn32
Datatable.PE32_dev = [PermEn32_devA.full PermEn32_devB.full]';
Datatable.PE32_std = [PermEn32_stdA.full PermEn32_stdB.full]';

% median surrogate PermEn32
Datatable.mspe32_dev = [mspe32_devA; mspe32_devB];
Datatable.mspe32_std = [mspe32_stdA; mspe32_stdB];

% add PermEn64
Datatable.PE64_dev = [PermEn64_devA.full PermEn64_devB.full]';
Datatable.PE64_std = [PermEn64_stdA.full PermEn64_stdB.full]';

% median surrogate PermEn64
Datatable.mspe64_dev = [mspe64_devA; mspe64_devB];
Datatable.mspe64_std = [mspe64_stdA; mspe64_stdB];


%%% MSE %%%
% add MSE (average across n timescales)
Datatable.MSE_dev = [mean(squeeze(MSE_devA.full)) mean(squeeze(MSE_devB.full))]';
Datatable.MSE_std = [mean(squeeze(MSE_stdA.full)) mean(squeeze(MSE_stdB.full))]';

% median surrogate MSE
Datatable.msmse_dev = [msmse_devA msmse_devB]';
Datatable.msmse_std = [msmse_stdA msmse_stdB]';

%%% SampEn %%%
% add SE (i.e., sample entropy, first timescale only)
Datatable.SE_dev = [squeeze(MSE_devA.full(1,:)) squeeze(MSE_devB.full(1,:))]';
Datatable.SE_std = [squeeze(MSE_stdA.full(1,:)) squeeze(MSE_stdB.full(1,:))]';

% median surrogate SE
Datatable.msse_dev = [msse_devA msse_devB]';
Datatable.msse_std = [msse_stdA msse_stdB]';


% truncation (chop off parts of signal to match the surrogate signal)

%%% CTW %%%
% add CTW
Datatable.CTW_devT = [CTW_devA.trunc CTW_devB.trunc]';
Datatable.CTW_stdT = [CTW_stdA.trunc CTW_stdB.trunc]';

%%% LZC %%%
% add LZC
Datatable.LZC_devT = [LZC_devA.trunc LZC_devB.trunc]';
Datatable.LZC_stdT = [LZC_stdA.trunc LZC_stdB.trunc]';

%%% PermEn %%%
% add PermEn32
Datatable.PE32_devT = [PermEn32_devA.trunc PermEn32_devB.trunc]';
Datatable.PE32_stdT = [PermEn32_stdA.trunc PermEn32_stdB.trunc]';

% add PermEn64
Datatable.PE64_devT = [PermEn64_devA.trunc PermEn64_devB.trunc]';
Datatable.PE64_stdT = [PermEn64_stdA.trunc PermEn64_stdB.trunc]';

%%% MSE %%%
% add MSE (average across n timescales)
Datatable.MSE_devT = [mean(squeeze(MSE_devA.trunc)) mean(squeeze(MSE_devB.trunc))]';
Datatable.MSE_stdT = [mean(squeeze(MSE_stdA.trunc)) mean(squeeze(MSE_stdB.trunc))]';

%%% SampEn %%%
% add SE (i.e., sample entropy, first timescale only)
Datatable.SE_devT = [squeeze(MSE_devA.trunc(1,:)) squeeze(MSE_devB.trunc(1,:))]';
Datatable.SE_stdT = [squeeze(MSE_stdA.trunc(1,:)) squeeze(MSE_stdB.trunc(1,:))]';


%%%%%%%%%%%%%% subtraction (subtract out the ERF template) %%%%%%%%%%%%%%%%
% NOTE - This was NOT done in the manuscript! 

%%% CTW %%%
% add CTW
Datatable.CTW_devS = [CTW_devA.subtr CTW_devB.subtr]';
Datatable.CTW_stdS = [CTW_stdA.subtr CTW_stdB.subtr]';

%%% LZC %%%
% add LZC
Datatable.LZC_devS = [LZC_devA.subtr LZC_devB.subtr]';
Datatable.LZC_stdS = [LZC_stdA.subtr LZC_stdB.subtr]';

%%% PermEn %%%
% add PermEn32
Datatable.PE32_devS = [PermEn32_devA.subtr PermEn32_devB.subtr]';
Datatable.PE32_stdS = [PermEn32_stdA.subtr PermEn32_stdB.subtr]';

% add PermEn64
Datatable.PE64_devS = [PermEn64_devA.subtr PermEn64_devB.subtr]';
Datatable.PE64_stdS = [PermEn64_stdA.subtr PermEn64_stdB.subtr]';

%%% MSE %%%
% add MSE (average across n timescales)
Datatable.MSE_devS = [mean(squeeze(MSE_devA.subtr)) mean(squeeze(MSE_devB.subtr))]';
Datatable.MSE_stdS = [mean(squeeze(MSE_stdA.subtr)) mean(squeeze(MSE_stdB.subtr))]';

%%% SampEn %%%
% add SE (i.e., sample entropy, first timescale only)
Datatable.SE_devS = [squeeze(MSE_devA.subtr(1,:)) squeeze(MSE_devB.subtr(1,:))]';
Datatable.SE_stdS = [squeeze(MSE_stdA.subtr(1,:)) squeeze(MSE_stdB.subtr(1,:))]';

%%%%%%%%%%%%% DETERMINISTIC/STOCHASTIC TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Stochasticity %%%
Datatable.Sto_dev = [Sto_devA.full Sto_devB.full]';
Datatable.Sto_std = [Sto_stdA.full Sto_stdB.full]';

% With ERF subtraction
Datatable.Sto_devS = [Sto_devA.subtr Sto_devB.subtr]';
Datatable.Sto_stdS = [Sto_stdA.subtr Sto_stdB.subtr]';

%%%%%%%%%%%%%%%%% Random amplitude PermEn32 %%%%%%%%%%%%%%%%%%%%%%%%%%%

Datatable.dev_amp = cell(size(Datatable,1),1);
Datatable.std_amp = cell(size(Datatable,1),1);
Datatable.dev_phi = cell(size(Datatable,1),1);
Datatable.std_phi = cell(size(Datatable,1),1);
Datatable.dev_ra = cell(size(Datatable,1),1);
Datatable.std_ra = cell(size(Datatable,1),1);

for irow = 1:size(Datatable,1)
    gsn = Datatable.global_standards_normalized{irow};
    Datatable.std_amp{irow} = abs(fft(gsn));
    Datatable.std_phi{irow} = angle(fft(gsn));
    gdn = Datatable.global_deviants_normalized{irow};
    Datatable.dev_amp{irow} = abs(fft(gdn));
    Datatable.dev_phi{irow} = angle(fft(gdn));
end

% Interaction effect: Delta^I (i) = H(A_i, Phi_i) - 0.5*[mean_j H(A_i, Phi_j) + mean_j H(A_j, Phi_i)]
% Amplitude effect: Delta^A (i) = mean_j H(A_i, Phi_j) - mean_{j,k} H(A_j, Phi_k)
% Phase effect: Delta^Phi (i) = mean_j H(A_j, Phi_i) - mean_{j,k} H(A_j, Phi_k)
% "Total" effect: Delta^T (i) = H(A_i, Phi_i) - mean_{j,k} H(A_j, Phi_k)

m = 3;
tol = 1e-9; % define a tolerance for the decomp
lag1 = 32;
tau1 = round((lag1/1000)/(1/fs));
lag2 = 64;
tau2 = round((lag2/1000)/(1/fs));

H = nan(size(Datatable,1)); % matrix of entropies, line of identity gives the "true" entropies
for irow = 1:size(H,1)
    for icol = 1:size(H,2)
        tmp = ifft(Datatable.dev_amp{irow}.*exp(1j.*Datatable.dev_phi{icol}), 'symmetric'); % devaint response
        Hdev(irow,icol) = PermEn(tmp,m,tau1);
        tmp = ifft(Datatable.std_amp{irow}.*exp(1j.*Datatable.std_phi{icol}), 'symmetric'); % standard response
        Hstd(irow,icol) = PermEn(tmp,m,tau1);
    end
end

for isub = 1:size(Datatable,1)
    assert(Datatable.PE32_std(isub) == Hstd(isub,isub),'Entropies don''t match') % sanity check
    assert(Datatable.PE32_dev(isub) == Hdev(isub,isub),'Entropies don''t match') % sanity check

    Datatable.PE32AC_dev(isub) = mean(Hdev(isub,:)) - mean(Hdev,[1 2]); % amplitude contribution
    Datatable.PE32PC_dev(isub) = mean(Hdev(:,isub)) - mean(Hdev,[1 2]); % phase contribution
    Datatable.PE32XC_dev(isub) = Hdev(isub,isub) - 0.5*(mean(Hdev(:,isub)) + mean(Hdev(isub,:))); % interaction contribution
    assert((Hdev(isub,isub) - mean(Hdev,[1 2])) - (Datatable.PE32XC_dev(isub) ...
        + (Datatable.PE32AC_dev(isub) + Datatable.PE32PC_dev(isub))/2) < tol,'Components don''t add up');


    Datatable.PE32AC_std(isub) = mean(Hstd(isub,:)) - mean(Hstd,[1 2]); % amplitude contribution
    Datatable.PE32PC_std(isub) = mean(Hstd(:,isub)) - mean(Hstd,[1 2]); % phase contribution
    Datatable.PE32XC_std(isub) = Hstd(isub,isub) - 0.5*(mean(Hstd(:,isub)) + mean(Hstd(isub,:))); % interaction contribution
    assert((Hstd(isub,isub) - mean(Hstd,[1 2])) - (Datatable.PE32XC_std(isub) ...
        + (Datatable.PE32AC_std(isub) + Datatable.PE32PC_std(isub))/2) < tol,'Components don''t add up');
end

assert(all([Datatable.ID] == [sub1'; sub2']),'Subject IDs don''t match')
assert(all([Datatable.ID] == [sub3'; sub4']),'Subject IDs don''t match')


%% Prepare data for linear mixed models

% This block of code prepares a large table (Bigtable) for a linear mixed model by stacking and organizing various features from the original data table (Datatable). Here's a summary:
%
% 1. Data Stacking:
% 
%     Creates a new table (Bigtable) by stacking two copies of the original data table (Datatable).
% 
% 2. Feature Concatenation:
% 
%     Concatenates various signal features for standards and deviants separately. Features include CTW, LZC, PE32, PE64, MSE, SE, and a binary indicator for stochasticity (Sto).
%     The data is organized to have the first half of rows corresponding to standards and the second half to deviants.
% 
% 3. Additional Features:
% 
%     Concatenates additional features, such as truncated (trn) and template subtracted (sbt) versions of CTW, LZC, PE32, PE64, MSE, and SE.
%     Includes columns for amplitude contributions (pe32ac), phase contributions (pe32pc), and interaction contributions (pe32xc).
% 
% 4. Additional Information:
% 
%     Concatenates mean signal values (mval), a binary indicator for rule violation (Violation), and a variable indicating the global rule (GlobalRule).
%     Creates a binary variable (Stimulus) based on whether there is a violation compared to the global rule.
% 
% 5. Log Transformation:
% 
%     Applies a log transformation to the SDNN (standard deviation of NN intervals) and the maximum peak amplitude (MPA).
% 
% 6. Adding T4 - T3 Difference:
% 
%     Adds a column (T4T3) representing the difference between T4 and T3 values from another data table (T). 
%     NOTE: we did not use this information in the manuscript. 
% 
% 7. Sanity Check:
% 
%     Asserts that the length of unique values in the T4T3 column is equal to the size of Bigtable for a sanity check.
% 
% 8. Removing Unnecessary Columns:
% 
%     Removes unnecessary columns from Bigtable to clean up and streamline the data for further analysis.


Bigtable = [Datatable; Datatable]; % stack tables
Bigtable.T4T3 = nan(size(Bigtable,1),1);

% First half of rows is standards, second half of rows are deviants
Bigtable.CTW = [Datatable.CTW_std; Datatable.CTW_dev];
Bigtable.LZC = [Datatable.LZC_std; Datatable.LZC_dev];
Bigtable.PE32 =  [Datatable.PE32_std; Datatable.PE32_dev];
Bigtable.PE64 =  [Datatable.PE64_std; Datatable.PE64_dev];
Bigtable.MSE = [Datatable.MSE_std; Datatable.MSE_dev];
Bigtable.SE  = [Datatable.SE_std;  Datatable.SE_dev];
Bigtable.Sto = strcmp([Datatable.Sto_std;  Datatable.Sto_dev],'stochastic');

Bigtable.CTWtrn = [Datatable.CTW_stdT; Datatable.CTW_devT];
Bigtable.LZCtrn = [Datatable.LZC_stdT; Datatable.LZC_devT];
Bigtable.PE32trn =  [Datatable.PE32_stdT; Datatable.PE32_devT];
Bigtable.PE64trn =  [Datatable.PE64_stdT; Datatable.PE64_devT];
Bigtable.MSEtrn = [Datatable.MSE_stdT; Datatable.MSE_devT];
Bigtable.SEtrn  = [Datatable.SE_stdT;  Datatable.SE_devT];

Bigtable.CTWsbt = [Datatable.CTW_stdS; Datatable.CTW_devS];
Bigtable.LZCsbt = [Datatable.LZC_stdS; Datatable.LZC_devS];
Bigtable.PE32sbt =  [Datatable.PE32_stdS; Datatable.PE32_devS];
Bigtable.PE64sbt =  [Datatable.PE64_stdS; Datatable.PE64_devS];
Bigtable.MSEsbt = [Datatable.MSE_stdS; Datatable.MSE_devS];
Bigtable.SEsbt  = [Datatable.SE_stdS;  Datatable.SE_devS];
Bigtable.Stosbt = strcmp([Datatable.Sto_stdS;  Datatable.Sto_devS],'stochastic');

Bigtable.mslzc  = [Datatable.mslzc_std;  Datatable.mslzc_dev];
Bigtable.msctw  = [Datatable.msctw_std;  Datatable.msctw_dev];
Bigtable.mspe32   = [Datatable.mspe32_std;  Datatable.mspe32_dev];
Bigtable.mspe64   = [Datatable.mspe64_std;  Datatable.mspe64_dev];
Bigtable.msmse  = [Datatable.msmse_std;  Datatable.msmse_dev];
Bigtable.msse  = [Datatable.msse_std;  Datatable.msse_dev];

Bigtable.pe32ac = [Datatable.PE32AC_std; Datatable.PE32AC_dev];
Bigtable.pe32pc = [Datatable.PE32PC_std; Datatable.PE32PC_dev];
Bigtable.pe32xc = [Datatable.PE32XC_std; Datatable.PE32XC_dev];

Bigtable.mval = [Datatable.mval_std; Datatable.mval_dev];

Bigtable.Violation = [zeros(size(Datatable,1),1); ones(size(Datatable,1),1)];
Bigtable.GlobalRule = Bigtable.condition-1; %0 = ssss block label, 1 = sssd block label
Bigtable.Stimulus = double(Bigtable.Violation ~= Bigtable.GlobalRule);

%log-scale SDNN (better than removing outliers, higher Shapario Wilk stat)
Bigtable.SDNN = log10(Bigtable.SDNN);
% Also produce a column that gives the log-scaled max amplitude
Bigtable.MPA = log10(Bigtable.maximum_peak_amplitude);

% Add T4 - T3 
for irow = 1:size(Bigtable,1)
    for jrow = 1:size(T,1)
        if T.ID(jrow) == Bigtable.ID(irow) & T.GA(jrow) == Bigtable.GA(irow) ...
                & strcmp(T.Block{jrow},Bigtable.Blocklabel(irow,1)) ...
                & T.StimCon(jrow)-1 == Bigtable.Stimulus(irow) & T.tonepart(jrow) == 2 % 350 - 650 ms window
            Bigtable.T4T3(irow) = T.value(jrow);
            break
        end
    end
end

assert(length(unique(Bigtable.T4T3))==size(Bigtable,1)); % sanity check

% Remove unnecessary columns
Bigtable = removevars(Bigtable,{'Blocklabel','all_trials','maximum_peak_amplitude', ...
    'all_trails_normalized','global_standards','global_standards_normalized',...
    'global_deviants','global_deviants_normalized','condition','CTW_dev','CTW_std',...
    'LZC_dev','LZC_std','MSE_dev','MSE_std','SE_dev','SE_std',...
    'CTW_devT','CTW_stdT','LZC_devT','LZC_stdT','MSE_devT','MSE_stdT','SE_devT','SE_stdT',...
    'CTW_devS','CTW_stdS','LZC_devS','LZC_stdS','MSE_devS','MSE_stdS','SE_devS','SE_stdS'});


% Remove NaNs and zscore data
TX = Bigtable(randperm(size(Bigtable,1)),[1:5 51:57]);
rmv = any(isnan(TX{:,:}),2);
TX(rmv,:) = [];
TX.SDNN = zscore(TX.SDNN);
TX.LZC = zscore(TX.LZC);
TX.CTW = zscore(TX.CTW);
TX.MSE = zscore(TX.MSE);
TX.SE = zscore(TX.SE);
TX.PE32 = zscore(TX.PE32);
TX.PE64 = zscore(TX.PE64);
TX.DD = zscore(TX.DD);
writetable(TX,'FetalFeatures.xlsx')


%% How many subjects had usable data (must also include SDNN)

% This code block performs several analyses to assess the impact of different predictors on entropy measures.
% 
% 1. Checking Usable Data:
% 
%     Initializes counters (scnt and rcnt) and an empty array (usable).
%     Iterates through subject IDs (IDs) and checks if there is unique data for each subject with non-NaN values in the SDNN column.
%     Updates counters and the usable array accordingly.
%     Prints the number of recordings and subjects with usable data.
% 
% 2. Adding Demographic Data:
% 
%     Reads demographic data from an Excel file (Overview_Participants.xlsx) into a table (Tdemo).
%     Extracts gestational age at birth (GAAB), mother's age, and birth weight from Tdemo.
%     Updates Bigtable with the extracted demographic data.
% 
% 3. Modeling:
% 
%     Performs stepwise linear mixed-effects modeling (stepwiselm) for each measure (CTW, LZC, PE32, PE64, MSE, SE) with various predictors, including gestational age (GA), standard deviation of NN intervals (SDNN), sex, global rule, stimulus, birth age, mother's age, birth weight, and BMI.
%     Compares models with and without random effects and prints whether random effects significantly improve the model.
% 
% 4. Model Comparison:
% 
%     Calculates p-values for the likelihood ratio test between models with and without random effects.
%     Outputs tables with t-statistics, p-values, and effect size estimates for gestational age, sex, and their interaction for each measure.
%     Writes the results to a CSV file (FetalModel.csv).
%     Converts the table into LaTeX format and writes it to a file (FetalModel.tex).
% 
% 5. Likelihood Ratio Test:
% 
%     Applies a false discovery rate correction to the p-values from the likelihood ratio test.
%     Outputs a table with measures, log-likelihood ratio statistics, and corrected p-values.
%     Writes the results to a CSV file (LogLikelihoodTable.csv).
%     Converts the table into LaTeX format and writes it to a file (LogLikelihoodTableLATEX.tex).

scnt = 0;
rcnt = 0;
usable = [];

for isub = 1:length(IDs)
    if ~isempty(unique(Bigtable.GA(Bigtable.ID==IDs(isub) & ~isnan(Bigtable.SDNN))))
        scnt = scnt + 1;
        usable = [usable IDs(isub)];
    end
    rcnt = rcnt + length(unique(Bigtable.GA(Bigtable.ID==IDs(isub) & ~isnan(Bigtable.SDNN))));
end

fprintf('%i recordings from %i fetal subjects with usable data (including SDNN)\n',rcnt,scnt)

measures = {'CTW','LZC','PE32','PE64','MSE','SE'};
Bigtable.Rule = xor(Bigtable.Stimulus,Bigtable.Violation);

% GA = gestational age
% SDNN = standard deviation of the interpeak interval in the MCG (heart rate variability)

RApvals = nan(length(measures),length(usable));

LRStats = nan(length(measures),1);
LR_pvals = nan(length(measures),1);
mm = cell(6,1);

% Update table with demographic data
Tdemo = readtable('./DataPaper/Overview_Participants.xlsx','Sheet',2);

% Uncomment below if you need to correct "ET" for 'expected time' to 40 weeks gestational age. 
% Tdemo.Var10{strcmp(Tdemo.Var10,'ET')} = '40+0'; % fix coding, 'expected time'

% add some info (birth age, mom's age)
for irow = 1:size(Bigtable,1)
    Tidx = find(Tdemo.Var1 == Bigtable.ID(irow));
    if isempty(Tidx), keyboard, end % stop script if needed for debugging
    try
        tmp = sscanf(Tdemo.Var10{Tidx},'%i+%i');
        GAAB = tmp(1)*7 + tmp(2)'; % gestational age at birth (in days)
    catch
        GAAB = nan;
    end

    Bigtable.BirthAge(irow) = GAAB;
    Bigtable.MomsAge(irow) = Tdemo.Var8(Tidx);
    Bigtable.BirthWeight(irow) = Tdemo.Var11(Tidx);
end

for i = 1:length(measures)
    % Stepwise model without random effects
    mdl = stepwiselm(Bigtable,sprintf('%s ~ 1',measures{i}),'ResponseVar',...
        measures{i},'PredictorVars',{'GA','SDNN','Sex','GlobalRule',...
        'Stimulus','BirthAge','MomsAge','BirthWeight','BMI'},...
        'CategoricalVar',{'Sex','GlobalRule','Stimulus'});
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

mdlstr0 = '';
cnt = 0;
for itok = 1:length(tokens)
    if Ntok(itok) >= length(measures)/2 % if this predictor is present for at least half the entropy features
        cnt = cnt + 1;
        if cnt == 1
            mdlstr0 = sprintf('%s ~ %s',mdlstr0,tokens{itok});
        else
            mdlstr0 = sprintf('%s + %s',mdlstr0,tokens{itok});
        end
    end
end
mdlstr = sprintf('%s + (1|ID)',mdlstr0)


% Uncomment below to use the alternative approach of choosing the modal
% model (NOTE: this is NOT what we did in the manuscript)
% mdlstr = char(mode(categorical(mm),[1 2])); % most common model

PvalsGA = nan(6,1);
PvalsSex = nan(6,1);
PvalsSexGA = nan(6,1);
PvalsMomsAge = nan(6,1);
tStatGA = nan(6,1);
tStatSex = nan(6,1);
tStatSexGA = nan(6,1);
tStatMomsAge = nan(6,1);
betaGA = nan(6,3);
betaSex = nan(6,3);
betaSexGA = nan(6,3);
betaMomsAge = nan(6,3);

% run the modal model
for i = 1:length(measures)
    lme_formula = sprintf('%s %s',measures{i},mdlstr);
    out = fitlme(Bigtable,lme_formula);
    PvalsGA(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'GA'));
    PvalsSex(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'Sex_F'));
    PvalsSexGA(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'GA:Sex_F'));
    PvalsMomsAge(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'MomsAge'));
    tStatGA(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'GA'));
    tStatSex(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'Sex_F'));
    tStatSexGA(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'GA:Sex_F'));
    tStatMomsAge(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'MomsAge'));
    betaGA(i,1) = out.Coefficients.Estimate(strcmp(out.Coefficients.Name,'GA'));
    betaSex(i,1) = out.Coefficients.Estimate(strcmp(out.Coefficients.Name,'Sex_F'));
    betaSexGA(i,1) = out.Coefficients.Estimate(strcmp(out.Coefficients.Name,'GA:Sex_F'));
    betaMomsAge(i,1) = out.Coefficients.Estimate(strcmp(out.Coefficients.Name,'MomsAge'));
    
    betaGA(i,2) = out.Coefficients.Lower(strcmp(out.Coefficients.Name,'GA'));
    betaSex(i,2) = out.Coefficients.Lower(strcmp(out.Coefficients.Name,'Sex_F'));
    betaSexGA(i,2) = out.Coefficients.Lower(strcmp(out.Coefficients.Name,'GA:Sex_F'));
    betaMomsAge(i,2) = out.Coefficients.Lower(strcmp(out.Coefficients.Name,'MomsAge'));

    betaGA(i,3) = out.Coefficients.Upper(strcmp(out.Coefficients.Name,'GA'));
    betaSex(i,3) = out.Coefficients.Upper(strcmp(out.Coefficients.Name,'Sex_F'));
    betaSexGA(i,3) = out.Coefficients.Upper(strcmp(out.Coefficients.Name,'GA:Sex_F'));
    betaMomsAge(i,3) = out.Coefficients.Upper(strcmp(out.Coefficients.Name,'MomsAge'));


    alt_formula = sprintf('%s %s',measures{i},mdlstr0);
    out0 = fitlme(Bigtable,alt_formula);

    if out.LogLikelihood  <  out0.LogLikelihood % if the simpler model without random effects has the higher log likelihood
        fprintf('Random effects don''t improve model for %s\n',measures{i})
        pause(1)
    else
        [results,siminfo] = compare(out0,out,'checknesting',true);
        LRStats(i) = results.LRStat(2);
        LR_pvals(i) = results.pValue(2);
        if results.pValue < 0.05
            fprintf('Random effects significantly (p = %1.1g) improve model for %s\n',LR_pvals(i),measures{i})
            pause(1)
        else
            fprintf('Random effects don''t significantly (p = %1.1g) improve model for %s\n',LR_pvals(i),measures{i})
            pause(1)
        end
    end
end

CRIT = 0.024; % The critical value as determined elsewhere using FDR

Tout = table(measures',tStatGA,PvalsGA,PvalsGA<CRIT,betaGA(:,1),betaGA(:,2),betaGA(:,3),...
    tStatSex,PvalsSex,PvalsSex<CRIT,betaSex(:,1),betaSex(:,2),betaSex(:,3),...
    tStatSexGA,PvalsSexGA,PvalsSexGA<CRIT,betaSexGA(:,1),betaSexGA(:,2),betaSexGA(:,3),...
    tStatMomsAge,PvalsMomsAge,PvalsMomsAge<CRIT,betaMomsAge(:,1),betaMomsAge(:,2),betaMomsAge(:,3))

writetable(Tout, 'FetalModel.csv')

tmp = table2cell(Tout);
Tltx = cell2table( table2cell(Tout)' );
Tltx.Measure = {'Measure','GA t-stat','GA P-value','GA survive FDR','GA \beta est.','GA \beta LB','GA \beta UB',...
    'Sex t-stat','Sex P-value','Sex survive FDR','Sex \beta est.','Sex \beta LB','Sex \beta UB',...
    'Sex x GA t-stat','Sex x GA P-value','Sex x GA survive FDR','Sex x GA \beta est.','Sex x GA \beta LB','Sex x GA \beta UB',...
    'MA t-stat','MA P-value','MA survive FDR','MA \beta est.','MA \beta LB','MA \beta UB'}'
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
table2latex(Tltx, 'FetalModel')


LR_qvals = mafdr(LR_pvals,'bhfdr',1);
Tlr = table(measures',round(LRStats,3,'significant'),round(LR_pvals,3,'significant'),'VariableNames',{'Measure','Log-liklihood ratio stat','Pval_RAW'});
writetable(Tlr,'LogLikelihoodTable.csv')
table2latex(Tlr,'LogLikelihoodTableLATEX')

% Clear unnecessary variables from the workspace.
clear out Tout

%% Test relationship between dynamics (stochastic versus deterministic) and other variables

% This code block fits a logistic regression model to assess the relationship between the dynamics and selected predictors. It extracts relevant statistical measures and produces a summary table for further analysis. The results are saved in both CSV and LaTeX formats.
% 1. Logistic Regression Model:
% 
%     Uses fitglme to perform a logistic regression model (out) with the predictor Sto (dynamics) and other predictors specified by the mdlstr formula.
%     The response variable is binary (Sto) are denotes whether model dynamics are stochastic (true/false).
% 
% 2. Extracting Statistical Measures:
% 
%     Extracts p-values and t-statistics for the predictors (GA, Sex_F, GA:Sex_F, MomsAge) from the logistic regression model.
% 
% 3. Creating a Summary Table:
% 
%     Creates a table (Tout) containing t-statistics and p-values for each predictor.
%     Writes the table to a CSV file (FetalModelDynamics.csv).
%     Converts the table into LaTeX format and writes it to a file (FetalModelDynamicsLATEX.tex).


out = fitglme(Bigtable,sprintf('Sto %s',mdlstr),'distribution','binomial','link','logit')
PvalsGA = nan(1,1);
PvalsSex = nan(1,1);
PvalsSexGA = nan(1,1);
PvalsMomsAge = nan(1,1);
tStatGA = nan(1,1);
tStatSex = nan(1,1);
tStatSexGA = nan(1,1);
tStatMomsAge = nan(1,1);

PvalsGA = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'GA'));
PvalsSex = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'Sex_F'));
PvalsSexGA = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'GA:Sex_F'));
PvalsMomsAge = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'MomsAge'));
tStatGA = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'GA'));
tStatSex = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'Sex_F'));
tStatSexGA = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'GA:Sex_F'));
tStatMomsAge = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'MomsAge'));


Tout = table(tStatGA,PvalsGA,...
    tStatSex,PvalsSex,tStatSexGA,PvalsSexGA,...
    tStatMomsAge,PvalsMomsAge,...
    'variablenames',{'GA_tstat','GA_Pvalue',...
    'Sex_tstat','Sex_Pvalue','SexGA_tstat','SexGA_Pvalue',...
    'MomsAge_tstat','MomsAge_Pvalue'})

writetable(Tout, 'FetalModelDyanmics.csv')
table2latex(Tout, 'FetalModelDyanmicsLATEX')

%% Make histograms of deterministic vs stochastic datasets

% This code block visually compares the distribution of gestational ages
% for global standards and global deviants based on their dynamics. 
% The resulting histograms help visualize differences between in dynamics
% between younger and older subjects.

% 1. Global Standards Histogram:
% 
%     Separates gestational ages into two groups: GA_sto_std (stochastic) and GA_dtm_std (deterministic).
%     Plots histograms for each group with blue and red colors, respectively.
%     Calculates the total count (N), count for stochastic (stoN), and count for deterministic (dtmN).
%     Adds a title indicating the percentages of stochastic and deterministic datasets.
%     Displays a legend and adjusts plot aesthetics.
%     Saves the figure as an SVG file (./Figures/DynamicsHistogramStdFull.svg).
% 
% 2. Global Deviants Histogram:
% 
%     Similar to the first part but focuses on global deviants (GA_sto_dev and GA_dtm_dev).
%     Plots histograms, calculates counts, adds a title, legend, and adjusts aesthetics.
%     Saves the figure as an SVG file (./Figures/DynamicsHistogramDevFull.svg).

myfigure2
GA_sto_std = Datatable.GA(strcmp(Datatable.Sto_std,'stochastic'));
GA_dtm_std = Datatable.GA(strcmp(Datatable.Sto_std,'deterministic'));
histogram(GA_sto_std,'facecolor','b')
histogram(GA_dtm_std,'facecolor','r')

N = length(GA_sto_std) + length(GA_dtm_std);
stoN = length(GA_sto_std);
dtmN = length(GA_dtm_std);

title(sprintf('Global standards (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
legend({'Stochastic','Deterministic'},'fontsize',16)
legend box off

xlabel('Gestational age')
ylabel('Count (datasets)')
makefighandsome
print('-dsvg','./Figures/DynamicsHistogramStdFull')

myfigure2
GA_sto_dev = Datatable.GA(strcmp(Datatable.Sto_dev,'stochastic'));
GA_dtm_dev = Datatable.GA(strcmp(Datatable.Sto_dev,'deterministic'));
histogram(GA_sto_dev,'facecolor','b')
histogram(GA_dtm_dev,'facecolor','r')

N = length(GA_sto_dev) + length(GA_dtm_dev);
stoN = length(GA_sto_dev);
dtmN = length(GA_dtm_dev);

title(sprintf('Global deviants (STO %i / %1.3f%%; DTM %i / %1.3f%%',...
    stoN,stoN/N*100,dtmN,dtmN/N*100),'fontsize',18)
legend({'Stochastic','Deterministic'},'fontsize',16)
legend box off

xlabel('Gestational age')
ylabel('Count (datasets)')
makefighandsome
print('-dsvg','./Figures/DynamicsHistogramDevFull')


%% Correlation between power and entropy

% Correlation Visualization:
% 
%     Interpolates correlation coefficients (R1, R2, R3, R4) at higher density (foi_hd) using spline interpolation.
%     Plots the correlation coefficients against the log2-transformed frequency.
%     Adds a mean line for better interpretation.
%     Sets axis labels, legend, and title.
%     Saves the figure in PNG and SVG formats.

foi_hd = 2.^linspace(0,log2(max(foi)),200);
R1hd = interp1(foi,R1,foi_hd,'spline');
R2hd = interp1(foi,R2,foi_hd,'spline');
R3hd = interp1(foi,R3,foi_hd,'spline');
R4hd = interp1(foi,R4,foi_hd,'spline');

R = [R1hd' R2hd' R3hd' R4hd']';

myfigure2
plot(log2(foi_hd),R,'linewidth',2)
plot(log2(foi_hd),mean(R),'k','linewidth',4)
xticks(0:3)
xticklabels(2.^[0:3])
xlabel('Frequency (Hz)')
ylabel('Correlation (r)')
legend({'sssd/ssss','ssss/sssd','ssss/ssss','sssd/sssd','Mean'},'fontsize',...
    16,'location','northeastoutside','autoupdate','off')
plot(log2(foi_hd),zeros(1,length(foi_hd)),'k:','linewidth',3)
title('Fetal log10(power) vs CTW correlation','fontsize',20)
xlim([0 log2(max(foi))])
ylim([-0.6 0.6])
yticks([-0.6:0.2:0.6])
legend boxoff
makefighandsome
print('-dpng','./Figures/CTW_Power_R.png')
print('-dsvg','./Figures/CTW_Power_R.svg')

%% Compute spectral power from the original signal stored in the CTW output

%  Compute spectral power:
% 
%     Initializes variables for storing power values (powDevA, powDevB, powStdA, powStdB).
%     Computes power values for each subject and condition.
%     Creates a table (P) containing the power values for each frequency.
%     Writes this table to a CSV file (FetalEntropyPower.csv).

cfg.oct_bw = 1;
% this will give an error if foi changes size
powDevA = nan(size(CTW_devA.erf,1),14);
powDevB = nan(size(CTW_devB.erf,1),14);
powStdA = nan(size(CTW_stdA.erf,1),14);
powStdB = nan(size(CTW_stdB.erf,1),14);

for isub = 1:size(CTW_devA.erf,1)
    % to make sure we return foi (frequencies of interest)
    if isub == 1, [~,~,~,~,~,~,foi] = ro_freq_wavelet(squeeze(CTW_devA.erf(isub,:)),cfg); end
    % compute power
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

P = Bigtable; % copy table

for ifreq = 1:size(powStdA,2)
    eval(sprintf('P.power_%i = log10([powStdA(:,%i); powStdB(:,%i); powDevA(:,%i); powDevB(:,%i)]);',...
        ifreq,ifreq,ifreq,ifreq,ifreq))
end

writetable(P, 'FetalEntropyPower.csv')

%% Look at beta coefficients on z-scored data from LMEs (should behave like correlation coefficient)

% Correlation Coefficients on Z-scored Data:
% 
%     Z-scores the data in the table (PZ).
%     For each complexity measure (cmpx) and frequency, fits a linear mixed-effects model (fitlme).
%     Extracts beta coefficients and checks if they are within the expected range [-1, 1].
%     Visualizes the multicollinearity pattern using a heatmap.
%     Saves the heatmap in PNG and SVG formats.

PZ = P; % copy table
for icol = 2:size(PZ,2) % for each column except the first (IDs)
    if isa(PZ{:,icol},'double') || isa(PZ{:,icol},'single') % if it's the appropriate class
        PZ{:,icol} = nanzscore(PZ{:,icol});
       hat assert(round(nanvar(PZ{:,icol}),3)==1,'Column does not have unit variance within given tolerance after z-scoring')
    end
end


cmpx = {'LZC','CTW','SE','MSE','PE32','PE64'};
betas = nan(length(cmpx),length(foi));
for i = 1:length(cmpx)
    for j = 1:length(foi)
        mdl = sprintf('%s ~ power_%i + (1|ID)',cmpx{i},j);
        model = fitlme(PZ,mdl);
        % Take the model beta that is NOT the intercept
        beta = model.Coefficients.Estimate(find(~strcmp(model.Coefficients.Name,'(Intercept)')));
        assert(round(beta,5) <= 1 & round(beta,5) >= -1,'Beta not bounded [-1, 1] inclusive')
        betas(i,j) = beta;
    end
end

myfigure2
mypcolor(flipud(betas))
title('Fetal multicolinearity','fontsize',18)
xlim([1 15])
xticks([1:16]+0.5)
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
print('-dpng','./Figures/Fetal_multicolinearity_betas.png')
print('-dsvg','./Figures/Fetal_multicolinearity_betas.svg')

%% Predict one entropy measure from another to compute betas
% 
% Predict One Entropy Measure from Another:
% 
%     Calculates beta coefficients for predicting one entropy measure from another.
%     Creates a heatmap visualizing these beta coefficients.
%     Saves the heatmap in PNG and SVG formats.
%     Generates additional visualizations showing differences and averages of beta coefficients.


betas2 = nan(length(cmpx),length(cmpx));
for i = 1:length(cmpx)
    for j = 1:length(cmpx)
        mdl = sprintf('%s ~ %s + (1|ID)',cmpx{i},cmpx{j});
        model = fitlme(PZ,mdl);
        % Take the model beta that is NOT the intercept
        beta = model.Coefficients.Estimate(find(~strcmp(model.Coefficients.Name,'(Intercept)')));
        assert(round(beta,5) <= 1 & round(beta,5) >= -1,'Beta not bounded [-1, 1] inclusive')
        betas2(i,j) = beta;
    end
end

plotme = rot90(betas2);
myfigure2
mypcolor(plotme)
xticks([1.5 2.5 3.5 4.5 5.5 6.5])
xticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
yticks([1.5 2.5 3.5 4.5 5.5 6.5])
yticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
makefigpretty
colormap(chrome)
xlim([1 7])
caxis([-1 1])
mycolorbar
print('-dpng','./Figures/EntropyMeasuresBetas.png')
print('-dsvg','./Figures/EntropyMeasuresBetas.svg')

myfigure2
mypcolor(plotme - rot90(flipud(plotme)))
xticks([1.5 2.5 3.5 4.5 5.5 6.5])
xticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
yticks([1.5 2.5 3.5 4.5 5.5 6.5])
yticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
makefigpretty
colormap(chrome)
xlim([1 7])
caxis([-1 1])
mycolorbar
print('-dpng','./Figures/EntropyMeasuresBetas_diff.png')
print('-dsvg','./Figures/EntropyMeasuresBetas_diff.svg')

myfigure2
mypcolor((plotme + rot90(flipud(plotme)))./2)
xticks([1.5 2.5 3.5 4.5 5.5 6.5])
xticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
yticks([1.5 2.5 3.5 4.5 5.5 6.5])
yticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
makefigpretty
colormap(chrome)
xlim([1 7])
caxis([-1 1])
mycolorbar
print('-dpng','./Figures/EntropyMeasuresBetas_avg_sym.png')
print('-dsvg','./Figures/EntropyMeasuresBetas_avg_sym.svg')

% Correlation Matrix of Entropy Measures:
% 
%     Computes the correlation matrix between different entropy measures.
%     Creates a heatmap visualizing the correlation matrix.
%     Saves the heatmap in PNG and SVG formats.
%     Generates additional visualizations showing differences between beta coefficients and correlation matrix.

r = corr([P.LZC P.CTW P.SE P.MSE P.PE32 P.PE64]);
myfigure2
mypcolor(rot90(r))
xticks([1.5 2.5 3.5 4.5 5.5 6.5])
xticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
yticks([1.5 2.5 3.5 4.5 5.5 6.5])
yticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
makefigpretty
colormap(chrome)
xlim([1 7])
caxis([-1 1])
mycolorbar
print('-dpng','./Figures/EntropyMeasuresRMatrix.png')
print('-dsvg','./Figures/EntropyMeasuresRMatrix.svg')

myfigure2
mypcolor(((plotme + rot90(flipud(plotme)))./2) - rot90(r))
xticks([1.5 2.5 3.5 4.5 5.5 6.5])
xticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
yticks([1.5 2.5 3.5 4.5 5.5 6.5])
yticklabels({'LZC','CTW','mSampEn','mMSE','PermEn32','PermEn64'})
makefigpretty
colormap(chrome)
xlim([1 7])
caxis([-1 1])
mycolorbar
print('-dpng','./Figures/EntropyMeasures_beta_minus_r.png')
print('-dsvg','./Figures/EntropyMeasures_beta_minus_r.svg')



%% Visualize correlations between entropy and age for both fetuses and newborns
% 
%     Mean centers data for the same participants.
%     Computes the correlation between entropy measures and age for both fetuses and newborns.

% Mean center data from the same participants
IDs = unique(P.ID);
I = find(contains(properties(P),{'CTW','LZC','SE','MSE','PE','ctw','lzc','se','mse','pe'}));
I = intersect(I,1:size(P,2)); % so that we restrict to properties that are table columns
for isub = 1:length(IDs)
    idx = P.ID == IDs(isub); % indices of data from this participant
    if length(idx) > 1 % only mean center if there is more than one recording (otherwise we just get zero!)
        mu = mean(P{idx,I});
        P{idx,I} = P{idx,I} - mu;
    end
end

cmpx = {'LZC','CTW','SE','MSE','PE32','PE64'};
R = nan(length(cmpx),length(foi));
for i = 1:length(cmpx)
    for j = 1:length(foi)
        eval(sprintf('R(i,j) = corr(P.%s,P.power_%i);',cmpx{i},j))
    end
end

%% Test whether entropy is noise using surrogate data
% 
%     Performs a statistical test to determine if entropy measures are significantly different from surrogate data.
%     Saves the results in a CSV file (FetalSurrogate.csv).

test = 'mixedmodels';
measures = {'CTW','LZC','SE','MSE','PE32','PE64'};

PvalsSur = nan(1,length(measures));
tStatsSur = nan(1,length(measures));

Biggertable = [Bigtable; Bigtable]; % stack two big tables
Biggertable.LZC = [Bigtable.LZCtrn; Bigtable.mslzc];
Biggertable.CTW = [Bigtable.CTWtrn; Bigtable.msctw];
Biggertable.PE32 = [Bigtable.PE32trn; Bigtable.mspe32];
Biggertable.PE64 = [Bigtable.PE64trn; Bigtable.mspe64];
Biggertable.MSE = [Bigtable.MSEtrn; Bigtable.msmse];
Biggertable.SE  = [Bigtable.SEtrn; Bigtable.msse];
Biggertable.Surrogate = logical([zeros(size(Bigtable,1),1); ones(size(Bigtable,1),1)]);
for i = 1:length(measures)
    eval(sprintf('model = fitlme(Biggertable,''%s ~ Surrogate + %s'')',measures{i},mdlstr(3:end)))
    SURidx = find(strcmp(model.Coefficients.Name,'Surrogate_1'));
    PvalsSur(i)   = model.Coefficients.pValue(SURidx);
    tStatsSur(i)   = model.Coefficients.tStat(SURidx);
end

Toutsur = table(measures',tStatsSur',PvalsSur',...
    'variablenames',{'Measure','Sur_tstat','Sur_Pvalue'})
writetable(Toutsur, 'FetalSurrogate.csv')

%% Find who has data pre and post birth and combine data
% 
% NOTE: We did NOT combine data in the manuscript -- it is unclear if
% combining fetal and neonatal data is valid given the very different
% nature of these reocrdings. 
%
%     Reads neonatal data from a CSV file.
%     Combines fetal and neonatal data into a single table (Alltable).
%     Adds information about whether the data is from a fetus or a newborn.

% read neonatal data
switch neoLP
    case false
        Neo = readtable('NeonatalTable.csv');
    case true
        Neo = readtable('NeonatalTable_10HzLP.csv');
end

Alltable = Bigtable;

Alltable.BornYet = zeros(size(Alltable,1),1);


% add newborn data
rowcnt = size(Alltable,1);
for irow = 1:size(Neo,1)
    rowcnt = rowcnt + 1;
    Tidx = find(Tdemo.Var1 == Neo.ID(irow));
    tmp = sscanf(Tdemo.Var10{Tidx},'%i+%i');
    GAAB = tmp(1)*7 + tmp(2)'; % gestational age at birth (in days)
    Alltable.ID(rowcnt) = Neo.ID(irow);
    Alltable.SDNN(rowcnt) = Neo.SDNN(irow);
    Alltable.Sex{rowcnt} = Neo.Sex{irow};
    Alltable.XChrom(rowcnt) = Neo.XChrom(irow);
    Alltable.PE32(rowcnt) = Neo.PermEn32(irow);
    Alltable.PE64(rowcnt) = Neo.PermEn64(irow);
    Alltable.MSE(rowcnt) = Neo.MSE(irow);
    Alltable.SE(rowcnt) = Neo.SE(irow);
    Alltable.LZC(rowcnt) = Neo.LZC(irow);
    Alltable.CTW(rowcnt) = Neo.CTW(irow);
    Alltable.Stimulus(rowcnt) = Neo.Stimulus(irow);
    Alltable.Violation(rowcnt) = Neo.Violation(irow);
    Alltable.GlobalRule(rowcnt) = xor(Neo.Stimulus(irow),Neo.Violation(irow));
    Alltable.BornYet(rowcnt) = 1;
    Alltable.GA(rowcnt) = floor((Neo.age(irow) + GAAB)/7); % postmenstral age in weeks
    Alltable.BirthAge(rowcnt) = GAAB;
    Alltable.MomsAge(rowcnt) = Tdemo.Var8(Tidx);
    Alltable.BirthWeight(rowcnt) = Tdemo.Var11(Tidx);
end

%% Scatter plot of average entropy at each age for FETUSES

% This code generates scatter plots of various measures (PermEn32, PermEn64, LZC, CTW, mSampEn, and mMSE) for male and female fetuses at different gestational ages. Here's a summary:
% 
%     Data Preparation:
%         The code starts by defining the age range (AGES) from 25 to 40 weeks.
%         Separate arrays are created to store the mean values for males (maleavg) and females (femaleavg) for each measure and gestational age.
% 
%     Data Separation:
%         The code separates the data for male (boys) and female (girls) samples, as well as for fetuses (fetuses) and newborns (newborns).
% 
%     Loop through Gestational Ages:
%         The code uses a loop to iterate over each gestational age in the specified range.
%         For each age, it calculates the median values of various measures (PE32, PE64, SE, MSE, LZC, CTW) for male and female fetuses separately.
% 
%     Combine Male and Female Averages:
%         The computed average values for males and females are then combined into arrays (avgPE32, avgPE64, etc.).
% 
%     Scatter Plots:
%         For each measure, the code creates scatter plots for male and female fetuses at different gestational ages.
%         It also calculates the correlation coefficient (r) and p-value for the relationship between gestational age and the corresponding measure.
% 
%     Connect data from the same individuals
%         The code connects the data points for each individual
% 
%     Plot Adjustments:
%         The code adjusts the axes, labels, titles, and legends for each plot.
% 
%     Saving Plots:
%         Finally, the code saves the generated scatter plots as PNG and SVG files.

AGES = 25:40;
maleavg = nan(2,length(AGES));
femaleavg = nan(2,length(AGES));

boys = strcmp(Alltable.Sex,'M');
girls = strcmp(Alltable.Sex,'F');

fetuses = Alltable.BornYet==0;
newborns = Alltable.BornYet==1;

cnt = 0;
for iage = 1:length(AGES)
    cnt = cnt + 1;
    midx = Alltable.GA == AGES(iage) & strcmp(Alltable.Sex,'M');
    maleavgPE32(iage) = median(Alltable.PE32(midx & fetuses));
    maleavgPE64(iage) = median(Alltable.PE64(midx & fetuses));
    maleavgSE(iage)   = median(Alltable.SE(midx & fetuses));
    maleavgMSE(iage)  = median(Alltable.MSE(midx & fetuses));
    maleavgLZC(iage)  = median(Alltable.LZC(midx & fetuses));
    maleavgCTW(iage)  = median(Alltable.CTW(midx & fetuses));

    fidx = Alltable.GA == AGES(iage) & strcmp(Alltable.Sex,'F');
    femaleavgPE32(iage) = median(Alltable.PE32(fidx & fetuses));
    femaleavgPE64(iage) = median(Alltable.PE64(fidx & fetuses));
    femaleavgSE(iage)   = median(Alltable.SE(fidx & fetuses));
    femaleavgMSE(iage)  = median(Alltable.MSE(fidx & fetuses));
    femaleavgLZC(iage)  = median(Alltable.LZC(fidx & fetuses));
    femaleavgCTW(iage)  = median(Alltable.CTW(fidx & fetuses));
end

avgPE32 = [maleavgPE32; femaleavgPE32];
avgPE64 = [maleavgPE64; femaleavgPE64];
avgSE = [maleavgSE; femaleavgSE];
avgMSE = [maleavgMSE; femaleavgMSE];
avgLZC = [maleavgLZC; femaleavgLZC];
avgCTW = [maleavgCTW; femaleavgCTW];

% PermEn32
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgPE32,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgPE32','rows','complete');
mylsline
scatter(Alltable.GA(boys&fetuses),Alltable.PE32(boys&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.PE32(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE32(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE32(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE32(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn32')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE32;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgPE32,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgPE32','rows','complete');
mylsline
scatter(Alltable.GA(girls&fetuses),Alltable.PE32(girls&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.PE32(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE32(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE32(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE32(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn32')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE32;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/PE32_scatter')
print('-dsvg','./Figures/2023/PE32_scatter')

% PermEn64
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgPE64,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgPE64','rows','complete');
mylsline
scatter(Alltable.GA(boys&fetuses),Alltable.PE64(boys&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.PE64(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE64(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE64(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE64(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn64')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE64;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgPE64,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgPE64','rows','complete');
mylsline
scatter(Alltable.GA(girls&fetuses),Alltable.PE64(girls&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.PE64(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE64(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE64(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE64(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn64')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE64;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/PE64_scatter')
print('-dsvg','./Figures/2023/PE64_scatter')

% LZC
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgLZC,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgLZC','rows','complete');
mylsline
scatter(Alltable.GA(boys&fetuses),Alltable.LZC(boys&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.LZC(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.LZC(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.LZC(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.LZC(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average LZC')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgLZC;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgLZC,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgLZC','rows','complete');
mylsline
scatter(Alltable.GA(girls&fetuses),Alltable.LZC(girls&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.LZC(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.LZC(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.LZC(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.LZC(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average LZC')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgLZC;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/LZC_scatter')
print('-dsvg','./Figures/2023/LZC_scatter')

% CTW
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgCTW,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgCTW','rows','complete');
mylsline
scatter(Alltable.GA(boys&fetuses),Alltable.CTW(boys&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.CTW(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.CTW(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.CTW(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.CTW(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average CTW')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgCTW;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgCTW,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgCTW','rows','complete');
mylsline
scatter(Alltable.GA(girls&fetuses),Alltable.CTW(girls&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.CTW(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.CTW(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.CTW(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.CTW(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average CTW')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgCTW;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/CTW_scatter')
print('-dsvg','./Figures/2023/CTW_scatter')

% mSampEn
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgSE','rows','complete');
mylsline
scatter(Alltable.GA(boys&fetuses),Alltable.SE(boys&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.SE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.SE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.SE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.SE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mSampEn')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgSE;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgSE','rows','complete');
mylsline
scatter(Alltable.GA(girls&fetuses),Alltable.SE(girls&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.SE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.SE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.SE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.SE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mSampEn')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgSE;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/SE_scatter')
print('-dsvg','./Figures/2023/SE_scatter')

% mMSE
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgMSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgMSE','rows','complete');
mylsline
scatter(Alltable.GA(boys&fetuses),Alltable.MSE(boys&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.MSE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.MSE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.MSE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.MSE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mMSE')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgMSE;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgMSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgMSE','rows','complete');
mylsline
scatter(Alltable.GA(girls&fetuses),Alltable.MSE(girls&fetuses),20,'+','markeredgecolor','r')

% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & ~Alltable.BornYet);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & ~Alltable.BornYet);
    plot(Alltable.GA(sidx1),Alltable.MSE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.MSE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.MSE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.MSE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mMSE')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgMSE;
axis( [24 40 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/MSE_scatter')
print('-dsvg','./Figures/2023/MSE_scatter')


%% Scatter plot of average entropy at each age for NEWBORNS
%
% This is just the same as the above block of code, but it generates
% scatter plots for neonatal data rather than fetal data.
%
AGES = 40:50;
maleavg = nan(2,length(AGES));
femaleavg = nan(2,length(AGES));

boys = strcmp(Alltable.Sex,'M');
girls = strcmp(Alltable.Sex,'F');

maleavgPE32 = nan(1,length(AGES)); 
maleavgPE64 = nan(1,length(AGES)); 
maleavgSE = nan(1,length(AGES)); 
maleavgMSE = nan(1,length(AGES)); 
maleavgLZC = nan(1,length(AGES)); 
maleavgCTW = nan(1,length(AGES)); 

femaleavgPE32 = nan(1,length(AGES));
femaleavgPE64 = nan(1,length(AGES)); 
femaleavgSE = nan(1,length(AGES));
femaleavgMSE = nan(1,length(AGES)); 
femaleavgLZC = nan(1,length(AGES));
femaleavgCTW = nan(1,length(AGES)); 

cnt = 0;
for iage = 1:length(AGES)
    cnt = cnt + 1;
    midx = Alltable.GA == AGES(iage) & strcmp(Alltable.Sex,'M');
    maleavgPE32(iage) = median(Alltable.PE32(midx & newborns));
    maleavgPE64(iage) = median(Alltable.PE64(midx & newborns));
    maleavgSE(iage)   = median(Alltable.SE(midx & newborns));
    maleavgMSE(iage)  = median(Alltable.MSE(midx & newborns));
    maleavgLZC(iage)  = median(Alltable.LZC(midx & newborns));
    maleavgCTW(iage)  = median(Alltable.CTW(midx & newborns));

    fidx = Alltable.GA == AGES(iage) & strcmp(Alltable.Sex,'F');
    femaleavgPE32(iage) = median(Alltable.PE32(fidx & newborns));
    femaleavgPE64(iage) = median(Alltable.PE64(fidx & newborns));
    femaleavgSE(iage)   = median(Alltable.SE(fidx & newborns));
    femaleavgMSE(iage)  = median(Alltable.MSE(fidx & newborns));
    femaleavgLZC(iage)  = median(Alltable.LZC(fidx & newborns));
    femaleavgCTW(iage)  = median(Alltable.CTW(fidx & newborns));
end


avgPE32 = [maleavgPE32; femaleavgPE32];
avgPE64 = [maleavgPE64; femaleavgPE64];
avgSE = [maleavgSE; femaleavgSE];
avgMSE = [maleavgMSE; femaleavgMSE];
avgLZC = [maleavgLZC; femaleavgLZC];
avgCTW = [maleavgCTW; femaleavgCTW];

% PermEn32
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgPE32,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgPE32','rows','complete');
mylsline
scatter(Alltable.GA(boys&newborns),Alltable.PE32(boys&newborns),20,'+','markeredgecolor','b')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.PE32(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE32(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE32(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE32(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn32')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE32;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgPE32,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgPE32','rows','complete');
mylsline

scatter(Alltable.GA(girls&newborns),Alltable.PE32(girls&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.PE32(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE32(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE32(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE32(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn32')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE32;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/PE32_scatter_newborn')
print('-dsvg','./Figures/2023/PE32_scatter_newborn')

% PermEn64
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgPE64,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgPE64','rows','complete');
mylsline
scatter(Alltable.GA(boys&newborns),Alltable.PE64(boys&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.PE64(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE64(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE64(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE64(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn64')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE64;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgPE64,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgPE64','rows','complete');
mylsline
scatter(Alltable.GA(girls&newborns),Alltable.PE64(girls&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.PE64(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.PE64(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.PE64(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.PE64(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average PermEn64')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgPE64;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/PE64_scatter_newborn')
print('-dsvg','./Figures/2023/PE64_scatter_newborn')

% LZC
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgLZC,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgLZC','rows','complete');
mylsline
scatter(Alltable.GA(boys&newborns),Alltable.LZC(boys&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.LZC(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.LZC(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.LZC(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.LZC(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average LZC')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgLZC;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgLZC,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgLZC','rows','complete');
mylsline
scatter(Alltable.GA(girls&newborns),Alltable.LZC(girls&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.LZC(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.LZC(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.LZC(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.LZC(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average LZC')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgLZC;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/LZC_scatter_newborn')
print('-dsvg','./Figures/2023/LZC_scatter_newborn')

% CTW

myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgCTW,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgCTW','rows','complete');
mylsline
scatter(Alltable.GA(boys&newborns),Alltable.CTW(boys&newborns),20,'+','markeredgecolor','b')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.CTW(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.CTW(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.CTW(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.CTW(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average CTW')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgCTW;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgCTW,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgCTW','rows','complete');
mylsline
scatter(Alltable.GA(girls&newborns),Alltable.CTW(girls&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.CTW(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.CTW(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.CTW(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.CTW(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average CTW')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgCTW;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/CTW_scatter_newborn')
print('-dsvg','./Figures/2023/CTW_scatter_newborn')

% mSampEn
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgSE','rows','complete');
mylsline
scatter(Alltable.GA(boys&newborns),Alltable.SE(boys&newborns),20,'+','markeredgecolor','b')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.SE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.SE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.SE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.SE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mSampEn')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgSE;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgSE','rows','complete');
mylsline
scatter(Alltable.GA(girls&newborns),Alltable.SE(girls&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.SE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.SE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.SE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.SE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mSampEn')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgSE;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/SE_scatter_newborn')
print('-dsvg','./Figures/2023/SE_scatter_newborn')

% mMSE
myfigure2
subplot(1,2,1)
hold on
scatter(AGES,maleavgMSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',maleavgMSE','rows','complete');
mylsline
scatter(Alltable.GA(boys&newborns),Alltable.MSE(boys&newborns),20,'+','markeredgecolor','b')

% connect the dots
subs = unique(Alltable.ID(boys));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.MSE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.MSE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.MSE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.MSE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mMSE')
title(sprintf('Male, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgMSE;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty

subplot(1,2,2)
hold on
scatter(AGES,femaleavgMSE,100,'fill','markerfacecolor','m')
[r,p] = corr(AGES',femaleavgMSE','rows','complete');
mylsline
scatter(Alltable.GA(girls&newborns),Alltable.MSE(girls&newborns),20,'+','markeredgecolor','b')


% connect the dots
subs = unique(Alltable.ID(girls));
for isub = 1:length(subs)
    sidx1 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx2 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & Alltable.Stimulus & newborns);
    sidx3 = find(Alltable.ID == subs(isub) & Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    sidx4 = find(Alltable.ID == subs(isub) & ~Alltable.GlobalRule & ~Alltable.Stimulus & newborns);
    plot(Alltable.GA(sidx1),Alltable.MSE(sidx1),'k:')
    if isub==2
        legend({'Average entropy per age bin','Least squares fit to averages','Fetal',...
            'Longitudinal'},...
            'AutoUpdate','Off','fontsize',10,'location','southwest')
        legend box off
    end
    plot(Alltable.GA(sidx2),Alltable.MSE(sidx2),'k:')
    plot(Alltable.GA(sidx3),Alltable.MSE(sidx3),'k:')
    plot(Alltable.GA(sidx4),Alltable.MSE(sidx4),'k:')
end
xlabel('Gestational age (weeks)')
ylabel('Average mMSE')
title(sprintf('Female, R^2 = %1.2f',r^2),'fontsize',15)
Hm = avgMSE;
axis( [38 50 min(Hm(:))-nanstd(Hm(:)) max(Hm(:))+nanstd(Hm(:)) ])
makefigpretty
print('-dpng','./Figures/2023/MSE_scatter_newborn')
print('-dsvg','./Figures/2023/MSE_scatter_newborn')


%% Compute ICCs

% This code block performs ICC analysis for different measures, compares ICC values between fetal and combined fetal-newborn datasets, and outputs the results in CSV and LaTeX formats. The results include ICC values, p-values, and significance levels for each measure.
%     Initialization:
%         subs: Unique IDs in the dataset.
%         Ns: Number of subjects.
%         nperm: Number of permutations for the significanc test.
% 
%     Data Preparation:
%         Matrices (Mlzc, Mctw, Mpe32, Mpe64, Mmse, Mse) are initialized to store the measures for each subject.
%         For each subject, the corresponding measures (LZC, CTW, PE32, PE64, MSE, SE) are added to the matrices.
% 
%     ICC Calculation (fetuses + newborns)
%         NOTE: the calculation for the combined data was NOT included in the manuscript! 
%         ICC is calculated for each measure using the myICC function.
%         Parallel computing is used for permutation testing to obtain a distribution of ICC values under the null hypothesis.
%         A permutation loop (parfor) shuffles the data and calculates ICC for each permutation.
% 
%     ICC Calculation (just fetuses):
%         Similar ICC calculations are performed using only fetal data.
%         These are the ICCs that we do report in the manuscript table
% 
%     Permutation Test:
%         Another permutation loop (parfor) is used for the second set of ICC calculations.
% 
%     P-Value Calculation:
%         P-values are calculated based on the observed ICC values compared to the null distribution from permutation testing.
% 
%     Results Output:
%         ICC values and p-values for both fetal and combined fetal-newborn datasets are printed.
%         A table (Ticc) is created with the ICC values and p-values.
%         The table is written to a CSV file (ICCs.csv) and a LaTeX file (ICCs_LaTeX).
% 

subs = unique(Alltable.ID);
Ns = length(subs);
nperm = 1000;

Mlzc = nan(20,Ns);
Mctw = nan(20,Ns);
Mpe32 = nan(20,Ns);
Mpe64 = nan(20,Ns);
Mmse = nan(20,Ns);
Mse = nan(20,Ns);


% ICC (fetuses + newborns, NOT included in the manuscript)
for isub = 1:Ns
    Tsub = Alltable(Alltable.ID == subs(isub),:);

    Mlzc(1:size(Tsub,1),isub) = Tsub.LZC;
    Mctw(1:size(Tsub,1),isub) = Tsub.CTW;
    Mpe32(1:size(Tsub,1),isub) = Tsub.PE32;
    Mpe64(1:size(Tsub,1),isub) = Tsub.PE64;
    Mmse(1:size(Tsub,1),isub) = Tsub.MSE;
    Mse(1:size(Tsub,1),isub) = Tsub.SE;
end

r_LZC = myICC(Mlzc', '1-1', 0.05);
r_CTW = myICC(Mctw', '1-1', 0.05);
r_PE32 = myICC(Mpe32', '1-1', 0.05);
r_PE64 = myICC(Mpe64', '1-1', 0.05);
r_MSE = myICC(Mmse', '1-1', 0.05);
r_SE = myICC(Mse', '1-1', 0.05);

% Permutation test
try
    parpool
catch
    fprintf('Looks like parpool is already open\n')
end

rp_LZC = nan(1,nperm);
rp_CTW = nan(1,nperm);
rp_PE32 = nan(1,nperm);
rp_PE64 = nan(1,nperm);
rp_MSE = nan(1,nperm);
rp_SE = nan(1,nperm);

parfor iperm = 1:nperm
    tmp = Mlzc(:);
    tmp2 = tmp(randperm(length(tmp)));
    Plzc = reshape(tmp2,20,Ns);
    rp_LZC(iperm) = myICC(Plzc', '1-1', 0.05);

    tmp = Mctw(:);
    tmp2 = tmp(randperm(length(tmp)));
    Pctw = reshape(tmp2,20,Ns);
    rp_CTW(iperm) = myICC(Pctw', '1-1', 0.05);

    tmp = Mpe32(:);
    tmp2 = tmp(randperm(length(tmp)));
    Ppe32 = reshape(tmp2,20,Ns);
    rp_PE32(iperm) = myICC(Ppe32', '1-1', 0.05);

    tmp = Mpe64(:);
    tmp2 = tmp(randperm(length(tmp)));
    Ppe64 = reshape(tmp2,20,Ns);
    rp_PE64(iperm) = myICC(Ppe64', '1-1', 0.05);

    tmp = Mmse(:);
    tmp2 = tmp(randperm(length(tmp)));
    Pmse = reshape(tmp2,20,Ns);
    rp_MSE(iperm) = myICC(Pmse', '1-1', 0.05);

    tmp = Mse(:);
    tmp2 = tmp(randperm(length(tmp)));
    Pse = reshape(tmp2,20,Ns);
    rp_SE(iperm) = myICC(Pse', '1-1', 0.05);
    fprintf('%1.1f%% finished \n', iperm/nperm*100)
end


Mlzc = nan(20,Ns);
Mctw = nan(20,Ns);
Mpe32 = nan(20,Ns);
Mpe64 = nan(20,Ns);
Mmse = nan(20,Ns);
Mse = nan(20,Ns);


% ICC (just fetuses, this is the one that we used in the manuscript)
for isub = 1:Ns
    Tsub = Alltable(Alltable.ID == subs(isub) & ~Alltable.BornYet,:);

    Mlzc(1:size(Tsub,1),isub) = Tsub.LZC;
    Mctw(1:size(Tsub,1),isub) = Tsub.CTW;
    Mpe32(1:size(Tsub,1),isub) = Tsub.PE32;
    Mpe64(1:size(Tsub,1),isub) = Tsub.PE64;
    Mmse(1:size(Tsub,1),isub) = Tsub.MSE;
    Mse(1:size(Tsub,1),isub) = Tsub.SE;
end

rf_LZC = myICC(Mlzc', '1-1', 0.05);
rf_CTW = myICC(Mctw', '1-1', 0.05);
rf_PE32 = myICC(Mpe32', '1-1', 0.05);
rf_PE64 = myICC(Mpe64', '1-1', 0.05);
rf_MSE = myICC(Mmse', '1-1', 0.05);
rf_SE = myICC(Mse', '1-1', 0.05);

rfp_LZC = nan(1,nperm);
rfp_CTW = nan(1,nperm);
rfp_PE32 = nan(1,nperm);
rfp_PE64 = nan(1,nperm);
rfp_MSE = nan(1,nperm);
rfp_SE = nan(1,nperm);

% another permutation test
parfor iperm = 1:nperm
    tmp = Mlzc(:);
    tmp2 = tmp(randperm(length(tmp)));
    Plzc = reshape(tmp2,20,Ns);
    rfp_LZC(iperm) = myICC(Plzc', '1-1', 0.05);

    tmp = Mctw(:);
    tmp2 = tmp(randperm(length(tmp)));
    Pctw = reshape(tmp2,20,Ns);
    rfp_CTW(iperm) = myICC(Pctw', '1-1', 0.05);

    tmp = Mpe32(:);
    tmp2 = tmp(randperm(length(tmp)));
    Ppe32 = reshape(tmp2,20,Ns);
    rfp_PE32(iperm) = myICC(Ppe32', '1-1', 0.05);

    tmp = Mpe64(:);
    tmp2 = tmp(randperm(length(tmp)));
    Ppe64 = reshape(tmp2,20,Ns);
    rfp_PE64(iperm) = myICC(Ppe64', '1-1', 0.05);

    tmp = Mmse(:);
    tmp2 = tmp(randperm(length(tmp)));
    Pmse = reshape(tmp2,20,Ns);
    rfp_MSE(iperm) = myICC(Pmse', '1-1', 0.05);

    tmp = Mse(:);
    tmp2 = tmp(randperm(length(tmp)));
    Pse = reshape(tmp2,20,Ns);
    rfp_SE(iperm) = myICC(Pse', '1-1', 0.05);
    fprintf('%1.1f%% finished \n', iperm/nperm*100)
end


Pval_lzc = (nperm-sum(r_LZC > rp_LZC))/nperm;
Pval_ctw = (nperm-sum(r_CTW > rp_CTW))/nperm;
Pval_pe32 = (nperm-sum(r_PE32 > rp_PE32))/nperm;
Pval_pe64 = (nperm-sum(r_PE64 > rp_PE64))/nperm;
Pval_mse = (nperm-sum(r_MSE > rp_MSE))/nperm;
Pval_se = (nperm-sum(r_SE > rp_SE))/nperm;

fPval_lzc = (nperm-sum(rf_LZC > rfp_LZC))/nperm;
fPval_ctw = (nperm-sum(rf_CTW > rfp_CTW))/nperm;
fPval_pe32 = (nperm-sum(rf_PE32 > rfp_PE32))/nperm;
fPval_pe64 = (nperm-sum(rf_PE64 > rfp_PE64))/nperm;
fPval_mse = (nperm-sum(rf_MSE > rfp_MSE))/nperm;
fPval_se = (nperm-sum(rf_SE > rfp_SE))/nperm;


fprintf('ICC for LZC = %1.3f (fetal, P = %1.3f), %1.3f (all data, P = %1.3f)\n',r_LZC,Pval_lzc,rf_LZC,fPval_lzc)
fprintf('ICC for CTW = %1.3f (fetal, P = %1.3f), %1.3f (all data, P = %1.3f)\n',r_CTW,Pval_ctw,rf_CTW,fPval_ctw)
fprintf('ICC for PE32 = %1.3f (fetal, P = %1.3f), %1.3f (all data, P = %1.3f)\n',r_PE32,Pval_pe32,rf_PE32,fPval_pe32)
fprintf('ICC for PE64 = %1.3f (fetal, P = %1.3f), %1.3f (all data, P = %1.3f)\n',r_PE64,Pval_pe64,rf_PE64,fPval_pe64)
fprintf('ICC for MSE = %1.3f (fetal, P = %1.3f), %1.3f (all data, P = %1.3f)\n',r_MSE,Pval_mse,rf_MSE,fPval_mse)
fprintf('ICC for SE = %1.3f (fetal, P = %1.3f), %1.3f (all data, P = %1.3f)\n',r_SE,Pval_se,rf_SE,fPval_se)

Ticc = table(measures',[r_CTW; r_LZC; r_SE; r_MSE; r_PE32; r_PE64],...
    [Pval_ctw; Pval_lzc; Pval_se; Pval_mse; Pval_pe32; Pval_pe64],...
    [rf_CTW; rf_LZC; rf_SE; rf_MSE; rf_PE32; rf_PE64], ...
    [fPval_ctw; fPval_lzc; fPval_se; fPval_mse; fPval_pe32; fPval_pe64],'VariableNames',...
    {'Measure','fetal_ICC','fetal_Pvalues','all_ICC','all_pvalues'});

Ticc = Ticc([1 2 5 6 4 3],:) % reorder table

writetable(Ticc,'./ICCs.csv')
Ticc{:,2} = round(Ticc{:,2},3,'significant');
Ticc{:,3} = round(Ticc{:,3},3,'significant');
table2latex(Ticc,'./ICCs_LaTeX')



%% Run new models, this time with all data pooled together (newborns and fetuses)

% This code block fits linear mixed-effects models for each measure, consolidates the significant predictors into a common formula, and extracts relevant statistics for GA, Sex, and BornYet. The results are saved in CSV and LaTeX formats.

% Data Preparation:
% 
%     Columns (GlobalRule, Stimulus, BornYet) are converted to logical type.
% 
% Model Fitting Loop:
% 
%     For each measure, a stepwise linear mixed-effects model (stepwiselm) is fitted without random effects.
%     The formula for the model includes various predictors such as gestational age (GA), standard deviation of NN intervals (SDNN), sex (Sex), birth status (BornYet), global rule (GlobalRule), stimulus (Stimulus), birth age (BirthAge), mother's age (MomsAge), and birth weight (BirthWeight).
% 
% Formula Consolidation:
% 
%     Unique predictors across all models are identified and consolidated into a common formula (mdlstrall).
%     The formula includes predictors present in at least half of the models for the entropy features.
% 
% Linear Mixed-Effects Model Fitting:
% 
%     The final model is fitted using the consolidated formula, including a random effect for subject ID ((1|ID)).
% 
% Results Extraction:
% 
%     For each measure, p-values and t-statistics are extracted for predictors GA, Sex, and BornYet.
% 
% Table Creation and Output:
% 
%     A table (Tout) is created with columns for measure, t-statistics, and p-values for GA, BornYet, and Sex.
%     The table is written to a CSV file (AllModel.csv) and a LaTeX file (AllModel).

Alltable.GlobalRule = logical(Alltable.GlobalRule);
Alltable.Stimulus = logical(Alltable.Stimulus);
Alltable.BornYet = logical(Alltable.BornYet);

for i = 1:length(measures)
    % Stepwise model without random effects
    fprintf('Formula: %s\n', sprintf('%s ~ 1+BornYet',measures{i}));
    mdl = stepwiselm(Alltable,sprintf('%s ~ 1+BornYet',measures{i}),'ResponseVar',...
        measures{i},'PredictorVars',{'GA','SDNN','Sex','BornYet','GlobalRule',...
        'Stimulus','BirthAge','MomsAge','BirthWeight'},...
        'CategoricalVar',{'Sex','GlobalRule','Stimulus','BornYet'},'Verbose',2);
    %sprintf('Using this formula: %s\n',mdl.Formula)
    formula1 = char(mdl.Formula);
    start = strfind(formula1,'~');
    mmall{i,1} = char(formula1(start+1:end));
end

allf = strcat(mmall{1},mmall{2},mmall{3},mmall{4},mmall{5},mmall{6});
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

mdlstrall = '';
cnt = 0;
for itok = 1:length(tokens)
    if Ntok(itok) >= length(measures)/2 % if this predictor is present for at least half the entropy features
        cnt = cnt + 1;
        if cnt == 1
            mdlstrall = sprintf('%s ~ %s',mdlstrall,tokens{itok});
        else
            mdlstrall = sprintf('%s + %s',mdlstrall,tokens{itok});
        end
    end
end
mdlstrall = sprintf('%s + (1|ID)',mdlstrall)


PvalsGA = nan(6,1);
PvalsSex = nan(6,1);
PvalsBornYet = nan(6,1);
tStatGA = nan(6,1);
tStatSex = nan(6,1);
tStatBornYet = nan(6,1);


% run the modal model
for i = 1:length(measures)
    lme_formula = sprintf('%s %s',measures{i},mdlstrall);
    out = fitlme(Alltable,lme_formula);
    PvalsGA(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'GA'));
    PvalsSex(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'Sex_F'));
    PvalsBornYet(i) = out.Coefficients.pValue(strcmp(out.Coefficients.Name,'BornYet_1'));
    tStatGA(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'GA'));
    tStatSex(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'Sex_F'));
    tStatBornYet(i) = out.Coefficients.tStat(strcmp(out.Coefficients.Name,'BornYet_1'));
end


Tout = table(measures',tStatGA,PvalsGA,...
    tStatBornYet,PvalsBornYet,tStatSex,PvalsSex,...
    'variablenames',{'Measure','GA_tstat','GA_Pvalue',...
    'BornYet_tstat','BornYet_Pvalue',...
    'Sex_tstat','Sex_Pvalue'})

writetable(Tout, 'AllModel.csv')
table2latex(Tout, 'AllModel')

clear out Tout % Clears unnecessary variables from the workspace


%% Test for difference between fetuses and newborns in proportion of
% stochastic vs determinstic dynamics

% This code block tests whether the proportion of stochastic dynamics significantly differs between fetuses and newborns and reports the results of the hypothesis test.
% 
%     Proportion Calculation:
%         The code calculates the proportions of stochastic dynamics for newborns (Neo) and fetuses (P).
% 
%     Chi^2 Test:
%         A Chi squaredtest is conducted using the prop_test function, comparing the proportions of stochastic dynamics in newborns and fetuses.
%         The null hypothesis (h) is that there is no significant difference between the proportions.
% 
%     Result Display:
%         Depending on the outcome of the test (h), the code prints a message indicating whether the proportions significantly differ or not.
%         The p-value (p) and chi-squared statistic (chi2stat) are also displayed.


[h,p, chi2stat,df] = prop_test([sum(Neo.IsSto) sum(P.Sto)] , [size(Neo,1) size(P,1)],0);

switch h
    case true
        fprintf('No subtraction: Proportions significantly differ, p = %1.3e, chi^2 = %1.3f\n',p,chi2stat)
    case false
        fprintf('No subtractioN: Proportions do not significantly differ, p = %1.3e, chi^2 = %1.3f\n',p,chi2stat)
end



%% Time-frequency results %%

% This code performs threshold-free cluster enhancement (TFCE) on the results of linear mixed-effects models applied to MEG data.

%     Setup:
%         Set random number generator seed (rng(111)).
%         Define parameters such as the number of permutations (Nperm), minimum and maximum values for lasso regression (lambda_min and lambda_max), and whether to allow interactions (interact).
%         Set up frequency transforms using the ro_freq_wavelet_TFT function.
% 
%     Data Preprocessing:
%         Apply frequency transforms to obtain time-frequency representations (powDevA, powStdA, powDevB, powStdB).
%         Ensure the time-frequency representations have matching sizes.
% 
%     Create a Table for Predictors:
%         Create a table (Tm) with predictors such as gestational age (GA), standard deviation of NN intervals (SDNN), sex, rule, stimulus, BMI, birth age, and birth weight.
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

rng(111)


Nperm = 200; % how many permutations to run
lambda_min = 0.02; % most permissive lasso regression to try
lambda_max = 2; % least permissive lasso regression to try
interact = false; % allow interactions? 

% frequency transforms
% These settings yield a 44 x 28 TFT, which we will downsample to 11 x 7
% for the pupose of model fitting before upscaling again
cfg = [];
cfg.fsample = fs;
cfg.foi_start = 1;
cfg.foi_end = 10.5; % make sure 10 Hz is included
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

Tm = table(Bigtable.ID,Bigtable.GA,Bigtable.SDNN,Bigtable.Stimulus,...
    Bigtable.GlobalRule,Bigtable.XChrom,Bigtable.MomsAge,...
    Bigtable.BMI,Bigtable.BirthAge,Bigtable.BirthWeight,...
    'VariableNames',{'ID','GA','SDNN','Stimulus','Rule','Sex','BMI','MomsAge','BirthAge','BirthWeight'});

kfold = 10;
POW = [powStdA; powStdB; powDevA; powDevB];
P = POW(:);
Tm2 = repmat(Tm,size(POW,2)*size(POW,3),1);
predictors = [Tm2.GA Tm2.SDNN Tm2.Sex Tm2.Rule Tm2.Stimulus Tm2.BMI Tm2.BirthAge Tm2.BirthWeight Tm2.MomsAge];
terms = {'GA','SDNN','Sex','Rule','Stimulus','BMI','BirthAge','BirthWeight','MomsAge'};

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

L = logspace(log10(lambda_min),log10(lambda_max),25);
[B0,FitInfo] =lassoglm(predictors, P,'normal','CV',kfold,'lambda',L,'Alpha',1);
[~,useme] = min(FitInfo.Deviance); % find optimal regularization
lambda = FitInfo.Lambda(useme); % regularization parameter
[B, FitInfo] = lassoglm(predictors, P, 'normal', 'Alpha', 1, 'Lambda', lambda);

IVs = newterms(find(abs(B)>0)) % select independent variables for our model

formula = '1';
for ivar = 1:length(IVs)
    formula = sprintf('%s + %s',formula,IVs{ivar});
end

save fetal_meg_stats_save_progress

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
%     The code checks for specific terms in the fitted model (e.g., 'MomsAge', 'BirthAge', 'GA', etc.) using ismember.
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
%     The results are saved in variables like GA_statsum_pos, AgeMA_statsum_pos, etc., which contain the TFCE results for positive values, and similarly for negative values.
%     The entire process's execution time is measured using tic and toc.
% 
% Saving to File:
% 
%     The final results are saved to a file named based on the current date (fetal_meg_stats_progress_%s).

PvalsBirthAge = nan(size(POW,2),size(POW,3));
PvalsAgeMA = nan(size(POW,2),size(POW,3));
PvalsMABA = nan(size(POW,2),size(POW,3));
PvalsMomsAge = nan(size(POW,2),size(POW,3));
PvalsGA = nan(size(POW,2),size(POW,3));
PvalsSex = nan(size(POW,2),size(POW,3));
PvalsSDNN = nan(size(POW,2),size(POW,3));
PvalsRule = nan(size(POW,2),size(POW,3));
PvalsStimulus = nan(size(POW,2),size(POW,3));
PvalsBMI = nan(size(POW,2),size(POW,3));

tStatBirthAge = nan(size(POW,2),size(POW,3));
tStatAgeMA = nan(size(POW,2),size(POW,3));
tStatMABA = nan(size(POW,2),size(POW,3));
tStatMomsAge = nan(size(POW,2),size(POW,3));
tStatGA = nan(size(POW,2),size(POW,3));
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


        if any(ismember(lmeout.Coefficients.Name,'MomsAge'))
            PvalsMomsAge(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'MomsAge'));
            tStatMomsAge(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'MomsAge'));
        end

        if any(ismember(lmeout.Coefficients.Name,'BirthAge'))
            PvalsBirthAge(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'BirthAge'));
            tStatBirthAge(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'BirthAge'));
        end

        if any(ismember(lmeout.Coefficients.Name,'GA'))
            PvalsGA(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'GA'));
            tStatGA(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'GA'));
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

[~,~,GA_statsum_pos] = TFCE(PvalsGA,tStatGA,'positive',0);
[~,~,AgeMA_statsum_pos] = TFCE(PvalsAgeMA,tStatAgeMA,'positive',0);
[~,~,MABA_statsum_pos] = TFCE(PvalsMABA,tStatMABA,'positive',0);
[~,~,MomsAge_statsum_pos] = TFCE(PvalsMomsAge,tStatMomsAge,'positive',0);
[~,~,BirthAge_statsum_pos] = TFCE(PvalsBirthAge,tStatBirthAge,'positive',0);
[~,~,Sex_statsum_pos] = TFCE(PvalsSex,tStatSex,'positive',0);
[~,~,SDNN_statsum_pos] = TFCE(PvalsSDNN,tStatSDNN,'positive',0);
[~,~,Rule_statsum_pos] = TFCE(PvalsRule,tStatRule,'positive',0);
[~,~,Stimulus_statsum_pos] = TFCE(PvalsStimulus,tStatStimulus,'positive',0);
[~,~,BMI_statsum_pos] = TFCE(PvalsBMI,tStatBMI,'positive',0);


[~,~,GA_statsum_neg] = TFCE(PvalsGA,tStatGA,'negative',0);
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
saveme = sprintf('fetal_meg_stats_progress_%s',date);
save(saveme)

%% now run permutations

% This code block performs permutations for threshold-free cluster enhancement (TFCE) for statistical analysis

% Initialization:
% 
%     Matrices (GA_ssPerm_neg, AgeMA_ssPerm_neg, ..., BMI_ssPerm_pos) are initialized to store the TFCE results for negative and positive values for different variables and model parameters. These matrices represent the results of permutation tests.
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
%     P-values and t-statistics for various coefficients are stored in matrices (PvalsGAPerm, PvalsAgeMAPerm, ..., tStatBMIPerm).
% 
% TFCE Application:
% 
%     TFCE is applied separately for negative and positive values for each variable and permutation. This involves enhancing clusters of significant values in the permutation results.
% 
% Results Storage:
% 
%     The TFCE results for negative and positive values are stored in the corresponding matrices (GA_ssPerm_neg, AgeMA_ssPerm_neg, ..., BMI_ssPerm_pos).
% 
% Progress Tracking:
% 
%     The code prints the completion percentage of the permutation loop to track the progress of the analysis.
% 
% Results Saving:
% 
%     After completing all permutations, the results are saved to a file with a name based on the current date (fetal_perm_stats_TFCEs_%s).


tic

GA_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
AgeMA_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
MABA_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
MomsAge_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
BirthAge_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
Sex_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
SDNN_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
Rule_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
Stimulus_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));
BMI_ssPerm_neg = nan(Nperm,size(powDevA,2),size(powDevA,3));


GA_ssPerm_pos = nan(Nperm,size(powDevA,2),size(powDevA,3));
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
    PvalsMABAPerm = nan(size(POW,2),size(POW,3));
    PvalsAgeMAPerm = nan(size(POW,2),size(POW,3));
    PvalsBirthAgePerm = nan(size(POW,2),size(POW,3));
    PvalsMomsAgePerm  = nan(size(POW,2),size(POW,3));
    PvalsGAPerm       = nan(size(POW,2),size(POW,3));
    PvalsSexPerm = nan(size(POW,2),size(POW,3));
    PvalsSDNNPerm = nan(size(POW,2),size(POW,3));
    PvalsRulePerm  = nan(size(POW,2),size(POW,3));
    PvalsStimulusPerm       = nan(size(POW,2),size(POW,3));
    PvalsBMIPerm       = nan(size(POW,2),size(POW,3));

    tStatMABAPerm = nan(size(POW,2),size(POW,3));
    tStatAgeMAPerm = nan(size(POW,2),size(POW,3));
    tStatBirthAgePerm = nan(size(POW,2),size(POW,3));
    tStatMomsAgePerm  = nan(size(POW,2),size(POW,3));
    tStatGAPerm       = nan(size(POW,2),size(POW,3));
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


            if any(ismember(lmeout.Coefficients.Name,'GA'))
                PvalsGAPerm(itime,ifreq) = lmeout.Coefficients.pValue(strcmp(lmeout.Coefficients.Name,'GA'));
                tStatGAPerm(itime,ifreq) = lmeout.Coefficients.tStat(strcmp(lmeout.Coefficients.Name,'GA'));
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

    [~,~,GA_ssPerm_neg(iperm,:,:)] = TFCE(PvalsGAPerm,tStatGAPerm,'negative',0);
    [~,~,AgeMA_ssPerm_neg(iperm,:,:)] = TFCE(PvalsAgeMAPerm,tStatAgeMAPerm,'negative',0);
    [~,~,MABA_ssPerm_neg(iperm,:,:)] = TFCE(PvalsMABAPerm,tStatMABAPerm,'negative',0);
    [~,~,MomsAge_ssPerm_neg(iperm,:,:)] = TFCE(PvalsMomsAgePerm,tStatMomsAgePerm,'negative',0);
    [~,~,BirthAge_ssPerm_neg(iperm,:,:)] = TFCE(PvalsBirthAgePerm,tStatBirthAgePerm,'negative',0);
    [~,~,Sex_ssPerm_neg(iperm,:,:)] = TFCE(PvalsSexPerm,tStatSexPerm,'negative',0);
    [~,~,SDNN_ssPerm_neg(iperm,:,:)] = TFCE(PvalsSDNNPerm,tStatSDNNPerm,'negative',0);
    [~,~,Rule_ssPerm_neg(iperm,:,:)] = TFCE(PvalsRulePerm,tStatRulePerm,'negative',0);
    [~,~,Stimulus_ssPerm_neg(iperm,:,:)] = TFCE(PvalsStimulusPerm,tStatStimulusPerm,'negative',0);
    [~,~,BMI_ssPerm_neg(iperm,:,:)] = TFCE(PvalsBMIPerm,tStatBMIPerm,'negative',0)

    [~,~,GA_ssPerm_pos(iperm,:,:)] = TFCE(PvalsGAPerm,tStatGAPerm,'positive',0);
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
saveme = sprintf('fetal_perm_stats_TFCEs_%s',date);
save(saveme)

%% Gather results

% Results Aggregation:
% 
%     Several matrices (GAposN, MAposN, ..., BMInegN) are initialized to store the count of significant clusters for each variable and permutation.
%     The matrices are initialized with zeros and will be incremented based on the significance of the TFCE clusters.
% 
% Counting Significant Clusters:
% 
%     Two separate loops are used to iterate over each permutation (iperm).
%     For each permutation, the code compares the TFCE results for each variable (GA_statsum_pos, MomsAge_statsum_pos, ..., BMI_statsum_neg) with the corresponding permuted results (GA_ssPerm_pos, MomsAge_ssPerm_pos, ..., BMI_ssPerm_neg).
%     The count matrices (GAposN, MAposN, ..., BMInegN) are incremented based on whether a cluster in the original data is more extreme than the clusters in the permuted data.

GAposN = zeros(size(POW,2),size(POW,3));
MAposN = zeros(size(POW,2),size(POW,3));
BAposN = zeros(size(POW,2),size(POW,3));
SEXposN = zeros(size(POW,2),size(POW,3));
SDNNposN = zeros(size(POW,2),size(POW,3));
STIMposN = zeros(size(POW,2),size(POW,3));
RuleposN = zeros(size(POW,2),size(POW,3));
MABAposN = zeros(size(POW,2),size(POW,3));
BMIposN = zeros(size(POW,2),size(POW,3));

GAnegN = zeros(size(POW,2),size(POW,3));
MAnegN = zeros(size(POW,2),size(POW,3));
BAnegN = zeros(size(POW,2),size(POW,3));
SEXnegN = zeros(size(POW,2),size(POW,3));
SDNNnegN = zeros(size(POW,2),size(POW,3));
STIMnegN = zeros(size(POW,2),size(POW,3));
RulenegN = zeros(size(POW,2),size(POW,3));
MABAnegN = zeros(size(POW,2),size(POW,3));
BMInegN = zeros(size(POW,2),size(POW,3));

for iperm = 1:Nperm
    GAposN  = GAposN +  (GA_statsum_pos > squeeze(GA_ssPerm_pos(iperm,:,:)));
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
    GAnegN = GAnegN +  (GA_statsum_neg > squeeze(GA_ssPerm_neg(iperm,:,:)));
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
%     Several variables (BAposPval, SDNNposPval, ..., BMInegPval) store the transformed p-values for positive and negative clusters.

tmp = (Nperm - BAposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
BAposPval = -log10(tmp);


tmp = (Nperm - SDNNposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
SDNNposPval = -log10(tmp);


tmp = (Nperm - GAposN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
GAposPval = -log10(tmp);


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
SDNNnegPval = -log10(tmp);


tmp = (Nperm - GAnegN)./Nperm;
tmp(tmp==0) = 1/(Nperm+1);
GAnegPval = -log10(tmp);


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
Pcrit = 0.05; % adjust this later
RWB = customcolormap_preset('red-white-blue');

%%%% SDNN/HRV %%%%
%     Plots figures related to SDNN/HRV (Standard Deviation of NN intervals/Heart Rate Variability) analysis.
%     Uses a custom colormap, adjusts color axis, and adds contours based on statistical significance.
%     Visualizes results over time and frequency.
%     Generates and saves the figure as an SVG file.

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
%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(SDNNtft,4,interpstyle))),[0.01 2],'linewidth',1.5,'color','k','linestyle',':')

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
makefighandsome
title('SDNN','fontsize',20)
print('-dsvg','./Figures/SDNNTFCEfetal')


%%%% Maternal AGE %%%%
%   ... do same as above, but for maternal age 

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

%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(MAtft,4,interpstyle))),[0.01 2],'linewidth',1.5,'color','k','linestyle',':')


plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
%set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Maternal age','fontsize',20)
print('-dsvg','./Figures/MATFCEfetal')


%%%% BMI %%%%
% ... now do plotting for maternal body mass index before pregnancy 

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

%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(BMItft,4,interpstyle))),[0.01 2],'linewidth',1.5,'color','k','linestyle',':')


plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
%set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Maternal BMI','fontsize',20)
print('-dsvg','./Figures/BMITFCEfetal')


%%% Birth AGE %%%%
% ... and so forth, you get the idea

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

%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(BAtft,4,interpstyle))),[0.01 2],'linewidth',1.5,'color','k','linestyle',':')


plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Birth age','fontsize',20)
print('-dsvg','./Figures/BATFCEfetal')


%%% GESTATIONAL AGE %%%%
myfigure2
plotme1 = flipud(rot90(interp2(GAposPval,4,interpstyle)));
plotme2 = flipud(rot90(interp2(GAnegPval,4,interpstyle)));
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


%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(GAtft,4,interpstyle))),[0.01 2],'linewidth',2,'color','k','linestyle',':')

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Age','fontsize',20)
print('-dsvg','./Figures/GATFCEfetal')


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


%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(AGEtft,4,interpstyle))),[0.01 2],'linewidth',2,'color','k','linestyle',':')

plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Mother''s age','fontsize',20)
print('-dsvg','./Figures/MATFCEfetal')


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

%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(SEXtft,4,interpstyle))),[0.01 2],'linewidth',1.5,'color','k','linestyle',':')


plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
%set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Sex','fontsize',20)
print('-dsvg','./Figures/SEXTFCEfetal')


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

%contour(plotme,[log10(Pcrit) -log10(Pcrit)],'linewidth',3,'color',	"#FF00FF")
contour(plotme,[2 -2],'linewidth',3,'color',	"#7E2F8E")
%contour(flipud(rot90(interp2(STIMtft,4,interpstyle))),[0.01 2],'linewidth',1.5,'color','k','linestyle',':')


plot(ones(1,size(plotme,1)).*xax(xidx(1)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(2)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(3)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
plot(ones(1,size(plotme,1)).*xax(xidx(4)),1:size(plotme,1),'k','linewidth',2,'linestyle','--')
colormap(RWB)
%set(h,'YMinorTick','on')
%h.YRuler.Scale = 'log';
makefighandsome
title('Stimulus','fontsize',20)
print('-dsvg','./Figures/STIMTFCEfetal')


