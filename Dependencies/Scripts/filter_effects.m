% Joel Frohlich
% University of Tuebingen
% This script generates the supplemental figure comparing different filter
% settings for neonatal MEG data. 

% The tables Neo1 and Neo2 have specific columns corresponding to the
% entropy measures being analyzed (CTW, LZC, PermEn32, etc.).

clearvars
close all

% Load two tables (Neo1 and Neo2) containing neonatal MEG data with 
% different lowpass filter cutoff frequencies 
Neo1 = readtable('NeonatalTable.csv');
Neo2 = readtable('NeonatalTable_10HzLP.csv');

% Calculate the Pearson correlation coefficients (rctw, rlzc, rp32, rp64,
% rmse, rse) between corresponding columns of Neo1 and Neo2 (CTW, LZC,
% PermEn32, PermEn64, MSE, SE).
rctw = corr(Neo1.CTW,Neo2.CTW);
rlzc = corr(Neo1.LZC,Neo2.LZC);
rp32 = corr(Neo1.PermEn32,Neo2.PermEn32);
rp64 = corr(Neo1.PermEn64,Neo2.PermEn64);
rmse = corr(Neo1.MSE,Neo2.MSE);
rse = corr(Neo1.SE,Neo2.SE);

%%% Individual Plots for Each Entropy Measure %%%

%% CTW

%     Scatter plot comparing CTW values between the two filter settings.
%     Linear regression line (mylsline) added to the plot.
%     Adjusted axis limits for better visualization.
%     Save the plot as 'CTW_filter_corr.png'.

myfigure
scatter(Neo2.CTW,Neo1.CTW);
mylsline
a = min([Neo2.CTW; Neo1.CTW]);
b = max([Neo2.CTW; Neo1.CTW]);
axis([0 b 0 b])
xlabel('1 - 10 Hz')
ylabel('1 - 15 Hz')
makefigpretty
title(sprintf('CTW, r = %1.2f',rctw),'fontsize',20)
print('-dpng','CTW_filter_corr')

%% LZC

% Do the same below for Lempel-Ziv

myfigure
scatter(Neo2.LZC,Neo1.LZC);
mylsline
a = min([Neo2.LZC; Neo1.LZC]);
b = max([Neo2.LZC; Neo1.LZC]);
axis([0 b 0 b])
xlabel('1 - 10 Hz')
ylabel('1 - 15 Hz')
makefigpretty
title(sprintf('LZC, r = %1.2f',rlzc),'fontsize',20)
print('-dpng','LZC_filter_corr')

%% PermEn32

% Do the same below for 32 ms permutation entropy

myfigure
scatter(Neo2.PermEn32,Neo1.PermEn32);
mylsline
a = min([Neo2.PermEn32; Neo1.PermEn32])*0.97;
b = max([Neo2.PermEn32; Neo1.PermEn32]);
axis([a b a b])
xlabel('1 - 10 Hz')
ylabel('1 - 15 Hz')
makefigpretty
title(sprintf('PermEn32, r = %1.2f',rp32),'fontsize',20)
print('-dpng','PermEn32_filter_corr')

%% PermEn64

% Do the same below for 64 ms permutation entropy

myfigure
scatter(Neo2.PermEn64,Neo1.PermEn64);
mylsline
a = min([Neo2.PermEn64; Neo1.PermEn64])*0.97;
b = max([Neo2.PermEn64; Neo1.PermEn64]);
axis([a b a b])
xlabel('1 - 10 Hz')
ylabel('1 - 15 Hz')
makefigpretty
title(sprintf('PermEn64, r = %1.2f',rp64),'fontsize',20)
print('-dpng','PermEn64_filter_corr')

%% mMSE

% Do the same for mMSE

myfigure
scatter(Neo2.MSE,Neo1.MSE);
mylsline
a = min([Neo2.MSE; Neo1.MSE]);
b = max([Neo2.MSE; Neo1.MSE]);
axis([0 b 0 b])
xlabel('1 - 10 Hz')
ylabel('1 - 15 Hz')
makefigpretty
title(sprintf('MSE, r = %1.2f',rmse),'fontsize',20)
print('-dpng','MSE_filter_corr')

%% mSampEn

% Do the same for mSampEn

myfigure
scatter(Neo2.SE,Neo1.SE);
mylsline
a = min([Neo2.SE; Neo1.SE]);
b = max([Neo2.SE; Neo1.SE]);
axis([0 b 0 b])
xlabel('1 - 10 Hz')
ylabel('1 - 15 Hz')
makefigpretty
title(sprintf('SE, r = %1.2f',rse),'fontsize',20)
print('-dpng','SE_filter_corr')



