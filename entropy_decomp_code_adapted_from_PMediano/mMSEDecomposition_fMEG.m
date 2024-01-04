function [ res ] = mMSEDecomposition_fMEG(x1, s1, x2, s2, nb_samples, fsample,group)
%% mMSEDECOMPOSITION Decomposes spectral, phasic, and interaction components
% of the mMSE difference between two datasets.
%
%   R = mMSEDECOMPOSITION(X1, X2), where X1,X2 are 1D cell arrays of 1-by-T time
%   series vectors, computes the average difference in mMSE(X1) - mMSE(X2) and
%   decomposes it in spectral, phasic, and interaction terms using phase
%   randomisation techniques.
%
%   R = mMSEDECOMPOSITION(X1, X2, N) uses N surrogate samples (default: 1000).
%
%   R = mMSEDECOMPOSITION(..., S) appends the results to pre-existing struct S.
%
% Results are returned in a struct R with all elements of the decomposition.
% These are:
%
%   - pe1: average mMSE of dataset X1
%   - pe2: average mMSE of dataset X2
%   - diff: total mMSE difference (X1 - X2)
%   - pai1: phase-amplitud component of mMSE in X1
%   - pai2: phase-amplitud component of mMSE in X2
%   - phase: phasic component of mMSE difference (X1 - X2)
%   - spec: spectral component of mMSE difference (X1 - X2)
%
% Reference:
%   Mediano, P. A., Rosas, F. E., Barrett, A. B., & Bor, D. (2020). Decomposing
%   spectral and phasic differences in non-linear features between datasets.
%   arXiv:2009.10015.
%   https://arxiv.org/pdf/2009.10015.pdf
%
% Adapted from code by Pedro Mediano, Oct 2019
% Modified by Joel Frohlich, August 2021
% Updated by Joel Frohlich for fetal MEG, June 2022

dbstop if error

if nargin < 7 
    group = 'fetal';
end

switch group
    case 'fetal'
        assert(s1 == s2,'Young and old data are from different subjects')
    otherwise
        % do nothing
end


%% Argument checks and parameter initialisation

% Note: in the version of this function that was written for Angelman, the
% last input argument didn't do anything; the data were downsampled to 125
% Hz using the get_sections helper routine.

if nargin < 3 || isempty(nb_samples)
    nb_samples = 200;
end
if nargin < 4  || isempty(fsample)
    fsample = 610.3516;
end

x1_valid = iscell(x1) && all(cellfun(@length, x1) == length(x1{1}));
x2_valid = iscell(x2) && all(cellfun(@length, x2) == length(x2{1}));
if ~(x1_valid && x2_valid)
    error('Inputs must be 1D cell arrays of equally-sized trials');
end
nb_windows_1 = length(x1);
nb_windows_2 = length(x2);

nscale = 20; % do 20 timescales

% Configuration for mMSE
cfg = [];
cfg.fs = fsample;
cfg.scales = 1:nscale;
cfg.type = 'Xie';
cfg.r = 0.15;
cfg.dynr = true;
cfg.gpu = false;
SE = @(u,v) MSE(u,v);



isoctave = exist('OCTAVE_VERSION', 'builtin');
if isoctave
    % We're in Octave -- disable parallelisation by default
    nb_workers = 0;
else
    % We're in Matlab -- enable parfor only if not run in a worker thread already
    isOnWorker = ~isempty(getCurrentTask());
    if isOnWorker
        nb_workers = 0;
    else
        nb_workers = 32;
    end
end

% Set function handle for symmetric IFFT depending on Octave/Matlab version
if isoctave
    symm_ifft = @(z) real(ifft(z));
else
    symm_ifft = @(z) ifft(z, 'symmetric');
end

%% loop through lag values


mse1 = nan;
mse2 = nan;

mse1_f1  = nan([nb_samples, 1]);
mse1_f12 = nan([nb_samples, 1]);
mse2_f2  = nan([nb_samples, 1]);
mse2_f12 = nan([nb_samples, 1]);

progcnt = 0; % progress counter


%% Pre-compute spectra, phases and true SE
psd1  = cellfun(@(v)   abs(fft(v)), x1, 'UniformOutput', false);
phi1  = cellfun(@(v) angle(fft(v)), x1, 'UniformOutput', false);
psd2  = cellfun(@(v)   abs(fft(v)), x2, 'UniformOutput', false);
phi2  = cellfun(@(v) angle(fft(v)), x2, 'UniformOutput', false);
% phi_all = [phi1; phi2];
phi_all = [phi1 phi2]; % JF edit - 07/02/21
assert(all(size(x1)==size(x2)))
CFG = repmat({cfg},size(x1,1),size(x1,2));
%out = cellfun(SE, x1, CFG);
mse1 = mean(cellfun(SE, x1, CFG));
mse2 = mean(cellfun(SE, x2, CFG));


%% Initialise arrays and compute surrogate SE values
%{
            Notation:
              x1     -- actual time series of condition 1
              x1_f1  -- time series with spectra from condition 1 and phases swapped within condition 1
              x1_f12 -- time series with spectra from condition 1 and phases pooled from both conditions
%}


parfor (isamp=1:nb_samples, nb_workers)
    %for isamp = 1:nb_samples
    % Swapped within condition

    rng(isamp,'twister') % we must include this line for rng to be set inside the parfor loop
    phi1_shuf = phi1(randperm(nb_windows_1));
    phi2_shuf = phi2(randperm(nb_windows_2));
    x1_f1   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
    CFG = repmat({cfg},size(x1_f1,1),size(x1_f1,2));
    mse1_f1(isamp) = mean(cellfun(SE, x1_f1, CFG));

    x2_f2   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
    CFG = repmat({cfg},size(x2_f2,1),size(x2_f2,2));
    mse2_f2(isamp) = mean(cellfun(SE, x2_f2, CFG));

    % Swapped across conditions
    shuf_idx = randperm(nb_windows_1 + nb_windows_2);
    phi1_shuf = phi_all(shuf_idx(1:nb_windows_1));
    phi2_shuf = phi_all(shuf_idx(nb_windows_1 + 1:end));

    x1_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
    CFG = repmat({cfg},size(x1_f12,1),size(x1_f12,2));
    mse1_f12(isamp) = mean(cellfun(SE, x1_f12, CFG));

    x2_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
    CFG = repmat({cfg},size(x2_f12,1),size(x2_f12,2));
    mse2_f12(isamp) = mean(cellfun(SE, x2_f12, CFG));
%     if isamp == 50
%         fprintf('     Getting there ...\n')
%     elseif isamp == 100
%         fprintf('     Halfway there ...\n')
%     elseif isamp == 150
%         fprintf('     Almost there!\n')
%     end
end



%% Put results together and return

% Uncomment below to average before saving in the results structure
mse1_f1  = squeeze(mean(mse1_f1));
mse1_f12 = squeeze(mean(mse1_f12));
mse2_f2  = squeeze(mean(mse2_f2));
mse2_f12 = squeeze(mean(mse2_f12));


res.mse1 = mse1';
res.mse2 = mse2';
res.diff = mse1' - mse2';
res.pai1 = mse1' - mse1_f1;
res.pai2 = mse2' - mse2_f2;
phamse1 = mse1_f1 - mse1_f12;
phamse2 = mse2_f2 - mse2_f12;
res.phase  = phamse1 - phamse2;
res.spec = mse1_f12 - mse2_f12; % the phases are pooled, so difference must be due to spectra
res.interact = res.pai1 - res.pai2;


end

function[out] = MSE(in,cfg)
% wrapper function for multiscale entropy that averages the sample
% entropies of each time scale

tmp = ro_mse(zscore(in),cfg); % make sure to zscore!
out = mean(tmp);

end
