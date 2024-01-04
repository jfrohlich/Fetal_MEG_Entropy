function [ res ] = PermEnDecomposition_fMEG(x1, s1, x2, s2, nb_samples, fsample, group)
%% PermEnDECOMPOSITION Decomposes spectral, phasic, and interaction components
% of the PermEn difference between two datasets.
%
%   R = PermEnDECOMPOSITION(X1, X2), where X1,X2 are 1D cell arrays of 1-by-T time
%   series vectors, computes the average difference in PermEn(X1) - PermEn(X2) and
%   decomposes it in spectral, phasic, and interaction terms using phase
%   randomisation techniques.
%
%   R = PermEnDECOMPOSITION(X1, X2, N) uses N surrogate samples (default: 1000).
%
%   R = PermEnDECOMPOSITION(..., S) appends the results to pre-existing struct S.
%
% Results are returned in a struct R with all elements of the decomposition.
% These are:
%
%   - pe1: average PermEn of dataset X1
%   - pe2: average PermEn of dataset X2
%   - diff: total PermEn difference (X1 - X2)
%   - pai1: phase-amplitud component of PermEn in X1
%   - pai2: phase-amplitud component of PermEn in X2
%   - phase: phasic component of PermEn difference (X1 - X2)
%   - spec: spectral component of PermEn difference (X1 - X2)
%
% Reference:
%   Mediano, P. A., Rosas, F. E., Barrett, A. B., & Bor, D. (2020). Decomposing
%   spectral and phasic differences in non-linear features between datasets.
%   arXiv:2009.10015.
%   https://arxiv.org/pdf/2009.10015.pdf
%
% Adapted from code by Pedro Mediano, Oct 2019
% Modified by Joel Frohlich, August 2021
% Updated by Joel Frohlich for fetal MEG, September 2022

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

PE = @(u,v,w) PermEn(u,v,w);



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

lags = [32 64]; % lags for permutation entropy in milliseconds
embed = [3]; % embedding dimension for PermEn

pe1 = nan(length(lags), 1);
pe2 = nan(length(lags), 1);

pe1_f1  = nan([nb_samples, length(lags)]);
pe1_f12 = nan([nb_samples, length(lags)]);
pe2_f2  = nan([nb_samples, length(lags)]);
pe2_f12 = nan([nb_samples, length(lags)]);

progcnt = 0; % progress counter

    for ilag = 1:length(lags)
        %% Pre-compute spectra, phases and true PE
        psd1  = cellfun(@(v)   abs(fft(v)), x1, 'UniformOutput', false);
        phi1  = cellfun(@(v) angle(fft(v)), x1, 'UniformOutput', false);
        psd2  = cellfun(@(v)   abs(fft(v)), x2, 'UniformOutput', false);
        phi2  = cellfun(@(v) angle(fft(v)), x2, 'UniformOutput', false);
        % phi_all = [phi1; phi2];
        phi_all = [phi1 phi2]; % JF edit - 07/02/21
        assert(all(size(x1)==size(x2)))
        tau = repmat({round((lags(ilag)/1000)/(1/fsample))},size(x1,1),size(x1,2));
        m   = repmat({embed},size(x1,1),size(x1,2));
        pe1(ilag) = mean(cellfun(PE, x1, m, tau));
        pe2(ilag) = mean(cellfun(PE, x2, m, tau));


        %% Initialise arrays and compute surrogate PE values
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
            tau = repmat({(lags(ilag)/1000)/(1/fsample)},size(x1_f1,1),size(x1_f1,2));
            m   = repmat({embed},size(x1_f1,1),size(x1_f1,2));
            pe1_f1(isamp,ilag) = mean(cellfun(PE, x1_f1, m, tau));

            x2_f2   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
            tau = repmat({(lags(ilag)/1000)/(1/fsample)},size(x2_f2,1),size(x2_f2,2));
            m   = repmat({embed},size(x2_f2,1),size(x2_f2,2));
            pe2_f2(isamp,ilag) = mean(cellfun(PE, x2_f2, m, tau));

            % Swapped across conditions
            shuf_idx = randperm(nb_windows_1 + nb_windows_2);
            phi1_shuf = phi_all(shuf_idx(1:nb_windows_1));
            phi2_shuf = phi_all(shuf_idx(nb_windows_1 + 1:end));

            x1_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
            tau = repmat({(lags(ilag)/1000)/(1/fsample)},size(x1_f12,1),size(x1_f12,2));
            m = repmat({embed},size(x1_f12,1),size(x1_f12,2));
            pe1_f12(isamp,ilag) = mean(cellfun(PE, x1_f12, m, tau));

            x2_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
            tau = repmat({(lags(ilag)/1000)/(1/fsample)},size(x2_f12,1),size(x2_f12,2));
            tau = repmat({embed},size(x2_f12,1),size(x2_f12,2));
            pe2_f12(isamp,ilag) = mean(cellfun(PE, x2_f12, m, tau));
%             if isamp == 50
%                 fprintf('     Getting there ...\n')
%             elseif isamp == 100
%                 fprintf('     Halfway there ...\n')
%             elseif isamp == 150
%                 fprintf('     Almost there!\n')
%             end
        end
        progcnt = progcnt + 1;
        %if mod(progcnt,2) == 0 % sets how often progress is reported
            fprintf('\n%2.1f%% complete for %s ...\n',...
                (progcnt/(length(lags)))*100, s1)
        %end
    end

%% Put results together and return

% Uncomment below to average before saving in the results structure
pe1_f1  = squeeze(mean(pe1_f1));
pe1_f12 = squeeze(mean(pe1_f12));
pe2_f2  = squeeze(mean(pe2_f2));
pe2_f12 = squeeze(mean(pe2_f12));


res.pe1 = pe1';
res.pe2 = pe2';
res.diff = pe1' - pe2';
res.pai1 = pe1' - pe1_f1;
res.pai2 = pe2' - pe2_f2;
phase1 = pe1_f1 - pe1_f12;
phase2 = pe2_f2 - pe2_f12;
res.phase  = phase1 - phase2;
res.spec = pe1_f12 - pe2_f12; % the phases are pooled, so difference must be due to spectra
res.interact = res.pai1 - res.pai2; 
res.lags = lags;


end

