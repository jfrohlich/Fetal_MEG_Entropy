function [ res ] = CTWDecomposition(x1, s1, x2, s2, nb_samples, fsample, group)
%% CTWDECOMPOSITION Decomposes spectral, phasic, and interaction components
% of the CTW difference between two datasets.
%
%   R = CTWDECOMPOSITION(X1, X2), where X1,X2 are 1D cell arrays of 1-by-T time
%   series vectors, computes the average difference in CTW(X1) - CTW(X2) and
%   decomposes it in spectral, phasic, and interaction terms using phase
%   randomisation techniques.
%
%   R = CTWDECOMPOSITION(X1, X2, N) uses N surrogate samples (default: 1000).
%
%   R = CTWDECOMPOSITION(..., S) appends the results to pre-existing struct S.
%
% Results are returned in a struct R with all elements of the decomposition.
% These are:
%
%   - ctw1: average CTW of dataset X1
%   - ctw2: average CTW of dataset X2
%   - diff: total CTW difference (X1 - X2)
%   - pai1: phase-amplitud component of CTW in X1
%   - pai2: phase-amplitud component of CTW in X2
%   - phase: phasic component of CTW difference (X1 - X2)
%   - spec: spectral component of CTW difference (X1 - X2)
%
% Reference:
%   Mediano, P. A., Rosas, F. E., Barrett, A. B., & Bor, D. (2020). Decomposing
%   spectral and phasic differences in non-linear features between datasets.
%   arXiv:2009.10015.
%
% Pedro Mediano, Oct 2019
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
if nargin < 5 || isempty(nb_samples)
  nb_samples = 1000;
end
x1_valid = iscell(x1) && all(cellfun(@length, x1) == length(x1{1}));
x2_valid = iscell(x2) && all(cellfun(@length, x2) == length(x2{1}));
if ~(x1_valid && x2_valid)
  error('Inputs must be 1D cell arrays of equally-sized trials');
end
nb_windows_1 = length(x1);
nb_windows_2 = length(x2);
CTW = @(v) CTWEntropyRate(v); % - JF edit, 07/05/21


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
    nb_workers = 4;
  end
end

% Set function handle for symmetric IFFT depending on Octave/Matlab version
if isoctave
    symm_ifft = @(z) real(ifft(z));
else
    symm_ifft = @(z) ifft(z, 'symmetric');
end


%% Pre-compute spectra, phases and true CTW
psd1  = cellfun(@(v)   abs(fft(v)), x1, 'UniformOutput', false);
phi1  = cellfun(@(v) angle(fft(v)), x1, 'UniformOutput', false);
psd2  = cellfun(@(v)   abs(fft(v)), x2, 'UniformOutput', false);
phi2  = cellfun(@(v) angle(fft(v)), x2, 'UniformOutput', false);
% phi_all = [phi1; phi2];
phi_all = [phi1 phi2]; % JF edit - 07/02/21

ctw1 = mean(cellfun(CTW, x1));
ctw2 = mean(cellfun(CTW, x2));


%% Initialise arrays and compute surrogate CTW values
%{
Notation:
  x1     -- actual time series of condition 1
  x1_f1  -- time series with spectra from condition 1 and phases swapped within condition 1
  x1_f12 -- time series with spectra from condition 1 and phases pooled from both conditions
%}
ctw1_f1  = zeros([nb_samples, 1]);
ctw1_f12 = zeros([nb_samples, 1]);
ctw2_f2  = zeros([nb_samples, 1]);
ctw2_f12 = zeros([nb_samples, 1]);

parfor (i=1:nb_samples, nb_workers)

  % Swapped within condition
  phi1_shuf = phi1(randperm(nb_windows_1));
  phi2_shuf = phi2(randperm(nb_windows_2));
  x1_f1   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
  ctw1_f1(i)  = mean(cellfun(CTW, x1_f1));

  x2_f2   = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
  ctw2_f2(i)  = mean(cellfun(CTW, x2_f2));

  % Swapped across conditions
  shuf_idx = randperm(nb_windows_1 + nb_windows_2);
  phi1_shuf = phi_all(shuf_idx(1:nb_windows_1));
  phi2_shuf = phi_all(shuf_idx(nb_windows_1 + 1:end));

  x1_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd1, phi1_shuf, 'UniformOutput', false);
  ctw1_f12(i) = mean(cellfun(CTW, x1_f12));

  x2_f12  = cellfun(@(s, p) symm_ifft(s.*exp(1j.*p)), psd2, phi2_shuf, 'UniformOutput', false);
  ctw2_f12(i) = mean(cellfun(CTW, x2_f12));
  if mod(i,10)==0, fprintf('.'), end
end


%% Put results together and return
ctw1_f1  = mean(ctw1_f1);
ctw1_f12 = mean(ctw1_f12);
ctw2_f2  = mean(ctw2_f2);
ctw2_f12 = mean(ctw2_f12);

res.ctw1 = ctw1;
res.ctw2 = ctw2;
res.diff = ctw1 - ctw2;
res.pai1 = ctw1 - ctw1_f1;
res.pai2 = ctw2 - ctw2_f2;
phase1 = ctw1_f1 - ctw1_f12;
phase2 = ctw2_f2 - ctw2_f12;
res.phase  = phase1 - phase2;
res.spec = ctw1_f12 - ctw2_f12; % the phases are pooled, so difference must be due to spectra
res.interact = res.pai1 - res.pai2;


%{
  Old results struct. Some of these quantities aren't of interest to most users,
  so have been removed of return struct. This is left here for now for possible
  future refactoring.
%}
% res.ctw1 = ctw1;
% res.ctw2 = ctw2;
% res.ctw1_f1 =  ctw1_f1;
% res.ctw2_f2 =  ctw2_f2;
% res.ctw1_f12 = ctw1_f12;
% res.ctw2_f12 = ctw2_f12;
% res.diff = ctw1 - ctw2;
% res.pai1 = ctw1 - ctw1_f1;
% res.pai2 = ctw2 - ctw2_f2;
% res.phase1 = ctw1_f1 - ctw1_f12;
% res.phase2 = ctw2_f2 - ctw2_f12;
% res.phase  = res.phase1 - res.phase2;
% res.spec = ctw1_f12 - ctw2_f12;

