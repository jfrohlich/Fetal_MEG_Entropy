% Joel Frohlich
% University of Tuebingen

function[trH,W] = timeresentropy_preissl(signal,window,overlap)
% Time resolved entropy using CTW
% signal - EEG time series
% window - the size of the time window in samples
% overlap - overlap between windows (proportion, i.e., 0.5 = 50%)

if nargin < 3, overlap = 0.95; end
if nargin < 2, window = 200;   end

assert(overlap<1,'Overlap is set to greater than 100%!')
step = window*(1-overlap);
W = 1:step:length(signal)-window;
N = length(W);
trH = nan(1,N); 

start = 1;
stop = window;

for istep = 1:N
    %fprintf('.')
    X = signal(start:stop);
    trH(istep) = CTWEntropyRate(X);
    start = start + step;
    stop = stop + step;
end
%fprintf('\n')

end


