% I downloaded this code by Laurie from Mathworks File Exchange (https://de.mathworks.com/matlabcentral/fileexchange/45966-compare-two-proportions-chi-square)

% Copyright (c) 2014, laurie
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
%
% * Neither the name of <organization> nor the names of its
%   contributors may be used to endorse or promote products derived from this
%   software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [h,p, chi2stat,df] = prop_test(X , N, correct)

% [h,p, chi2stat,df] = prop_test(X , N, correct)
% 
%
% A simple Chi-square test to compare two proportions
% It is a 2 sided test with alpha=0.05
%
% Input:
% X = vector with number of success for each sample (e.g. [20 22])
% N = vector of total counts for each sample (e.g. [48 29])
% correct = true/false : Yates continuity correction for small samples?
%
% Output:
% h = hypothesis (H1/H0)
% p = p value
% chi2stat= Chi-square value
% df = degrees of freedom (always equal to 1: 2 samples)
%
% Needs chi2cdf from the Statistics toolbox
% Inspired by prop.test() in "R" but much more basic
%
% Example: [h,p,chi]=prop_test([20 22],[48 29], true)
% The above example tests if 20/48 differs from 22/29 using Yate's correction

if (length(X)~= 2)||(length(X)~=length(N))
    disp('Error: bad vector length')
elseif (X(1)>N(1))|| (X(2)>N(2))
    disp('Error: bad counts (X>N)')
else  
    df=1; % 2 samples
    
    % Observed data
    n1 = X(1);
    n2 = X(2);
    N1 = N(1);
    N2 = N(2);
    
    % Pooled estimate of proportion
    p0 = (n1+n2) / (N1+N2);
    
    % Expected counts under H0 (null hypothesis)
    n10 = N1 * p0;
    n20 = N2 * p0;
    observed = [n1 N1-n1 n2 N2-n2];
    expected = [n10 N1-n10 n20 N2-n20];
    
    if correct == false
        % Standard Chi-square test
        chi2stat = sum((observed-expected).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
    else
        % Yates continuity correction        
        chi2stat = sum((abs(observed - expected) - 0.5).^2 ./ expected);
        p = 1 - chi2cdf(chi2stat,1);
    end
    
    h=0;
    if p<0.05
        h=1;
    end
end
end
