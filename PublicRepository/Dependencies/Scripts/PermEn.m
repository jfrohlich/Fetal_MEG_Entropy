
function [PermEn,PermEn_unnormed] = PermEn(y,m,t)
%  Calculate the permutation entropy
%  Input:   y: time series;
%           m: order of permuation entropy
%           t: delay time of permuation entropy, 
% Output: 
%           pe:    permuation entropy
%           hist:  the histogram for the order distribution
% Edited by Joel Frohlich 25.02.2022
% I modified this from code downloaded here https://de.mathworks.com/matlabcentral/fileexchange/37289-permutation-entropy
%Ref: G Ouyang, J Li, X Liu, X Li, Dynamic Characteristics of Absence EEG Recordings with Multiscale Permutation %     %                             Entropy Analysis, Epilepsy Research, doi: 10.1016/j.eplepsyres.2012.11.003
%     X Li, G Ouyang, D Richards, Predictability analysis of absence seizures with permutation entropy, Epilepsy %     %                            Research,  Vol. 77pp. 70-74, 2007

warning('off', 'all')

ly = length(y);
permlist = perms(1:m);
c = nan(1,length(permlist));
    
 for j=1:ly-t*(m-1)
     [~,iv]=sort(y(j:t:j+t*(m-1)));
     for jj=1:length(permlist)
         if (abs(permlist(jj,:)-iv))==0
             if isnan(c(jj))
                 c(jj) = 1;
             else
                 c(jj) = c(jj) + 1;
             end
         end
     end
 end
%hist = c;
 
p = c/nansum(c);
pe = -nansum(p .* log(p));
PermEn = pe/log(factorial(m));
PermEn_unnormed = pe;

end