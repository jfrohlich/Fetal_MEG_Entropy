% This script fixes errors in decoding subject IDs

load time_traces_data_manuscript
% Correct 0s written as letter o
data_tab.ID(find(data_tab.ID=="Co12-2-1A")) = "C012-2-1A";

N = size(data_tab,1) - 1; % leave out last row, doens't contain ID

old = nan(1,N); % old IDs
new = nan(1,N); % new (correct) IDs

for i = 1:N
    try
        ID = sscanf(data_tab.ID{i},'C00%i');
        assert(~isempty(ID),'Failed to extract ID')
        new(i) = ID;
    catch
        ID = sscanf(data_tab.ID{i},'C0%i');
        assert(~isempty(ID),'Failed to extract ID')
        new(i) = ID;
    end
    old(i) = sscanf(data_tab.ID{i},'C%i');
end

T = table(old',new','VariableNames',{'OldID','NewID'});

writetable(T, 'NeonateSubjectDecoder.csv')

