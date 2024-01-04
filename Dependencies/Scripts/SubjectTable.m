% Joel Frohlich
% University of Tuebingen

%%% SUMMARY %%% 
% %     Generate a comprehensive summary table for all subjects using data from 'Alltable.mat'.
% 
% %     Data Loading: Load the dataset 'Alltable.mat'.
% 
% %     Subject Identification:
% %         Extract unique subject IDs from the dataset.
% 
% %     Data Processing:
% %         Iterate through each subject ID.
% %         For each subject, gather relevant information (number of visits, postmenstrual ages, mother's age, mother's BMI, sex, and type).
% 
% %     Type Classification:
% %         Categorize subjects as Fetal, Neonatal, or Both based on the 'BornYet' variable.
% 
% %     Descriptive stats:
% %         Calculate the number of visits, average time between visits, and mean values for postmenstrual ages, mother's age, and mother's BMI.
% 
% %     Table Creation: Assemble a summary table ('Tout') with organized data and appropriate variable names.
% 
% %     Table Export: Save the generated table as both a CSV and an Excel file named 'SubjectTable.csv' and 'SubjectTable.xlsx', respectively.


%% Make supplmental table of subjects
load Alltable.mat

IDs = unique(Alltable.ID);
N = nan(1,length(IDs));
ages = cell(1, length(IDs));
type = cell(1, length(IDs));
MA = nan(1,length(IDs));
Mbmi = nan(1,length(IDs));
sex = cell(1, length(IDs));

for isub = 1:length(IDs)
    idx = Alltable.ID == IDs(isub);
    T = Alltable(idx,:);
    N(isub) = length(unique(T.GA));
    ages{isub} = unique(T.GA)';
    MA(isub) = unique(T.MomsAge);
    bmi = unique(T.BMI);
    bmi(bmi==0) = [];
    try
        Mbmi(isub) = bmi;
    catch
        Mbmi(isub) = nan;
    end
    
    switch unique(T.XChrom)
        case 1
            sex{isub} = 'M';
        case 2 
            sex{isub} = 'F';
        otherwise
            error('Bad sex coding')
    end
    
    tmp = length(unique(T.BornYet));
    switch tmp
        case 1
            if unique(T.BornYet) == 0
                type{isub} = 'Fetal';
            elseif unique(T.BornYet) == 1
                type{isub} = 'Neonatal';
            else
                error('Bad coding')
            end
        case 2
            type{isub} = 'Both';
        otherwise
            error('Bad coding')
    end
end


avgdif = nan(1,length(ages));
% get mean time between visits
for irow = 1:length(ages)
    avgdif(irow) = mean(diff(ages{irow}));
end

Tout = table(IDs,N',avgdif',type',sex',ages',MA',Mbmi','VariableNames',...
    {'Subject ID','Number of visits','Average time between visits (weeks)','Type','Sex','Postmenstrual ages at recordings (weeks)','Mother''s age (years)','Mother''s BMI'})

writetable(Tout,'./SubjectTable.csv')
writetable(Tout,'./SubjectTable.xlsx')




