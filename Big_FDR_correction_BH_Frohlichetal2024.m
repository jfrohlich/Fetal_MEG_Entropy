% Joel Frohlich
% University of Tuebingen
% This is the script that loads everything and does one huge FDR correction

% Clear all variables
clearvars

%% Begin by loading all data from CSV files into tables. Then, separate
% the data into different tables based on fetal and neonatal studies.

%%% Load fetal tables containing p-values 

Tf1 = readtable('FetalModel.csv');
Tf2 = readtable('FetalModelDyanmics.csv');
Tf3 = readtable('FetalSurrogate.csv');
Tf4 = readtable('fMEGEntropyDecomp.csv');

Pf1_GA   = Tf1.GA_Pvalue;
Pf1_Sex  = Tf1.Sex_Pvalue;
Pf1_SexGA = Tf1.SexGA_Pvalue;
Pf1_MA = Tf1.MomsAge_Pvalue;
Pf2_GA   = Tf2.GA_Pvalue;
Pf2_Sex  = Tf2.Sex_Pvalue;
Pf2_SexGA = Tf2.SexGA_Pvalue;
Pf3      = Tf3.Sur_Pvalue;
Pf4      = Tf4.P_values;


%%% Load neonatal tables containing p-values

Tn1 = readtable('NeonatalModel.csv');
Tn2 = readtable('NeonatalSurrogate.csv');
Tn3 = readtable('neonate_decomp.csv');
Tn4 = readtable('NeonatalModelDyanmics.csv');

Pn1_age = Tn1.PvalsAGE;
Pn2 = Tn2.Sur_Pvalue;
Pn3 = Tn3.P_values;
Pn4_age = Tn4.age_Pvalue;

%% P-values from different tables are combined into a single vector P. 

P = [Pf1_GA; Pf1_Sex; Pf1_SexGA; Pf1_MA; Pf2_GA; Pf3; Pf4; ...
    Pn1_age; Pn2; Pn3; Pn4_age];

%% Call R from MATLAB to perform a Benjamini-Hochberg False Discovery Rate 
% (FDR) correction on the combined P-values (P). We compare the FDR 
% correction results obtained from R (p_adjusted1) and MATLAB's mafdr 
% function (p_adjusted2) to ensure agreement.

% Call R from Matlab
x_str = sprintf('%.2f,', P);
x_str = x_str(1:end-1);

cmd = sprintf('Rscript -e "x <- c(%s); p_adjusted <- p.adjust(x, method = \\"BH\\"); write.table(p_adjusted, stdout(), row.names = FALSE, col.names = FALSE, quote = FALSE)"', x_str);
[status, result] = system(cmd);
p_adjusted1 = str2num(result);

p_adjusted2 = mafdr(P,'bhfdr',1); 

assert(all(p_adjusted1<0.05==p_adjusted2<0.05),'R and Matlab don''t agree')

%% Displaying results

% The combined P-values, FDR-corrected results, and logical indicator of significance are displayed using the vpa function.
% The critical P-value for significance (critP) is determined.
% Logical indicators of survival after FDR correction are added to the original tables.
% The tables are converted to LaTeX format using functions (table2latex and roundtable).

vpa([P p_adjusted1<0.05],10)

idx = p_adjusted1 < 0.05;
critP = max(P(idx))


Ta.GA_xsurvive_FDR = logical2str(Ta.GA_Pvalue <= critP);
Ta.Sex_xsurvive_FDR = logical2str(Ta.Sex_Pvalue <= critP);
Ta.BornYet_xsurvive_FDR = logical2str(Ta.BornYet_Pvalue <= critP);
table2latex(roundtable(Ta), 'AllModelLaTeX')

Tf1.GA_xsurvive_FDR = logical2str(Tf1.GA_Pvalue <= critP);
Tf1.Sex_xsurvive_FDR = logical2str(Tf1.Sex_Pvalue <= critP);
Tf1.SexGA_xsurvive_FDR = logical2str(Tf1.SexGA_Pvalue <= critP);
Tf1.MomsAge_xsurvive_FDR = logical2str(Tf1.MomsAge_Pvalue <= critP);
table2latex(roundtable(Tf1), 'FetalModelLaTeX')

Tf2.GA_xsurvive_FDR = Tf2.GA_Pvalue <= critP;
Tf2.Sex_xsurvive_FDR = Tf2.Sex_Pvalue <= critP;
Tf2.SexGA_xsurvive_FDR = Tf2.SexGA_Pvalue <= critP;
table2latex(roundtable(Tf2), 'FetalDynamicsLaTeX')

Tf3.survive_FDR = Tf3.Sur_Pvalue <= critP;
table2latex(roundtable(Tf3), 'FetalSurrogateLaTeX')

Tf4.survive_FDR = Tf4.P_values <= critP;
table2latex(roundtable(Tf4), 'FetalEntropyDecompLaTeX')

Tn1.age_xsurvive_FDR = logical2str(Tn1.age_Pvalue <= critP);
table2latex(roundtable(Tn1), 'NeonatalModelLaTeX')

Tn2.Sur_xsurvive_FDR = logical2str(Tn2.Sur_Pvalue <= critP);
table2latex(roundtable(Tn2), 'NeonatalSurrogateLaTeX')

Tn3.survive_FDR = logical2str(Tn3.P_values <= critP);
table2latex(roundtable([Tn3]), 'NeonatalEntropyDecompLaTeX')

Tn4.age_xsurvive_FDR = logical2str(Tn4.age_Pvalue <= critP);
table2latex(roundtable(Tn4), 'NeonatalModelDynamicsLaTeX')


%% Helper functions

% Here we have two helper functions:
% 
%     roundtable: Rounds numerical values in a table and sorts columns.
%     logical2str: Converts logical values to string representation ('TRUE' or 'FALSE').

function[Tout] = roundtable(Tin)

% This function processes a table by rounding numeric values, sorting 
% columns, transposing the table (if applicable), and modifying string 
% entries for better readability.

% Input:
% 
%     Tin: Input table to be processed.
% 
% Output:
% 
%     Tout: Processed table with rounded numerical values, sorted columns,
%     transposed if necessary, and modified string entries.
% 
% Steps:
% 
%     Initialization:
%         Tout is initially set equal to Tin.
% 
%     Rounding Numeric Values:
%         Iterate through each column of Tin.
%         If the column contains numeric values (either double or single), 
%         round each element to 3 significant digits and update the 
%         corresponding column in Tout.
% 
%     Sort Columns:
%         Sort columns of Tout (excluding the first column) based on 
%         variable names.
% 
%     Transpose Table (if applicable):
%         If Tout has more than one row, transpose the table using the 
%         rows2vars function.
% 
%     String Replacement:
%         Replace underscores with spaces in the first column of Tout.
%         Replace specific string patterns ('xsurvive FDR', 'SexGA', 
%         'MomsAge', 'BornYet') with more readable alternatives.
% 
%     Output:
%         Return the processed table Tout.

Tout = Tin;

for icol = 1:size(Tin,2)
    if strcmp(class(Tin{:,icol}),'double') || strcmp(class(Tin{:,icol}),'single')
        Tout{:,icol} = round(Tin{:,icol},3,'significant');
    end
end


% sort the columns as desired
[~,NDX] = sort(Tout.Properties.VariableNames(2:end));
Tout = Tout(:,[1 NDX+1]);

if size(Tout,1) > 1
    % transpose table
    tmp = rows2vars(Tout);
    varnames = tmp{1,1:end};
    tmp(1,:) = [];
    tmp.Properties.VariableNames = varnames;
    Tout = tmp;
else
    Tout = rows2vars(Tout);
end

% Fix strings
for irow = 1:size(Tout,1)
    Tout{irow,1} = strrep(Tout{irow,1}, '_', ' '); 
    Tout{irow,1} = strrep(Tout{irow,1}, 'xsurvive FDR', 'survive FDR'); 
    Tout{irow,1} = strrep(Tout{irow,1}, 'SexGA', 'Sex x GA');
    Tout{irow,1} = strrep(Tout{irow,1}, 'MomsAge', 'MA');
    Tout{irow,1} = strrep(Tout{irow,1}, 'BornYet', 'Pre/Post Birth');
end


end



function[str] = logical2str(x)
% 
% This function takes an input array of logical values (x) and converts 
% each element to a corresponding string value. The output is a cell array 
% of strings (str), where each element represents the logical value in 
% string form.
% 
%     Input:
%         x: Array of logical values.
% 
%     Output:
%         str: Cell array of strings representing the logical values in x.
% 
%     Processing:
%         The function initializes a cell array str with the same size as the input x.
%         It then iterates through each element of x.
%         For each element, if it is true, the corresponding string in str is set to 'TRUE'. If it is false, the string is set to 'FALSE'.
%         The final cell array str contains string representations of the logical values in x.

str = cell(size(x));
for i = 1:numel(x)
    if x(i)
        str{i} = 'TRUE';
    else
        str{i} = 'FALSE';
    end
end
end
