% Joel Frohlich 
% University of Tuebingen
% Generates the ERF response template for newborns

% Okay, what is the script doing in this repo? After all, we didn't use its
% results in the manuscript.
%
% Once upon a time, we had the idea to try looking at the entropy of
% signals after computing a grand-averaged ERF from all newborns and then
% subtracting this template from each recording. The idea was to use this
% approach to look at the background signal rather than the evoked
% response. It turns out this was a bad idea and we abandoned it. But this
% code is still included in the code because its output got grandfathered
% into later code as input.
%
% If nothing else, you may still find this script useful for plotting
% grand-averaged ERFs across recordings ...

%%% SUMMARY

%     Loading Data:
%         Close all figures, clear variables, and load the MEG data from 'time_traces_data_manuscript.mat'.
% 
%     Data Filtering:
%         Initialize empty arrays ERFdev, ERFstd, IDdev, and IDstd.
%         Iterate through the rows of 'data_tab' and select rows based on conditions.
%         Repair a string in the 'ID' column if it starts with 'Co12'.
%         Populate ERFdev, ERFstd, IDdev, and IDstd arrays with corresponding data from 'data_tab'.
%         Replace strings in IDdev and IDstd with unique numerical values.
% 
%     Calculate ERF Differences:
%         Compute the difference array ERFdif by subtracting ERFstd from ERFdev.
% 
%     Statistical Analysis:
%         Set transparency transp to 0.3.
%         Create a time vector t from -200 to 3000 ms with the length of ERFdev.
%         Perform t-tests to obtain confidence intervals (CIdev, CIstd, CIdif) for ERFdev, ERFstd, and ERFdif.
% 
%     Plotting and Visualization:
%         Create three separate plots for Deviant Response, Standard Response, and Deviant - Standard Response.
%         Use nestedmean to calculate the average ERF for each condition.
%         Plot the average ERF and fill the area between confidence intervals.
%         Add dashed lines at specific time points and set axis labels, legends, and titles for each plot.
%         Adjust the appearance of the figures using custom functions (myfigure2, makefighandsome).
% 
%     Save Figures:
%         Save each figure as a PNG and SVG file in the './Figures/' directory.
% 
%     Save Results:
%         Save ERFdev, ERFstd, ERFdif, IDdev, and IDstd in a file named 'Newborn_ERFs.mat'.


close all
clearvars
load time_traces_data_manuscript.mat

ERFdev = [];
ERFstd = [];
IDdev = [];
IDstd = [];

for irow = 1:size(data_tab,1)
    if strcmp('Co12',data_tab.ID{irow}(1:4))
        data_tab.ID{irow}(1:4) = 'C012'; % repair string
    end
    if contains(data_tab.condition{irow},'global_dev')
        ERFdev = [ERFdev; data_tab.data(irow,:)];
        IDdev = [IDdev; data_tab.ID{irow}(1:4)];
    elseif ~strcmp(data_tab.condition{irow},'time')
        ERFstd = [ERFstd; data_tab.data(irow,:)];
        IDstd = [IDstd; data_tab.ID{irow}(1:4)];
    end
end

% Replace strings with unique numbers
[~,~,IDdev] = unique(IDdev,'rows');
[~,~,IDstd] = unique(IDstd,'rows');

ERFdif = ERFdev-ERFstd;

transp = 0.3;
t = linspace(-200,3000,length(ERFdev));
[~,~,CIdev] = ttest(ERFdev);
[~,~,CIstd] = ttest(ERFstd);
[~,~,CIdif] = ttest(ERFdif);


% Deviant response

myfigure2
plot(t,nestedmean(ERFdev,IDdev),'k','linewidth',2)
patch([t fliplr(t)],[CIdev(2,:) fliplr(CIdev(1,:))],'k','facealpha',transp,'EdgeColor','none')

xticks([0 1000 2000 3000])
xlabel('Time (ms)')
ylabel('Normalized amplitude (AUs)')

plot(zeros(1,100),linspace(-100,100,100),'k--')
legend({'Average ERF','95% CI','Tone'},'autoupdate','off','fontsize',16)
legend box off

plot(ones(1,100).*600,linspace(0,200,100),'k--')
plot(ones(1,100).*1200,linspace(0,200,100),'k--')
plot(ones(1,100).*1800,linspace(0,200,100),'k--')
title('ERF (Deviant), Newborns','fontsize',18)
xlim([-200 3000])
ylim([0 7e-14])

makefighandsome
print('-dpng','./Figures/GrandAverageNeonatalDeviantERF.png')
print('-dsvg','./Figures/GrandAverageNeonatalDeviantERF.svg')

% Standard response

myfigure2
plot(t,nestedmean(ERFstd,IDstd),'k','linewidth',2)
patch([t fliplr(t)],[CIstd(2,:) fliplr(CIstd(1,:))],'k','facealpha',transp,'EdgeColor','none')

xticks([0 1000 2000 3000])
xlabel('Time (ms)')
ylabel('Normalized amplitude (AUs)')

plot(zeros(1,100),linspace(-100,100,100),'k--')
legend({'Average ERF','95% CI','Tone'},'autoupdate','off','fontsize',16)
legend box off

plot(ones(1,100).*600,linspace(0,200,100),'k--')
plot(ones(1,100).*1200,linspace(0,200,100),'k--')
plot(ones(1,100).*1800,linspace(0,200,100),'k--')
title('ERF (Standard), Newborns','fontsize',18)
xlim([-200 3000])
ylim([0 7e-14])

makefighandsome
print('-dpng','./Figures/GrandAverageNeonatalStandardERF.png')
print('-dsvg','./Figures/GrandAverageNeonatalStandardERF.svg')

% Dev - Std response

myfigure2
plot(t,nestedmean(ERFdif),'k','linewidth',2)
patch([t fliplr(t)],[CIdif(2,:) fliplr(CIdif(1,:))],'k','facealpha',transp,'EdgeColor','none')

xticks([0 1000 2000 3000])
xlabel('Time (ms)')
ylabel('Normalized amplitude (AUs)')

plot(zeros(1,100),linspace(-100,100,100),'k--')
legend({'Average ERF','95% CI','Tone'},'autoupdate','off','fontsize',16)
legend box off

plot(ones(1,100).*600,linspace(-100,100,100),'k--')
plot(ones(1,100).*1200,linspace(-100,100,100),'k--')
plot(ones(1,100).*1800,linspace(-100,100,100),'k--')
title('ERF (Deviant - Standard), Newborns','fontsize',18)
xlim([-200 3000])
ylim([-3e-14 3e-14])

makefighandsome
print('-dpng','./Figures/GrandAverageNeonatalDifferenceERF.png')
print('-dsvg','./Figures/GrandAverageNeonatalDifferenceERF.svg')

save Newborn_ERFs ERFdev ERFstd ERFdif IDdev IDstd