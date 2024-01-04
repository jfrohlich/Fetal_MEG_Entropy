% Joel Frohlich 
% University of Tuebingen
% Generates the ERF response template for fetuses n weeks and older

% Okay, what is the script doing in this repo? After all, we didn't use its
% results in the manuscript. 
%
% Once upon a time, we had the idea to try looking at the entropy of
% signals after computing a grand-averaged ERF from the older
% fetuses and then subtracting this template from each recording. (Why 
% older fetuses? Because they show a more reliable P300 like response, see
% Moser et al. 2021). The idea was to use this approach to look at the
% background signal rather than the evoked response. It turns out this was
% a bad idea and we abandoned it. But this code is still included in the
% code because its output got grandfathered into later code as input.
%
% If nothing else, you may still find this script useful for plotting
% grand-averaged ERFs across recordings ... 

% SUMMARY

%     Loading Data:
%         Close all figures, clear variables, and load the MEG data from 'Data_traces.mat'.
% 
%     Data Filtering:
%         Initialize empty arrays ERFdev, ERFstd, and ID.
%         Set an age threshold agethresh to 32 weeks (can be adjusted as
%         desired or necessary).
%         Iterate through the rows of 'Datatable' and select rows where the gestational age (GA) is greater than or equal to the age threshold.
%         Populate ERFdev, ERFstd, and ID arrays with corresponding data from 'Datatable'.
% 
%     Calculate Differences:
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
%         Save ERFdev, ERFstd, ERFdif, and ID in a file named 'grand_average_ERFs.mat'.

close all
clearvars
load Data_traces.mat

ERFdev = [];
ERFstd = [];
ID = [];
agethresh = 32;

for irow = 1:size(Datatable,1)
    if Datatable.GA(irow) >= agethresh
        ERFdev = [ERFdev; Datatable.global_deviants_normalized{irow}];
        ERFstd = [ERFstd; Datatable.global_standards_normalized{irow}];
        ID     = [ID; Datatable.ID(irow)];
    end
end

ERFdif = ERFdev-ERFstd;

transp = 0.3;
t = linspace(-200,3000,length(ERFdev));
[~,~,CIdev] = ttest(ERFdev);
[~,~,CIstd] = ttest(ERFstd);
[~,~,CIdif] = ttest(ERFdif);


% Deviant response

myfigure2
plot(t,nestedmean(ERFdev,ID),'k','linewidth',2)
patch([t fliplr(t)],[CIdev(2,:) fliplr(CIdev(1,:))],'k','facealpha',transp,'EdgeColor','none')

xticks([0 1000 2000 3000])
xlabel('Time (ms)')
ylabel('Normalized amplitude (AUs)')

plot(zeros(1,100),linspace(0,200,100),'k--')
legend({'Average ERF','95% CI','Tone'},'autoupdate','off','fontsize',16)
legend box off

plot(ones(1,100).*600,linspace(0,200,100),'k--')
plot(ones(1,100).*1200,linspace(0,200,100),'k--')
plot(ones(1,100).*1800,linspace(0,200,100),'k--')
title(sprintf('ERF (Deviant), GA >= %i months',agethresh),'fontsize',18)
xlim([-200 3000])
ylim([40 160])

makefighandsome
print('-dpng','./Figures/GrandAverageDeviantERF35AndUp.png')
print('-dsvg','./Figures/GrandAverageDeviantERF35AndUp.svg')

% Standard response

myfigure2
plot(t,nestedmean(ERFstd,ID),'k','linewidth',2)
patch([t fliplr(t)],[CIstd(2,:) fliplr(CIstd(1,:))],'k','facealpha',transp,'EdgeColor','none')

xticks([0 1000 2000 3000])
xlabel('Time (ms)')
ylabel('Normalized amplitude (AUs)')

plot(zeros(1,100),linspace(0,200,100),'k--')
legend({'Average ERF','95% CI','Tone'},'autoupdate','off','fontsize',16)
legend box off

plot(ones(1,100).*600,linspace(0,200,100),'k--')
plot(ones(1,100).*1200,linspace(0,200,100),'k--')
plot(ones(1,100).*1800,linspace(0,200,100),'k--')
title(sprintf('ERF (Standard), GA >= %i months',agethresh),'fontsize',18)
xlim([-200 3000])
ylim([40 160])

makefighandsome
print('-dpng','./Figures/GrandAverageStandardERF32AndUp.png')
print('-dsvg','./Figures/GrandAverageStandardERF32AndUp.svg')

% Dev - Std response

myfigure2
plot(t,nestedmean(ERFdif,ID),'k','linewidth',2)
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
title(sprintf('ERF (Deviant - Standard), GA >= %i months',agethresh),'fontsize',18)
xlim([-200 3000])
ylim([-50 50])

makefighandsome
print('-dpng','./Figures/GrandAverageDifferenceERF32AndUp.png')
print('-dsvg','./Figures/GrandAverageDifferenceERF32AndUp.svg')

save grand_average_ERFs ERFdev ERFstd ERFdif ID
