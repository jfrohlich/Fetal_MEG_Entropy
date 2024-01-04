
% edited 06.21.2019
% Joel Frohlich
% 05.23.19
% This is a wrapper function that makes it easier to use the function
% patch() to plot confidence interval highlights

% First input param: matrix of data that needs 95% CI plotting (this
% function will take the variance across rows)

% Second input param: RGB color vector (random color used if not specified)

% This input param: transparency of CI highlight (default = 30%) 

function[] = CI_highlight(X,Y,colorx,colory,transp)

% if RGB color vector not specified, do random color

if nargin < 4
    transp = 0.3;
end
if nargin < 3
    colorx = rand(1,3);
    colory = rand(1,3);
end

% get confidence interval 
[~,~,CIx] = ttest(X);
[~,~,CIy] = ttest(Y);

figure('color','w')
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.3 0.8]);
% draw patch
plot(mean(X),'color',colorx,'linewidth',2)
hold on
patch([1:size(X,2) size(X,2):-1:1],[CIx(2,:) fliplr(CIx(1,:))],colorx,'facealpha', ...
    transp,'EdgeColor','none')
plot(mean(Y),'color',colory,'linewidth',2)
patch([1:size(Y,2) size(Y,2):-1:1],[CIy(2,:) fliplr(CIy(1,:))],colory,'facealpha', ...
    transp,'EdgeColor','none')

makefighandsome
    
end