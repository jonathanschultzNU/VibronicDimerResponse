function MDplot(xdata, ydata, zdata,t2val)

%	Author: Jonathan D. Schultz
%	Email: jonathanschultz2022@u.northwestern.edu
%	Last revision date: February 1st, 2021
%
%	Copyright: Jonathan D. Schultz, 2021

%   Please see readme file for information about this package

cmap = cmap2d(10);

% Create contour
contourf(xdata,ydata,zdata,'LineWidth',0.25,'LevelStep',0.1);
colormap(cmap)
hold on

% Create line
line('XData',[min(ydata) max(ydata)],'YData',[min(ydata) max(ydata)],'Linewidth',1.5);

% Create ylabel
ylabel('\omega_{3} (cm^{-1})');
% ylabel('\lambda_{3} (nm)');

% Create xlabel
xlabel('\omega_{1} (cm^{-1})');
% xlabel('\lambda_{1} (nm)');

% Create title
title(t2val);

box on
set(gca,'BoxStyle','full','CLim',[-1 1],'DataAspectRatio',[1 1 1],...
    'FontSize',18,'Layer','top');

xlim([min(xdata) max(xdata)])
ylim([min(ydata) max(ydata)])

