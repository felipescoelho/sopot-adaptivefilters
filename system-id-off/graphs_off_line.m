% Figure generator for the System Identification
%
% Luiz Felipe da S. Coelho - lfscoelho@ieee.org
% mar 22, 2021
%

clc
clear
close all

%% Load everything
load('SystemID_off.mat')
load('meta.mat')
% Simulation info:
M = 32;  % For the last run
lambda = 0.97;
delta = 0.42;
mu_LMS = 0.1;
mu_NLMS = 0.098;
gamma = 1e-12;
K = 1000;
N = 4;
normH = 1;
SNR = [-3 3 10 15];

%% Norm w
y_label = '$\left|\hspace{-.05cm}\left|\mathbf{w} (k)\right|\hspace{-.05cm}\right|_2$';
x_label = 'Iterations, $k$';
plot1(1:K, norm_w_RLS_avg(2:end, :), norm_w_RLS_SPT_avg(2:end, :, :), K,...
    'norm_w_RLS', x_label, y_label)
plot1(1:K, norm_w_LMS_avg(2:end, :), norm_w_LMS_SPT_avg(2:end, :, :), K,...
    'norm_w_LMS', x_label, y_label)
plot1(1:K, norm_w_NLMS_avg(2:end, :), norm_w_NLMS_SPT_avg(2:end, :, :),...
    K, 'norm_w_NLMS', x_label, y_label)

%% MSE
y_label = 'MSE';
x_label = 'Iterations, $k$';
plot2(1:K, 10*log10(MSE_RLS_avg), 10*log10(MSE_RLS_SPT_avg), 'MSE_RLS',...
    x_label, y_label)
plot2(1:K, 10*log10(MSE_LMS_avg), 10*log10(MSE_LMS_SPT_avg), 'MSE_LMS',...
    x_label, y_label)
plot2(1:K, 10*log10(MSE_NLMS_avg), 10*log10(MSE_NLMS_SPT_avg),...
    'MSE_NLMS', x_label, y_label)

%% MSD
y_label = 'MSD';
x_label = 'Iterations, $k$';
plot2(1:K, 10*log10(MSD_RLS_avg(2:end, :)),...
    10*log10(MSD_RLS_SPT_avg(2:end, :, :)), 'MSD_RLS', x_label, y_label)
plot2(1:K, 10*log10(MSD_LMS_avg(2:end, :)),...
    10*log10(MSD_LMS_SPT_avg(2:end, :, :)), 'MSD_LMS', x_label, y_label)
plot2(1:K, 10*log10(MSD_NLMS_avg(2:end, :)),...
    10*log10(MSD_NLMS_SPT_avg(2:end, :, :)), 'MSD_NLMS', x_label, y_label)


%% Functions
function plot1(x, y, y_SPT, K, name, x_label, y_label)
% Definitions

% Colors
flame_scarlet = [205 33 42]./255;
saffron = [255 165 0]./255;
classic_blue = [15 76 129]./255;
biscay_green = [86 198 169]./255;

% IEEE one column width = 8.99 com
% IEEE two columns width = 18.2 cm
gold = (1 + sqrt(5))./2;  % golden ration, for the relation height width
WIDTH = 3.5*8.99;  % figure width
HEIGHT = WIDTH./gold;  % figure height
fontSize = 25;  % font size
fontName = 'Times New Roman';

fig = figure;
fig.Units = 'centimeters';
fig.Position = [2 2 WIDTH HEIGHT];
fig.Color = 'w';
fig.Name = name;
xDist = 3;
yDist = 2.5;
plotHeight = HEIGHT/2 - yDist + 0.3;
plotWidth = WIDTH - xDist - 1.1;
lineWidth = 2;

color1 = flame_scarlet;
color2 = saffron;
color3 = classic_blue;
color4 = biscay_green;

% generate axes
ax0 = axes('units', 'centimeters', 'position', [xDist-1, yDist,...
    WIDTH-xDist, HEIGHT-yDist]);
ax1 = axes('units', 'centimeters', 'position', [xDist, HEIGHT/2 + yDist/2,...
    plotWidth, plotHeight]);  % upper axis
ax2 = axes('units', 'centimeters', 'position', [xDist, yDist, plotWidth,...
    plotHeight]);  % lower axis

subplot(ax0)  % plot axis 0, just for the cetered y label
set(ax0, 'XColor', 'w')
set(ax0, 'YColor', 'w')
ylabel(y_label, 'interpreter', 'latex', 'FontSize', fontSize, 'Color',...
    'k', 'FontName', fontName)
% ------------------------------------------------------------------------
% upper axis
subplot(ax1)
plot(x, ones(1, K), '-.k', 'linewidth', lineWidth), hold on
plot(x, y_SPT(:, end-12, 1), 'linewidth', lineWidth, 'Color', color1)
plot(x, y_SPT(:, end-8, 1), 'linewidth', lineWidth, 'Color', color2)
plot(x, y_SPT(:, end, 1), 'linewidth', lineWidth, 'Color', color3)
plot(x, y(:, 1), '--', 'linewidth', lineWidth, 'Color', color4)
hold off, grid on
set(ax1, 'FontSize', fontSize)
set(ax1, 'XTickLabel', [])
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'linewidth', 2)
set(ax1, 'box', 'off')
ylim([0 1.5])
ylabel('(a)', 'interpreter', 'latex', 'FontSize', fontSize, 'rotation',...
    0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'Color', 'k')
% ------------------------------------------------------------------------
% lower axis
subplot(ax2)
plot(x, ones(1, K), '-.k', 'linewidth', lineWidth), hold on
plot(x, y_SPT(:, end-12, 3), 'linewidth', lineWidth, 'Color', color1)
plot(x, y_SPT(:, end-8, 3), 'linewidth', lineWidth, 'Color', color2)
plot(x, y_SPT(:, end, 3), 'linewidth', lineWidth, 'Color', color3)
plot(x, y(:, 3), '--', 'linewidth', lineWidth, 'Color', color4)
ylim([0 1.5])
rearange = get(ax2, 'Children');
set(ax2, 'Children', rearange([1 5 2 3 4]))
lgnd = legend('$M/N=5$', '$M/N=6$', '$M/N=8$',...
    'Unknown System, $H$', 'Full Precision');
lgnd.Units = 'centimeters';
lgnd.Position(1) = 8.4;
lgnd.Position(2) = 3.4;
lgnd.NumColumns = 2;
lgnd.Interpreter = 'latex';
lgnd.FontSize = fontSize;
lgnd.FontName = fontName;
set(ax2, 'FontSize', fontSize)
set(ax2, 'TickLabelInterpreter', 'latex')
set(ax2, 'linewidth', 2)
set(ax2, 'box', 'off')
ax2.XAxis.Exponent = 3;
ylabel('(b)', 'interpreter', 'latex', 'FontSize', fontSize, 'rotation',...
    0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'FontName', fontName, 'Color', 'k')
xlabel(x_label, 'interpreter', 'latex', 'FontSize', fontSize,...
    'FontName', fontName, 'Color', 'k')
hold off, grid on


saveas(fig, name, 'epsc')

end

function plot2(x, y, y_SPT, name, x_label, y_label)
% Definitions

% Colors
flame_scarlet = [205 33 42]./255;
saffron = [255 165 0]./255;
classic_blue = [15 76 129]./255;
biscay_green = [86 198 169]./255;

% IEEE one column width = 8.99 com
% IEEE two columns width = 18.2 cm
gold = (1 + sqrt(5))./2;  % golden ration, for the relation height width
WIDTH = 3.5*8.99;  % figure width
HEIGHT = WIDTH./gold;  % figure height
fontSize = 25;  % font size
fontName = 'Times New Roman';

fig = figure;
fig.Units = 'centimeters';
fig.Position = [2 2 WIDTH HEIGHT];
fig.Color = 'w';
fig.Name = name;
xDist = 3;
yDist = 2.5;
plotHeight = HEIGHT/2 - yDist + 0.3;
plotWidth = WIDTH - xDist - 1.1;
lineWidth = 2;

color1 = flame_scarlet;
color2 = saffron;
color3 = classic_blue;
color4 = biscay_green;

% generate axes
ax0 = axes('units', 'centimeters', 'position', [xDist-1, yDist,...
    WIDTH-xDist, HEIGHT-yDist]);
ax1 = axes('units', 'centimeters', 'position', [xDist, HEIGHT/2 + yDist/2,...
    plotWidth, plotHeight]);  % upper axis
ax2 = axes('units', 'centimeters', 'position', [xDist, yDist, plotWidth,...
    plotHeight]);  % lower axis

subplot(ax0)  % plot axis 0, just for the cetered y label
set(ax0, 'XColor', 'w')
set(ax0, 'YColor', 'w')
ylabel(y_label, 'interpreter', 'latex', 'FontSize', fontSize, 'Color',...
    'k', 'FontName', fontName)
% ------------------------------------------------------------------------
% upper axis
subplot(ax1)
plot(x, y_SPT(:, end-12, 1), 'linewidth', lineWidth, 'Color', color1), hold on
plot(x, y_SPT(:, end-8, 1), 'linewidth', lineWidth, 'Color', color2)
plot(x, y_SPT(:, end, 1), 'linewidth', lineWidth, 'Color', color3)
plot(x, y(:, 1), '--', 'linewidth', lineWidth, 'Color', color4)
hold off, grid on
set(ax1, 'FontSize', fontSize)
set(ax1, 'XTickLabel', [])
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'linewidth', 2)
set(ax1, 'box', 'off')
% ylim([0 1.5])
ylabel('(a)', 'interpreter', 'latex', 'FontSize', fontSize, 'rotation',...
    0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'Color', 'k')
% ------------------------------------------------------------------------
% lower axis
subplot(ax2)
plot(x, y_SPT(:, end-12, 3), 'linewidth', lineWidth, 'Color', color1), hold on
plot(x, y_SPT(:, end-8, 3), 'linewidth', lineWidth, 'Color', color2)
plot(x, y_SPT(:, end, 3), 'linewidth', lineWidth, 'Color', color3)
plot(x, y(:, 3), '--', 'linewidth', lineWidth, 'Color', color4)
% ylim([0 1.5])
rearange = get(ax2, 'Children');
set(ax2, 'Children', rearange([1 2 3 4]))
lgnd = legend('$M/N=5$', '$M/N=6$', '$M/N=8$', 'Full Precision');
lgnd.Units = 'centimeters';
lgnd.Position(1) = 8.4;
lgnd.Position(2) = 3.4;
lgnd.NumColumns = 2;
lgnd.Interpreter = 'latex';
lgnd.FontSize = fontSize;
lgnd.FontName = fontName;
set(ax2, 'FontSize', fontSize)
set(ax2, 'TickLabelInterpreter', 'latex')
set(ax2, 'linewidth', 2)
set(ax2, 'box', 'off')
ax2.XAxis.Exponent = 3;
ylabel('(b)', 'interpreter', 'latex', 'FontSize', fontSize, 'rotation',...
    0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'FontName', fontName, 'Color', 'k')
xlabel(x_label, 'interpreter', 'latex', 'FontSize', fontSize,...
    'FontName', fontName, 'Color', 'k')
hold off, grid on


saveas(fig, name, 'epsc')

end

function plot3(x, y, y_SPT, name, x_label, y_label)
% Definitions

% Colors
flame_scarlet = [205 33 42]./255;
saffron = [255 165 0]./255;
classic_blue = [15 76 129]./255;
biscay_green = [86 198 169]./255;

% IEEE one column width = 8.99 com
% IEEE two columns width = 18.2 cm
gold = (1 + sqrt(5))./2;  % golden ration, for the relation height width
WIDTH = 3.5*8.99;  % figure width
HEIGHT = WIDTH./gold;  % figure height
fontSize = 25;  % font size
fontName = 'Times New Roman';

fig = figure;
fig.Units = 'centimeters';
fig.Position = [2 2 WIDTH HEIGHT];
fig.Color = 'w';
fig.Name = name;
xDist = 3;
yDist = 2.5;
plotHeight = HEIGHT/2 - yDist + 0.3;
plotWidth = WIDTH - xDist - 1.1;
lineWidth = 2;

color1 = flame_scarlet;
color2 = saffron;
color3 = classic_blue;
color4 = biscay_green;

% generate axes
ax0 = axes('units', 'centimeters', 'position', [xDist-1, yDist,...
    WIDTH-xDist, HEIGHT-yDist]);
ax1 = axes('units', 'centimeters', 'position', [xDist, HEIGHT/2 + yDist/2,...
    plotWidth, plotHeight]);  % upper axis
ax2 = axes('units', 'centimeters', 'position', [xDist, yDist, plotWidth,...
    plotHeight]);  % lower axis
axZoom1 = axes('Units', 'centimeters', 'position', [WIDTH-zoomPlotWidth-.35,...
    HEIGHT-zoomPlotHeight-.1, zoomPlotWidth, zoomPlotHeight]);
axZoom2 = axes('Units', 'centimeters', 'position', [WIDTH-zoomPlotWidth-.35,...
    HEIGHT/2-zoomPlotHeight+1.1, zoomPlotWidth, zoomPlotHeight]);

subplot(ax0)  % plot axis 0, just for the cetered y label
set(ax0, 'XColor', 'w')
set(ax0, 'YColor', 'w')
ylabel(y_label, 'interpreter', 'latex', 'FontSize', fontSize, 'Color',...
    'k', 'FontName', fontName)
% ------------------------------------------------------------------------
% upper axis
subplot(ax1)
plot(x, y_SPT(:, end-12, 1), 'linewidth', lineWidth, 'Color', color1), hold on
plot(x, y_SPT(:, end-8, 1), 'linewidth', lineWidth, 'Color', color2)
plot(x, y_SPT(:, end, 1), 'linewidth', lineWidth, 'Color', color3)
plot(x, y(:, 1), '--', 'linewidth', lineWidth, 'Color', color4)
hold off, grid on
set(ax1, 'FontSize', fontSize)
set(ax1, 'XTickLabel', [])
set(ax1, 'TickLabelInterpreter', 'latex')
set(ax1, 'linewidth', 2)
set(ax1, 'box', 'off')
% ylim([0 1.5])
ylabel('(a)', 'interpreter', 'latex', 'FontSize', fontSize, 'rotation',...
    0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'Color', 'k')
% ------------------------------------------------------------------------
% upper axis zoom
subplot(axZoom1)
plot(zoomStart:zoomEnd, 10*log10(y_SPT(zoomStart+1:zoomEnd+1, end-12, 1)), 'LineWidth', lineWidth, 'Color', color1), hold on
plot(zoomStart:zoomEnd, 10*log10(y_SPT(zoomStart+1:zoomEnd+1, end-8, 1)), 'linewidth', lineWidth, 'Color', color2)
plot(zoomStart:zoomEnd, 10*log10(y_SPT(zoomStart+1:zoomEnd+1, end, 1)), 'linewidth', lineWidth, 'Color', color3)
plot(zoomStart:zoomEnd, 10*log10(y(zoomStart+1:zoomEnd+1, 1)), '--', 'linewidth', lineWidth, 'Color', color4), hold off, grid on
xlim([zoomStart zoomEnd])
xticks([leftTick centerTick rightTick])
ylim([-11.5 -7.3])
yticks([-11.5 -10 -9])
set(axZoom1, 'TickLabelInterpreter', 'latex')
set(axZoom1, 'FontSize', fontSize-5)
set(axZoom1, 'FontName', fontName)
set(axZoom1, 'linewidth', 2)
set(axZoom1, 'XTickLabel', tickVector)
% ------------------------------------------------------------------------
% lower axis
subplot(ax2)
plot(x, y_SPT(:, end-12, 3), 'linewidth', lineWidth, 'Color', color1), hold on
plot(x, y_SPT(:, end-8, 3), 'linewidth', lineWidth, 'Color', color2)
plot(x, y_SPT(:, end, 3), 'linewidth', lineWidth, 'Color', color3)
plot(x, y(:, 3), '--', 'linewidth', lineWidth, 'Color', color4)
% ylim([0 1.5])
rearange = get(ax2, 'Children');
set(ax2, 'Children', rearange([1 2 3 4]))
lgnd = legend('$M/N=5$', '$M/N=6$', '$M/N=8$', 'Full Precision');
lgnd.Units = 'centimeters';
lgnd.Position(1) = 8.4;
lgnd.Position(2) = 3.4;
lgnd.NumColumns = 2;
lgnd.Interpreter = 'latex';
lgnd.FontSize = fontSize;
lgnd.FontName = fontName;
set(ax2, 'FontSize', fontSize)
set(ax2, 'TickLabelInterpreter', 'latex')
set(ax2, 'linewidth', 2)
set(ax2, 'box', 'off')
ax2.XAxis.Exponent = 3;
ylabel('(b)', 'interpreter', 'latex', 'FontSize', fontSize, 'rotation',...
    0, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right',...
    'FontName', fontName, 'Color', 'k')
xlabel(x_label, 'interpreter', 'latex', 'FontSize', fontSize,...
    'FontName', fontName, 'Color', 'k')
hold off, grid on
% ------------------------------------------------------------------------
% lower axis zoom
subplot(axZoom2)
plot(zoomStart:zoomEnd, 10*log10(y_SPT(zoomStart+1:zoomEnd+1,...
    end-5, 3)), 'LineWidth', lineWidth, 'Color', color1), hold on
plot(zoomStart:zoomEnd, 10*log10(av_MSD_sopot(zoomStart+1:zoomEnd+1,...
    end-3, 3)), 'linewidth', lineWidth, 'Color', color2)
plot(zoomStart:zoomEnd, 10*log10(av_MSD_sopot(zoomStart+1:zoomEnd+1,...
    end, 3)), 'linewidth', lineWidth, 'Color', color3)
plot(zoomStart:zoomEnd, 10*log10(av_MSD_fi(zoomStart+1:zoomEnd+1, end,...
    3)), '--', 'linewidth', lineWidth, 'Color', color4), hold off, grid on
xlim([zoomStart zoomEnd])
xticks([leftTick centerTick rightTick])
ylim([-25 -16])
yticks([-23 -19])
set(axZoom2, 'TickLabelInterpreter', 'latex')
set(axZoom2, 'FontSize', fontSize-5)
set(axZoom2, 'FontName', fontName)
set(axZoom2, 'linewidth', 2)
set(axZoom2, 'XTickLabel', tickVector)


saveas(fig, name, 'epsc')
end


% EoF
