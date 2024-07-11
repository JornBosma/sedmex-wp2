%% Initialisation
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;
% fontsize = 26; % ultra-wide screen

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];


%% Load sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20201202.csv'];
% dataPath = [folderPath 'grainsizes' filesep 'GS_20221026.csv'];

opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

% Option 1
GS_20201202 = readtable(dataPath, opts);
sieveSizes = [8000, 4000, 2800, 2000, 1700, 1400, 1180, 1000, 850, 710,...
    600, 500, 425, 355, 300, 250, 212, 180, 150, 125, 90, 75, 63, 45, 0];
massRetained = GS_20201202{23, 17:end};  % T08S01

% Option 2
% GS_20221026 = readtable(dataPath, opts);
% sieveSizes = [8000, 4000, 2000, 1000, 710, 500, 425, 355, 300, 250, 180, 125, 63, 0];
% massRetained = GS_20221026{1, 16:end};


%% Computations
cumulativeMass = cumsum(massRetained);
totalMass = cumulativeMass(end);

normalizedMass = cumulativeMass / totalMass;
normalizedMass = 1-normalizedMass;

frequency = [0, diff(cumulativeMass) / totalMass * 100];


%% Visualisation
f1 = figure('Position',[2095, 559, 874, 504]);

ax1 = axes;
plot(ax1, sieveSizes, normalizedMass*100, '-', 'Color',cbf.redpurp,...
    'LineWidth',4); hold on
scatter(ax1, sieveSizes, normalizedMass*100, 150, 'MarkerEdgeColor',cbf.blue,...
    'MarkerFaceColor',cbf.blue, 'Marker','square')

D90 = interp1(normalizedMass(1:end-2)*100, sieveSizes(1:end-2), 90, 'pchip');
D50 = interp1(normalizedMass(1:end-2)*100, sieveSizes(1:end-2), 50, 'pchip');
D10 = interp1(normalizedMass(1:end-2)*100, sieveSizes(1:end-2), 10, 'pchip');
text(280, 85, ['D_{90} = ',num2str(D90, 4),' μm'], 'FontSize',fontsize*.9)
text(120, 45, ['D_{50} = ',num2str(D50, 3),' μm'], 'FontSize',fontsize*.9)
text(80, 15, ['D_{10} = ',num2str(D10, 3),' μm'], 'FontSize',fontsize*.9)

line([D90, D90], [20, 90], 'Color','k', 'LineStyle','--', 'LineWidth',2)
line([D50, D50], [0, 50], 'Color','k', 'LineStyle','--', 'LineWidth',2)
line([D10, D10], [0, 10], 'Color','k', 'LineStyle','--', 'LineWidth',2)

line([sieveSizes(end-1), D90], [90, 90], 'Color','k', 'LineStyle','--', 'LineWidth',2)
line([sieveSizes(end-1), D50], [50, 50], 'Color','k', 'LineStyle','--', 'LineWidth',2)
line([sieveSizes(end-1), D10], [10, 10], 'Color','k', 'LineStyle','--', 'LineWidth',2)

ax1.XScale = 'log';
xlim(ax1, [min(sieveSizes) max(sieveSizes)])

xlabel('particle diameter (μm)') % tex
ylabel('cumulative mass retained (%)') % tex

grid on
grid minor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rectangle(ax1, 'Position', [500 2 7000 54], 'FaceColor','white', 'EdgeColor','none')
position = [0.5 0.25 0.32 0.32];  % [x y width height]
ax2 = axes('Position', position);

boundaries = [fliplr(sieveSizes)];
h = stairs(boundaries, [fliplr(frequency)], 'LineStyle','none'); hold on

x = repelem(h.XData(2:end), 2);
y = repelem(h.YData(1:end-1), 2);
x(end) = []; y(1) = [];
fill([x, fliplr(x)], [y, min(h.YData)*ones(size(y))], cbf.blue, 'LineWidth',1)

for n = 1:length(h.XData)
    line([h.XData(n), h.XData(n)], [0, h.YData(n)], 'Color','k', 'LineWidth',1)
end

ax2.XScale = 'log';
ax2.XAxisLocation = 'bottom';
ax2.YAxisLocation = 'right';
ax2.Color = [0.9 0.9 0.9];
ax2.FontSize = fontsize*.9;

% xlabel('particle diameter (μm)', 'FontSize',fontsize*.9)
ylabel('class weight (%)', 'FontSize',fontsize*.9)

axis padded


