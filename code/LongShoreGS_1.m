%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];
locsY = {'SL', 'S', 'L2' 'L3.5', 'L4', 'T', 'L6'};
locsYnew = {'Lgn', 'L0', 'L2', 'L3.5', 'L4', 'Tmb', 'L6'};


%% GS_L2

% Initialise table
dataPath = [folderPath 'grainsizes' filesep 'GS_20210920.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_L2 = readtable(dataPath, opts);
GS_L2 = GS_L2(:, 1:15);
GS_L2 = GS_L2(startsWith(GS_L2.Sample_Identity, ["L2C3","L2C5"]), :);

sampleDate = ['20210928'; '20210930'; '20211001'; '20211006';...
    '20211007'; '20211008'; '20211011'; '20211013'; '20211015'];

% Load remaining L2 sediment data
for i = 1:length(sampleDate)
    dataPath = [folderPath 'grainsizes' filesep 'GS_' sampleDate(i,:) '.csv'];
    opts = detectImportOptions(dataPath);
    opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
    GS_L2_temp = readtable(dataPath, opts);
    GS_L2_temp = GS_L2_temp(:, 1:15);

    GS_L2 = [GS_L2; GS_L2_temp];
end

GS_L2 = GS_L2(startsWith(GS_L2.Sample_Identity, ["L2C3_","L2C5W","L2C5_"]), :);
GS_L2.Sample_Number = nan(height(GS_L2), 1); % add dummy column

% Loop through each row of the table
for i = 1:height(GS_L2)
    % Get the Sample_Identity value for the current row
    sampleIdentity = GS_L2.Sample_Identity{i};
    
    % Check if Sample_Identity starts with 'L2C3'
    if startsWith(sampleIdentity, 'L2C3')
        GS_L2.Sample_Number(i) = 1;
    % Check if Sample_Identity starts with 'L2C5'
    elseif startsWith(sampleIdentity, 'L2C5')
        GS_L2.Sample_Number(i) = 3;
    end
end

GS_L2.Sample_Identity = regexprep(GS_L2.Sample_Identity, ["_.*", "W"], '');

% Group the table by 'Sample_Identity' and 'Date_ddMMyyyy' and calculate the mean of other columns
groupedTable = groupsummary(GS_L2, {'Sample_Identity', 'Date_ddMMyyyy'}, 'mean');
groupedTable.GroupCount = []; % remove extra column

% Set column names of the second table to be the same as the first table
groupedTable.Properties.VariableNames = GS_L2.Properties.VariableNames;

GS_L2 = groupedTable;
GS_L2.Sample_Identity = regexprep(GS_L2.Sample_Identity, "C.*", '');

clearvars GS_L2_temp sampleIdentity sampleDate groupedTable dataPath


%% GS_20211008/9

% [R10; R08; R07; R05; R02; R01]
% [L6; T; L4; L3.5; L2; S]

dataPath{1} = [folderPath 'grainsizes' filesep 'GS_20211008.csv'];
dataPath{2} = [folderPath 'grainsizes' filesep 'GS_20211009.csv'];

opts = detectImportOptions(dataPath{1});
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20211008 = readtable(dataPath{1}, opts);
GS_20211009 = readtable(dataPath{2}, opts);
GS_20211008 = [GS_20211009; GS_20211008];
GS_20211008 = GS_20211008(:, 1:15);

GS_20211008 = GS_20211008(startsWith(GS_20211008.Sample_Identity,...
    ["R01E","R02E","R05E","R07E","R08E","R10E"]), :);
GS_20211008.Sample_Number = repelem(1, 6)';

GS_20211008.Sample_Identity = {'S', 'L2', 'L3.5', 'L4', 'T', 'L6'}';

clearvars GS_20211009 dataPath


%% GS_20210921

% Load sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210921.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_20210921 = readtable(dataPath, opts);

% Prepare table
GS_20210921 = GS_20210921(~startsWith(GS_20210921.Sample_Identity, 'L2C'), :);
GS_20210921.Sample_Number = [1; 3; 1; 2; 3; 1; 3; 1; 3; 1; 3; 1; 3];
GS_20210921.Sample_Identity = regexprep(GS_20210921.Sample_Identity, '_.*', '');

% Visualisation
f1 = figureRH;
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(7,[5 5])
h1 = heatmap(GS_20210921, 'Sample_Number', 'Sample_Identity', 'ColorVariable','Mean_mu');
colormap(h1, crameri('lajolla'))
clim([0, 2000])

h1.Title = [];
h1.XDisplayLabels = {'0 m','-0.5 m','-0.75 m'};
h1.YDisplayData = locsY;
h1.YDisplayLabels = locsYnew;
h1.XLabel = '';
h1.YLabel = '';
h1.FontSize = fontsize;
h1.CellLabelFormat = '%.0f';
h1.ColorbarVisible = 'off';
h1.GridVisible = 'off';
h1.MissingDataColor = 'w';
h1.MissingDataLabel = 'no data';

nexttile(6)
text(0, .5, '2021-09-21', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
axis off

% heatdata = h1.ColorDisplayData(any(~isnan(h1.ColorDisplayData), 2), :);
heatdata = flipud(h1.ColorDisplayData);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
nexttile(1,[1 5])
errorbar(1.5:3.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 4])
% ylim([0 1000])
ylabel('µm')
xticks([])

meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');
nexttile(12,[5 1])
errorbar(meanT, linspace(1.5, 6.5, 7), stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 7])
% xlim([0 1000])
xlabel('µm')
yticks([])


%% GS_20210928

% Load sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210928.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_20210928 = readtable(dataPath, opts);

% Prepare table
GS_20210928 = GS_20210928(~startsWith(GS_20210928.Sample_Identity, 'L2C'), :);
GS_20210928.Sample_Number = [1; 3; 1; 2; 3; 1; 3; 1; 3; 1; 3; 1; 3];
GS_20210928.Sample_Identity = regexprep(GS_20210928.Sample_Identity, '_.*', '');

% Visualisation
f2 = figureRH;
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(7,[5 5])
h2 = heatmap(GS_20210928, 'Sample_Number', 'Sample_Identity', 'ColorVariable','Mean_mu');
colormap(h2, crameri('lajolla'))
clim([0, 2000])

h2.Title = [];
h2.XDisplayLabels = {'0 m','-0.5 m','-0.75 m'};
h2.YDisplayData = locsY;
h2.YDisplayLabels = locsYnew;
h2.XLabel = '';
h2.YLabel = '';
h2.FontSize = fontsize;
h2.CellLabelFormat = '%.0f';
h2.ColorbarVisible = 'off';
h2.GridVisible = 'off';
h2.MissingDataColor = 'w';
h2.MissingDataLabel = 'no data';

nexttile(6)
text(0, .5, '2021-09-28', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
axis off

heatdata = flipud(h2.ColorDisplayData);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
nexttile(1,[1 5])
errorbar(1.5:3.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 4])
% ylim([0 1000])
ylabel('µm')
xticks([])

meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');
nexttile(12,[5 1])
errorbar(meanT, linspace(1.5, 6.5, 7), stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 7])
% xlim([0 1000])
xlabel('µm')
yticks([])


%% GS_20211003

% Load sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20211003.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_20211003 = readtable(dataPath, opts);

% Prepare table
GS_20211003.Time_hhmmss = [];
GS_20211003(strcmp(GS_20211003{:,1}, 'TB'), :) = [];  % not using the transverse bar sample
GS_20211003.Sample_Number = [1; 3; 1; 2; 3; 1; 3; 1; 3];
GS_20211003.Sample_Identity = regexprep(GS_20211003.Sample_Identity, '_.*', '');

% Visualisation
f3 = figureRH;
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(7,[5 5])
h3 = heatmap(GS_20211003, 'Sample_Number', 'Sample_Identity', 'ColorVariable','Mean_mu');
colormap(h3, crameri('lajolla'))
clim([0, 2000])

h3.Title = [];
h3.XDisplayLabels = {'0 m','-0.5 m','-0.75 m'};
h3.YDisplayData = locsY;
h3.YDisplayLabels = locsYnew;
h3.XLabel = '';
h3.YLabel = '';
h3.FontSize = fontsize;
h3.CellLabelFormat = '%.0f';
h3.ColorbarVisible = 'off';
h3.GridVisible = 'off';
h3.MissingDataColor = 'w';
h3.MissingDataLabel = 'no data';

nexttile(6)
text(0, .5, '2021-10-03', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
axis off

heatdata = flipud(h3.ColorDisplayData);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
nexttile(1,[1 5])
errorbar(1.5:3.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 4])
% ylim([0 1000])
ylabel('µm')
xticks([])

meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');
nexttile(12,[5 1])
errorbar(meanT, linspace(1.5, 6.5, 7), stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 7])
% xlim([0 1000])
xlabel('µm')
yticks([])


%% Longshore complete NAP -x m

% Sample_Number == 1 --> NAP + 0.00 m
% Sample_Number == 2 --> NAP - 0.50 m
% Sample_Number == 3 --> NAP - 0.75 m

% Merge tables
Longshore = [GS_20210921(:, [1:15, 31]); GS_20210928(:, [1:15, 31]);...
    GS_20211003(:, [1:15, 31]); GS_20211008; GS_L2];
LongshoreN = Longshore(Longshore.Sample_Number == 1, :);

% Find rows with dates equal to OLD_DATE and replace them with NEW_DATE
indices = LongshoreN.Date_ddMMyyyy == datetime('09/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN.Date_ddMMyyyy(indices) = datetime('08/10/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN.Date_ddMMyyyy == datetime('20/09/2021', 'Format', 'dd/MM/yyyy');
LongshoreN.Date_ddMMyyyy(indices) = datetime('21/09/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN.Date_ddMMyyyy == datetime('01/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN.Date_ddMMyyyy(indices) = datetime('03/10/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN.Date_ddMMyyyy == datetime('07/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN.Date_ddMMyyyy(indices) = datetime('08/10/2021', 'Format', 'dd/MM/yyyy');


%% Visualisation: longshore mean (detailed)
f4 = figureRH;
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(7,[5 5])
h4 = heatmap(LongshoreN, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Mean_mu');
colormap(h4, crameri('lajolla'))
% clim([0, 2000])

h4.Title = [];
h4.XDisplayLabels = datetime(h4.XDisplayLabels, 'Format','dd/MM');
h4.YDisplayData = locsY(2:end);
h4.YDisplayLabels = locsYnew(2:end);
h4.XLabel = '';
h4.YLabel = '';
h4.FontSize = fontsize;
h4.CellLabelFormat = '%.0f';
h4.ColorbarVisible = 'off';
h4.GridVisible = 'off';
h4.MissingDataColor = 'w';
h4.MissingDataLabel = 'no data';

nexttile(6)
% text(0, .5, 'NAP -0.00 m', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(0, .5, ['NAP ', num2str(mean(LongshoreN.zNAP_m), '%.1f'), ' m'], 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
axis off

heatdata = flipud(h4.ColorDisplayData);
meanT = mean(heatdata, 1, 'omitmissing');
stdT = std(heatdata, 0, 1, 'omitmissing');
nexttile(1,[1 5])
errorbar(1.5:9.5, meanT, stdT, '-ok', 'LineWidth',3)
yline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 10])
% ylim([0 1000])
ylabel('µm')
xticks([])

meanY = mean(heatdata, 2, 'omitmissing');
stdY = std(heatdata, 0, 2, 'omitmissing');
nexttile(12,[5 1])
errorbar(meanY, 1.5:6.5, stdY, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanY, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 7])
% xlim([0 1000])
xlabel('µm')
yticks([])

