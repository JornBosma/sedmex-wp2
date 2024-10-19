%% Initialisation
close all
clear
clc

[~, fontsize, cbf, PHZ, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep...
    'DataDescriptor' filesep 'grainsizes' filesep];

% List all files in the folder
fileList = dir(fullfile(folderPath, 'GS_*.csv'));

S = struct();

for n = 1:length(fileList)
    fileName = fileList(n).name;
    dataPath{1} = [folderPath filesep fileName];

    opts = detectImportOptions(dataPath{1});
    opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy'); 

    [~, fileNam, ~] = fileparts(fileName); % Extract the name without extension
    T = readtable(dataPath{1}, opts);
    S.(fileNam) = T;
end


%% Load sediment data
dataPath{2} = [folderPath 'GS_20211008.csv'];
dataPath{3} = [folderPath 'GS_20211009.csv'];
dataPath{4} = [folderPath 'GS_20221026.csv'];

opts = detectImportOptions(dataPath{2});
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20211008 = readtable(dataPath{2}, opts);
GS_20211009 = readtable(dataPath{3}, opts);
GS_2021 = [GS_20211009; GS_20211008];
GS_2022 = readtable(dataPath{4}, opts);


%% Sort table by track number and sample elevation
% Initialize an empty array to store the extracted numbers
Track_Number = zeros(size(GS_2021, 1), 1);

% Loop through each row in the 'name' column and extract the numbers
for row = 1:size(GS_2021, 1)
    name = GS_2021.Sample_Identity{row}; % Get the name from the 'name' column
    numbers = regexp(name, '\d+', 'match'); % Extract all numbers using regular expression
    if ~isempty(numbers)
        Track_Number(row) = str2double(numbers{1}); % Convert the first extracted number to a double
    end
end

% Convert the cell array of numbers to a column in your table
GS_2021.Track_Number = Track_Number;

GS_2021 = sortrows(GS_2021, 'Track_Number', 'ascend');

% Assign sample numbers by elevation
Sample_Number = repelem([1 2 3 4 5], 10, 1)';
Sample_Number = Sample_Number(:);
GS_2021.Sample_Number = Sample_Number;

% Remove all characters after the last number in each name
GS_2021.Sample_Identity = regexprep(GS_2021.Sample_Identity, '[^0-9]*$', '');

% Create new column Mean_mm
GS_2021.Mean_mm = GS_2021.Mean_mu/1e3;

% Clear temp vars
clearvars Track_Number Sample_Number numbers name row dataPath opts GS_20211008 GS_20211009


%% Area-wide (mean) 8/9 Oct 2021
f1 = figure('Position',[730, 1846, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h1 = heatmap(GS_2021, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Mean_mm');
% colormap(h1, crameri('lajolla'))
colormap(h1, brewermap([],"YlOrRd"))
clim([.3, 1.4])

h1.Title = [];
h1.XDisplayData = flipud(h1.XDisplayData);
% h1.XDisplayLabels = NaN(size(flipud(h1.XDisplayLabels)));
h1.XDisplayLabels = flipud(h1.XDisplayLabels);
h1.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h1.XLabel = '';
h1.YLabel = '';
h1.FontSize = fontsize*.8;
h1.CellLabelFormat = '%0.2f';
h1.ColorbarVisible = 'off';
h1.GridVisible = 'off';
h1.MissingDataColor = 'w';
h1.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2021-10-08/09', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.2, .4, 'M_{G} (mm)', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h1.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

nexttile(1,[2 5])
errorbar(1.5:10.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
ylim([.25 1.15])
% ylabel('M_{G} (µm)')
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
xlim([.25 1.15])
% xlabel('M_{G} (µm)')
yticks([])


%% Area-wide (std) 8/9 Oct 2021
f2 = figure('Position',[730, 1375, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h2 = heatmap(GS_2021, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Sorting');
% colormap(h2, crameri('imola'))
colormap(h2, brewermap([],"YlGn"))
clim([1.5, 4])

h2.Title = [];
h2.XDisplayData = flipud(h2.XDisplayData);
h2.XDisplayLabels = flipud(h2.XDisplayLabels);
h2.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h2.XLabel = '';
h2.YLabel = '';
h2.FontSize = fontsize*.8;
h2.CellLabelFormat = '%0.2f';
h2.ColorbarVisible = 'off';
h2.GridVisible = 'off';
h2.MissingDataColor = 'w';
h2.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2021-10-08/09', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.4, .4, 'σ_{G}', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h2.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

nexttile(1,[2 5])
errorbar(1.5:10.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
ylim([1.5 3.2])
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
xlim([1.5 3.2])
yticks([])


%% Area-wide (Sk) 8/9 Oct 2021
f3 = figure('Position',[-550, 1846, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h3 = heatmap(GS_2021, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Skewness');
colormap(h3, crameri('tokyo'))
clim([-.2, .6])

h3.Title = [];
h3.XDisplayData = flipud(h3.XDisplayData);
h3.XDisplayLabels = flipud(h3.XDisplayLabels);
h3.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h3.XLabel = '';
h3.YLabel = '';
h3.FontSize = fontsize*.8;
h3.CellLabelFormat = '%0.2f';
h3.ColorbarVisible = 'off';
h3.GridVisible = 'off';
h3.MissingDataColor = 'w';
h3.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2021-10-08/09', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.4, .4, 'Sk_{G}', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h3.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

nexttile(1,[2 5])
errorbar(1.5:10.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
ylim([-.2 .6])
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
xlim([-.2 .6])
yticks([])


%% Area-wide (K) 8/9 Oct 2021
f4 = figure('Position',[-550, 1375, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h4 = heatmap(GS_2021, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Kurtosis');
colormap(h4, crameri('-nuuk'))
% clim([-.2, .6])

h4.Title = [];
h4.XDisplayData = flipud(h4.XDisplayData);
h4.XDisplayLabels = flipud(h4.XDisplayLabels);
h4.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h4.XLabel = '';
h4.YLabel = '';
h4.FontSize = fontsize*.8;
h4.CellLabelFormat = '%0.2f';
h4.ColorbarVisible = 'off';
h4.GridVisible = 'off';
h4.MissingDataColor = 'w';
h4.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2021-10-08/09', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.4, .4, 'K_{G}', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h4.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

nexttile(1,[2 5])
errorbar(1.5:10.5, meanS, stdS, '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
% ylim([.6 1.7])
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
% xlim([.6 1.7])
yticks([])


%% Complete table with missing samples
GS_2022 = [GS_2022(1:5,:); GS_2022(5,:);...                 % R1
    GS_2022(6:10,:); GS_2022(10,:);...                      % R2
    GS_2022(11:14,:); GS_2022(14,:); GS_2022(14,:);...      % R3
    GS_2022(35:end,:);...                                   % R4
    GS_2022(15:18,:); GS_2022(18,:); GS_2022(18,:);...      % R5
    GS_2022(35:end,:);...                                   % R6
    GS_2022(19:23,:); GS_2022(23,:);...                     % R7
    GS_2022(24:28,:); GS_2022(28,:);...                     % R8
    GS_2022(29:end,:)];                                     % R9 & R10
GS_2022{[6 12 17:24 29:36 42 48], 3:end} = NaN;
GS_2022.Sample_Identity(6) = {'R1'};
GS_2022.Sample_Identity(12) = {'R2'};
GS_2022.Sample_Identity(17:18) = {'R3'};
GS_2022.Sample_Identity(19:24) = {'R4'};
GS_2022.Sample_Identity(29:30) = {'R5'};
GS_2022.Sample_Identity(31:36) = {'R6'};
GS_2022.Sample_Identity(42) = {'R7'};
GS_2022.Sample_Identity(48) = {'R8'};


%% Sort table by track number and sample elevation
% Initialize an empty array to store the extracted numbers
Track_Number = zeros(size(GS_2022, 1), 1);

% Loop through each row in the 'name' column and extract the numbers
for row = 1:size(GS_2022, 1)
    name = GS_2022.Sample_Identity{row}; % Get the name from the 'name' column
    numbers = regexp(name, '\d+', 'match'); % Extract all numbers using regular expression
    if ~isempty(numbers)
        Track_Number(row) = str2double(numbers{1}); % Convert the first extracted number to a double
    end
end

% Convert the cell array of numbers to a column in your table
GS_2022.Track_Number = Track_Number;

% Assign sample numbers by elevation
Sample_Number = repelem([1 2 3 4 5 6], 10, 1)';
Sample_Number = Sample_Number(:);
GS_2022.Sample_Number = Sample_Number;

% Remove all characters after the last number in each name
GS_2022.Sample_Identity = regexprep(GS_2022.Sample_Identity, '[^0-9]*$', '');

% Add a '0' before numbers < 10
GS_2022.Sample_Identity = regexprep(GS_2022.Sample_Identity, 'R(\d)$', 'R0$1');

% Remove samples < NAP+0m
GS_2022(6:6:end,:) = [];

% Create new column Mean_mm
GS_2022.Mean_mm = GS_2022.Mean_mu/1e3;

% Clear temp vars
clearvars Track_Number Sample_Number numbers name row


%% Area-wide (mean) 26 Oct 2022
f5 = figure('Position',[1910, 475, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h5 = heatmap(GS_2022, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Mean_mm');
colormap(h5, crameri('lajolla'))
clim([.3, 1.4])

h5.Title = [];
h5.XDisplayData = flipud(h5.XDisplayData);
% h3.XDisplayLabels = NaN(size(flipud(h1.XDisplayLabels)));
h5.XDisplayLabels = flipud(h5.XDisplayLabels);
h5.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h5.XLabel = '';
h5.YLabel = '';
h5.FontSize = fontsize*.8;
h5.CellLabelFormat = '%0.2f';
h5.ColorbarVisible = 'off';
h5.GridVisible = 'off';
h5.MissingDataColor = 'w';
h5.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2022-10-26', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.2, .4, 'M_{G} (mm)', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h5.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

noNaN = ~isnan(meanS);
x = 1.5:10.5;

nexttile(1,[2 5])
errorbar(x(noNaN), meanS(noNaN), stdS(noNaN), '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
ylim([.25 1.15])
% ylabel('M_{G} (µm)')
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
xlim([.25 1.15])
% xlabel('M_{G} (µm)')
yticks([])


%% Calculate difference maps
GS_diff = GS_2022;
GS_diff{:, 3:15} = GS_diff{:, 3:15}-GS_2021{:, 3:15};
GS_diff.Mean_mm = GS_diff.Mean_mu/1e3;


%% Area-wide (mean) 26 Oct 2022 - 8/9 Oct 2021
f6 = figure('Position',[1910, 4, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h6 = heatmap(GS_diff, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Mean_mm');
colormap(h6, crameri('bam'))
% clim([-.5, .8])

h6.Title = [];
h6.XDisplayData = flipud(h6.XDisplayData);
% h4.XDisplayLabels = NaN(size(flipud(h1.XDisplayLabels)));
h6.XDisplayLabels = flipud(h6.XDisplayLabels);
h6.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h6.XLabel = '';
h6.YLabel = '';
h6.FontSize = fontsize*.8;
h6.CellLabelFormat = '%0.2f';
h6.ColorbarVisible = 'off';
h6.GridVisible = 'off';
h6.MissingDataColor = 'w';
h6.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2022-10-26', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.1, .4, '\DeltaM_{G} (mm)', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h6.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

noNaN = ~isnan(meanS);
x = 1.5:10.5;

nexttile(1,[2 5])
errorbar(x(noNaN), meanS(noNaN), stdS(noNaN), '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
ylim([-.6 .6])
% ylabel('M_{G} (µm)')
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
xlim([-.6 .6])
% xlabel('M_{G} (µm)')
yticks([])


%% Area-wide (std) 26 Oct 2022 - 8/9 Oct 2021
f7 = figure('Position',[1910, 4, 1279, 391]);
tiledlayout(6, 6, 'TileSpacing','tight', 'Padding','compact')

nexttile(13,[4 5])
h7 = heatmap(GS_diff, 'Sample_Identity', 'Sample_Number', 'ColorVariable','Sorting');
colormap(h7, crameri('bam'))
% clim([-.5, .8])

h7.Title = [];
h7.XDisplayData = flipud(h7.XDisplayData);
% h4.XDisplayLabels = NaN(size(flipud(h1.XDisplayLabels)));
h7.XDisplayLabels = flipud(h7.XDisplayLabels);
h7.YDisplayLabels = {'+1.00 m', '+0.75 m', '+0.50 m', '+0.25 m', '+0.00 m'};
h7.XLabel = '';
h7.YLabel = '';
h7.FontSize = fontsize*.8;
h7.CellLabelFormat = '%0.2f';
h7.ColorbarVisible = 'off';
h7.GridVisible = 'off';
h7.MissingDataColor = 'w';
h7.MissingDataLabel = 'no data';

nexttile(6,[2 1])
% text(0, .5, '2022-10-26', 'FontSize',fontsize, 'FontWeight','bold', 'EdgeColor','k', 'Margin',6)
text(.3, .4, '\Deltaσ_{G}', 'FontSize',fontsize*.8, 'FontWeight','bold', 'EdgeColor','none', 'Margin',6)
axis off

heatdata = rot90(h7.ColorData, 2);
meanS = mean(heatdata, 1, 'omitmissing');
stdS = std(heatdata, 0, 1, 'omitmissing');
meanT = mean(heatdata, 2, 'omitmissing');
stdT = std(heatdata, 0, 2, 'omitmissing');

noNaN = ~isnan(meanS);
x = 1.5:10.5;

nexttile(1,[2 5])
errorbar(x(noNaN), meanS(noNaN), stdS(noNaN), '-ok', 'LineWidth',3)
yline(mean(meanS, 'omitmissing'), '--k', 'LineWidth',2)
xlim([1 11])
ylim([-1 .8])
% ylabel('M_{G} (µm)')
xticks([])

nexttile(18,[4 1])
errorbar(meanT, 1.5:5.5, stdT, 'horizontal', '-ok', 'LineWidth',3)
xline(mean(meanT, 'omitmissing'), '--k', 'LineWidth',2)
ylim([1 6])
xlim([-1 .8])
% xlabel('M_{G} (µm)')
yticks([])


%% Apply mask
A = load('PHZ_2022_Q2','-mat');
% B = load('PHZ_2020_Q3','-mat');

pgns = getPgons;

mask_scope = inpolygon(A.DEM.X, A.DEM.Y, pgns.scope(:,1), pgns.scope(:,2));
A.DEM.Z(~mask_scope) = NaN;

contourMatrixB = contourc(A.DEM.X(1,:), A.DEM.Y(:,1), A.DEM.Z, [0 0]);
x_0 = contourMatrixB(1, 2:end);
y_0 = contourMatrixB(2, 2:end);
x_0(x_0<min(A.DEM.X(1,:))) = NaN;
y_0(y_0<min(A.DEM.Y(:,1))) = NaN;

contourMatrixC = contourc(A.DEM.X(1,:), A.DEM.Y(:,1), A.DEM.Z, [1 1]);
x_1 = contourMatrixC(1, 2:end);
y_1 = contourMatrixC(2, 2:end);
x_1(x_1<min(A.DEM.X(1,:))) = NaN;
y_1(y_1<min(A.DEM.Y(:,1))) = NaN;

loclabels = string(flipud(h7.XDisplayLabels));
loclabels = [loclabels(5:end); loclabels(1:4)];

%% Horizontal sample locations
X = [S.GS_20211008.xRD_m([1, 5, 9, 14, 19, 24]); S.GS_20211009.xRD_m([1, 6, 11, 16])];
Y = [S.GS_20211008.yRD_m([1, 5, 9, 14, 19, 24]); S.GS_20211009.yRD_m([1, 6, 11, 16])];

X(4) = X(4)+25;  % Move marker downward (i.e. eastward)

% Locations
L1C2 = [117445, 560045];  % L1C2OSSI
L2C9 = [117222, 559793];  % L2C9OSSI
L4C3 = [116125, 558917];  % L4C3OSSI
L5C2 = [115716, 558560];  % L5C2OSSI
L6C2 = [115470, 558176];  % L6C2OSSI

L1C1 = [117421,	560054];  % L1C1VEC
L2C5 = [117198, 559815];  % L2C5SONTEK1
L3C1 = [116839, 559536];  % L3C1VEC
L4C1 = [116103,	558945];  % L4C1VEC
L5C1 = [115670,	558604];  % L5C1VEC
L6C1 = [115402,	558225];  % L6C1VEC

L2C1 = [117158, 559855];  % L2C1KELLER
L2C2 = [117193, 559822];  % L2C2VEC
L2C4 = [117197, 559818];  % L2C4VEC
L2C10 = [117235, 559781]; % L2C10VEC

OSSI = [L6C2; L5C2; L4C3; L2C9; L1C2];
OSSI_names = {"L6", "L5", "L4", "L2", "L1"};

ADV = [L6C1; L5C1; L4C1; L3C1; L2C5; L1C1];
ADV_names = {"L6", "L5", "L4", "L3", "L2", "L1"};


%% Visualisation
f8 = figure('Position',[740, 957, 1719, 1336]);

hold on
% plot(x_0, y_0, '-k', 'LineWidth', 2)
% plot(x_1, y_1, '-k', 'LineWidth', 2)
[M,c] = contourf(A.DEM.X, A.DEM.Y, A.DEM.Z, [0 1]);
colormap([0.9804, 0.9216, 0.8431; 1 1 1]);
% contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [-1.4 -1.4], 'k')

% legend(c,{"sample area"}, 'Location','northwest')

% Instrument locations
scatter(ADV(:,1), ADV(:,2), 150, 'filled', 'MarkerEdgeColor','k')
text(ADV(:,1)+80, ADV(:,2), ADV_names, 'FontSize',fontsize*.8)

% L2 transect
% scatter(L2_loc(1), L2_loc(2), 200, 'LineWidth',5)
% plot([L2C1(1), L2C10(1)], [L2C1(2), L2C10(2)], 'LineWidth',3, 'Color',cbf.bluegreen, 'LineStyle',':')
% text(L2C4(1)-130, L2C4(2)+50, 'L2', 'FontSize',fontsize*.6)

% Sampling locations
scatter(X, Y);
% text(X+60, Y-80, loclabels(1:end), 'FontSize',fontsize*.6)
% text(X-100, Y+20, loclabels(1:end), 'FontSize',fontsize*.6)
text(X-80, Y+20, loclabels(1:end), 'FontSize',fontsize*.8)

% Specify the size of the rectangle
width = 15;   % Width of the rectangle
height = 40;  % Height of the rectangle

% Loop through each data point
for i = 1:length(X)

    if i == 5 || i == 6 || i == 7
        angle = 60;   % Rotation angle in degrees
    else
        angle = 46;
    end
        
        % Coordinates of the rectangle's corners
        rectX = [-width/2, width/2, width/2, -width/2];
        rectY = [-height/2, -height/2, height/2, height/2];
        
        % Rotation matrix
        R = [cosd(angle) -sind(angle); sind(angle) cosd(angle)];
        
        % Rotate and translate the rectangle
        rect = R * [rectX; rectY];
        rect(1, :) = rect(1, :) + X(i);
        rect(2, :) = rect(2, :) + Y(i);
        
        % Draw the rotated rectangle
        patch(rect(1, :), rect(2, :), cbf.bluegreen, 'EdgeColor','k');
end
hold off

xlim(PHZ.xLim);
ylim(PHZ.yLim);

% axis off vis3d
view(46, 90);

axis off equal
view(40, 90)

% Narrow(fontsize)

