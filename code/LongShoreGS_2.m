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
GS_L2 = GS_L2(:, 1:16);
GS_L2 = GS_L2(startsWith(GS_L2.Sample_Identity, ["L2C3","L2C5"]), :);

sampleDate = ['20210928'; '20210930'; '20211001'; '20211006';...
    '20211007'; '20211008'; '20211011'; '20211013'; '20211015'];

% Load remaining L2 sediment data
for i = 1:length(sampleDate)
    dataPath = [folderPath 'grainsizes' filesep 'GS_' sampleDate(i,:) '.csv'];
    opts = detectImportOptions(dataPath);
    opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
    GS_L2_temp = readtable(dataPath, opts);
    GS_L2_temp = GS_L2_temp(:, 1:16);

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
GS_20211008 = GS_20211008(:, 1:16);

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


%% Longshore complete NAP -x m

% Sample_Number == 1 --> NAP + 0.00 m
% Sample_Number == 2 --> NAP - 0.50 m
% Sample_Number == 3 --> NAP - 0.75 m

% Merge tables
Longshore1 = [GS_20210921(:, [1:16, 31]); GS_20210928(:, [1:16, 31]);...
    GS_20211003(:, [1:16, 31]); GS_20211008; GS_L2];

% Convert units of grain diameter from um to mm
Longshore1.Mean_mm = Longshore1.Mean_mu/1e3;

% Select isobath of interest
LongshoreN1 = Longshore1(Longshore1.Sample_Number == 1, :);

% Find rows with dates equal to OLD_DATE and replace them with NEW_DATE
indices = LongshoreN1.Date_ddMMyyyy == datetime('09/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN1.Date_ddMMyyyy(indices) = datetime('08/10/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN1.Date_ddMMyyyy == datetime('20/09/2021', 'Format', 'dd/MM/yyyy');
LongshoreN1.Date_ddMMyyyy(indices) = datetime('21/09/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN1.Date_ddMMyyyy == datetime('01/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN1.Date_ddMMyyyy(indices) = datetime('03/10/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN1.Date_ddMMyyyy == datetime('07/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN1.Date_ddMMyyyy(indices) = datetime('08/10/2021', 'Format', 'dd/MM/yyyy');


%% Visualisation: longshore mean
f1 = figure('Position', [944, 1467, 757, 430]);

h1 = heatmap(LongshoreN1, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Mean_mm');
colormap(h1, brewermap([],"YlOrRd"))
clim([.25, 1.8])

% h1.Title = ['Mean: NAP ', num2str(mean(LongshoreN.zNAP_m), '%.1f'), ' m'];
h1.Title = [];
h1.YDisplayData = locsY(1:end);
h1.YDisplayLabels = locsYnew(1:end);
h1.XDisplayLabels = datetime(h1.XDisplayLabels, 'Format','dd/MM');
h1.XLabel = '';
h1.YLabel = '';
h1.FontSize = fontsize;
h1.CellLabelColor = 'none';
h1.GridVisible = 'off';
h1.MissingDataColor = 'w';
h1.MissingDataLabel = 'no data';

hs1 = struct(h1);
ylabel(hs1.Colorbar, 'M_G (mm)');


%% Visualisation: longshore std
f2 = figure('Position', [944, 957, 757, 430]);

h2 = heatmap(LongshoreN1, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Sorting');
colormap(h2, brewermap([],"YlGnBu"))
clim([1.4, 3])

% h2.Title = ['Sorting: NAP ', num2str(mean(LongshoreN.zNAP_m), '%.1f'), ' m'];
h2.Title = [];
h2.YDisplayData = locsY(1:end);
h2.YDisplayLabels = locsYnew(1:end);
h2.XDisplayLabels = datetime(h2.XDisplayLabels, 'Format','dd/MM');
h2.XLabel = '';
h2.YLabel = '';
h2.FontSize = fontsize;
h2.CellLabelColor = 'none';
h2.GridVisible = 'off';
h2.MissingDataColor = 'w';
h2.MissingDataLabel = 'no data';

hs2 = struct(h2);
ylabel(hs2.Colorbar, '\sigma_G');


%% Compute temporal means
LS_mean = array2table(h1.ColorDisplayData, 'RowNames',h1.YDisplayLabels, 'VariableNames',h1.XDisplayLabels);
LS_std = array2table(h2.ColorDisplayData, 'RowNames',h2.YDisplayLabels, 'VariableNames',h2.XDisplayLabels);

temporalMean = mean(LS_mean{:, ["21/09", "28/09", "03/10", "08/10"]}, 2, 'omitmissing');
temporalMean = array2table(temporalMean, 'RowNames',h1.YDisplayLabels);

% temporalMean = temporalMean./1e3;


%% Longshore complete NAP -x m

% Sample_Number == 1 --> NAP + 0.00 m
% Sample_Number == 2 --> NAP - 0.50 m
% Sample_Number == 3 --> NAP - 0.75 m

% Merge tables
Longshore2 = [GS_20210921(:, [1:16, 31]); GS_20210928(:, [1:16, 31]);...
    GS_20211003(:, [1:16, 31]); GS_20211008; GS_L2];

% Convert units of grain diameter from um to mm
Longshore2.Mean_mm = Longshore2.Mean_mu/1e3;

% Select isobath of interest
LongshoreN2 = Longshore2(Longshore2.Sample_Number == 3, :);

% Find rows with dates equal to OLD_DATE and replace them with NEW_DATE
indices = LongshoreN2.Date_ddMMyyyy == datetime('09/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN2.Date_ddMMyyyy(indices) = datetime('08/10/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN2.Date_ddMMyyyy == datetime('20/09/2021', 'Format', 'dd/MM/yyyy');
LongshoreN2.Date_ddMMyyyy(indices) = datetime('21/09/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN2.Date_ddMMyyyy == datetime('01/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN2.Date_ddMMyyyy(indices) = datetime('03/10/2021', 'Format', 'dd/MM/yyyy');
indices = LongshoreN2.Date_ddMMyyyy == datetime('07/10/2021', 'Format', 'dd/MM/yyyy');
LongshoreN2.Date_ddMMyyyy(indices) = datetime('08/10/2021', 'Format', 'dd/MM/yyyy');


%% Visualisation: longshore mean
f3 = figure('Position', [1702, 1467, 757, 430]);

h3 = heatmap(LongshoreN2, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Mean_mm');
colormap(h3, brewermap([],"YlOrRd"))
clim([.25, 1.8])

% h3.Title = ['Mean: NAP ', num2str(mean(LongshoreN.zNAP_m), '%.1f'), ' m'];
h3.Title = [];
h3.YDisplayData = locsY(1:end);
h3.YDisplayLabels = locsYnew(1:end);
h3.XDisplayLabels = datetime(h3.XDisplayLabels, 'Format','dd/MM');
h3.XLabel = '';
h3.YLabel = '';
h3.FontSize = fontsize;
h3.CellLabelColor = 'none';
h3.GridVisible = 'off';
h3.MissingDataColor = 'w';
h3.MissingDataLabel = 'no data';

hs3 = struct(h3);
ylabel(hs3.Colorbar, 'M_G (mm)');


%% Visualisation: longshore std
f4 = figure('Position', [1702, 957, 757, 430]);

h4 = heatmap(LongshoreN2, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Sorting');
colormap(h4, brewermap([],"YlGnBu"))
clim([1.4, 3])

% h4.Title = ['Sorting: NAP ', num2str(mean(LongshoreN.zNAP_m), '%.1f'), ' m'];
h4.Title = [];
h4.YDisplayData = locsY(1:end);
h4.YDisplayLabels = locsYnew(1:end);
h4.XDisplayLabels = datetime(h4.XDisplayLabels, 'Format','dd/MM');
h4.XLabel = '';
h4.YLabel = '';
h4.FontSize = fontsize;
h4.CellLabelColor = 'none';
h4.GridVisible = 'off';
h4.MissingDataColor = 'w';
h4.MissingDataLabel = 'no data';

hs4 = struct(h4);
ylabel(hs4.Colorbar, '\sigma_G');


%% Visualisation: longshore mean
f5 = figure('Position', [699, 1665, 1279, 628]);
tiledlayout(2,2, 'TileSpacing','compact')

nexttile
h5a = heatmap(LongshoreN1, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Mean_mm');
colormap(h5a, brewermap([],"YlOrRd"))
clim([.25, 1.8])

h5a.Title = [];
h5a.YDisplayData = locsY(1:end);
h5a.YDisplayLabels = locsYnew(1:end);
h5a.XDisplayLabels = repmat({''}, 1, length(h5a.XDisplayData));
h5a.XLabel = '';
h5a.YLabel = '';
h5a.FontSize = fontsize;
h5a.CellLabelColor = 'none';
h5a.GridVisible = 'off';
h5a.MissingDataColor = 'w';
h5a.ColorbarVisible = 'off';


nexttile
h5b = heatmap(LongshoreN2, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Mean_mm');
colormap(h5b, brewermap([],"YlOrRd"))
clim([.25, 1.8])

h5b.Title = [];
h5b.YDisplayData = locsY(1:end);
h5b.YDisplayLabels = repmat({''}, 1, length(h5b.YDisplayData));
h5b.XDisplayLabels = repmat({''}, 1, length(h5b.XDisplayData));
h5b.XLabel = '';
h5b.YLabel = '';
h5b.FontSize = fontsize;
h5b.CellLabelColor = 'none';
h5b.GridVisible = 'off';
h5b.MissingDataColor = 'w';
h5b.MissingDataLabel = 'no data';

hs5b = struct(h5b);
ylabel(hs5b.Colorbar, 'M_G (mm)');


nexttile
h5c = heatmap(LongshoreN1, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Sorting');
colormap(h5c, brewermap([],"YlGnBu"))
clim([1.4, 3])

h5c.Title = [];
h5c.YDisplayData = locsY(1:end);
h5c.YDisplayLabels = locsYnew(1:end);
h5c.XDisplayLabels = datetime(h5c.XDisplayLabels, 'Format','dd/MM');
h5c.XLabel = '';
h5c.YLabel = '';
h5c.FontSize = fontsize;
h5c.CellLabelColor = 'none';
h5c.GridVisible = 'off';
h5c.MissingDataColor = 'w';
h5c.ColorbarVisible = 'off';


nexttile
h5d = heatmap(LongshoreN2, 'Date_ddMMyyyy', 'Sample_Identity', 'ColorVariable','Sorting');
colormap(h5d, brewermap([],"YlGnBu"))
clim([1.4, 3])

h5d.Title = [];
h5d.YDisplayData = locsY(1:end);
h5d.YDisplayLabels = repmat({''}, 1, length(h5d.YDisplayData));
h5d.XDisplayLabels = datetime(h5d.XDisplayLabels, 'Format','dd/MM');
h5d.XLabel = '';
h5d.YLabel = '';
h5d.FontSize = fontsize;
h5d.CellLabelColor = 'none';
h5d.GridVisible = 'off';
h5d.MissingDataColor = 'w';
h5d.MissingDataLabel = 'no data';

hs5d = struct(h5d);
ylabel(hs5d.Colorbar, '\sigma_G');