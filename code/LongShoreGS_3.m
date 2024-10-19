%% Initialisation
close all
clear
clc

[~, ~, cbf, ~, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];


%% GS_20210921

% Load sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210921.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_20210921 = readtable(dataPath, opts);

% Prepare table
GS_20210921 = GS_20210921(~startsWith(GS_20210921.Sample_Identity, 'L2C'), :);
% GS_20210921.Sample_Number = [1; 3; 1; 2; 3; 1; 3; 1; 3; 1; 3; 1; 3];
GS_20210921.Sample_Identity = regexprep(GS_20210921.Sample_Identity, '_.*', '');

clearvars dataPath opts


%% GS_20210928

% Load sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210928.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_20210928 = readtable(dataPath, opts);

% Prepare table
GS_20210928 = GS_20210928(~startsWith(GS_20210928.Sample_Identity, 'L2C'), :);
% GS_20210928.Sample_Number = [1; 3; 1; 2; 3; 1; 3; 1; 3; 1; 3; 1; 3];
GS_20210928.Sample_Identity = regexprep(GS_20210928.Sample_Identity, '_.*', '');

clearvars dataPath opts


%% GS_L2

% Initialise table
dataPath = [folderPath 'grainsizes' filesep 'GS_20210920.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_L2 = readtable(dataPath, opts);
GS_L2 = GS_L2(startsWith(GS_L2.Sample_Identity, ["L2C3","L2C5"]), :);

% Reduce sieve tower
GS_L2.Sieve_2000mu_g = GS_L2.Sieve_2000mu_g + GS_L2.Sieve_2800mu_g;
GS_L2.Sieve_2800mu_g = [];

GS_L2.Sieve_1180mu_g = GS_L2.Sieve_1180mu_g + GS_L2.Sieve_1400mu_g + GS_L2.Sieve_1700mu_g;
GS_L2.Sieve_1400mu_g = [];
GS_L2.Sieve_1700mu_g = [];

GS_L2.Sieve_710mu_g = GS_L2.Sieve_710mu_g + GS_L2.Sieve_850mu_g;
GS_L2.Sieve_850mu_g = [];

GS_L2.Sieve_500mu_g = GS_L2.Sieve_500mu_g + GS_L2.Sieve_600mu_g;
GS_L2.Sieve_600mu_g = [];

GS_L2.Sieve_180mu_g = GS_L2.Sieve_180mu_g + GS_L2.Sieve_212mu_g;
GS_L2.Sieve_212mu_g = [];

GS_L2.Sieve_125mu_g = GS_L2.Sieve_125mu_g + GS_L2.Sieve_150mu_g;
GS_L2.Sieve_150mu_g = [];

GS_L2.Sieve_63mu_g = GS_L2.Sieve_63mu_g + GS_L2.Sieve_75mu_g + GS_L2.Sieve_90mu_g;
GS_L2.Sieve_75mu_g = [];
GS_L2.Sieve_90mu_g = [];

GS_L2.Sieve_Pan_g = GS_L2.Sieve_Pan_g + GS_L2.Sieve_53mu_g;
GS_L2.Sieve_53mu_g = [];

% Load remaining L2 sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210928.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_L2_temp = readtable(dataPath, opts);

GS_L2 = [GS_L2; GS_L2_temp];
GS_L2 = GS_L2(startsWith(GS_L2.Sample_Identity, ["L2C3_","L2C5W","L2C5_"]), :);
GS_L2.Sample_Identity = regexprep(GS_L2.Sample_Identity, ["_.*", "W"], '');

% Group the table by 'Sample_Identity' and 'Date_ddMMyyyy' and calculate the mean of other columns
groupedTable = groupsummary(GS_L2, {'Sample_Identity', 'Date_ddMMyyyy'}, 'mean');
groupedTable.GroupCount = []; % remove extra column

% Set column names of the second table to be the same as the first table
groupedTable.Properties.VariableNames = GS_L2.Properties.VariableNames;

GS_L2 = groupedTable;
GS_L2.Sample_Identity = regexprep(GS_L2.Sample_Identity, "C.*", '');

clearvars GS_L2_temp sampleIdentity groupedTable dataPath i opts


%% Merge data
GS_20210921 = [GS_20210921(1:9, :); GS_L2([1, 3], :); GS_20210921(10:end, :)];
GS_20210928 = [GS_20210928(1:9, :); GS_L2([2, 4], :); GS_20210928(10:end, :)];

% Extract useful data
sieveSizes_20210921 = fliplr(extractSieveNumbers(GS_20210921));
massRetained_20210921 = fliplr(GS_20210921{:, 17:end});
sieveSizes_20210928 = fliplr(extractSieveNumbers(GS_20210928));
massRetained_20210928 = fliplr(GS_20210928{:, 17:end});

tables = {GS_20210921, GS_20210928}; % Collection of tables
sieveSizes = {sieveSizes_20210921, sieveSizes_20210928}; % Corresponding sieve sizes
massRetained = {massRetained_20210921, massRetained_20210928};

clearvars GS_20210921 GS_20210928 GS_L2 sieveSizes_20210921 sieveSizes_20210928 sieveSizes_L2 massRetained_20210921 massRetained_20210928 massRetained_L2 


%% Visualisation
% n = 1; % Select the table index
% f1 = gobjects(height(massRetained{n}), 1);
% 
% for i = 1%:height(tables{n}) % Loop through each row of the selected table
% 
%     sieveDiam = [sieveSizes{n}.*1e-3, 16]; % in mm
%     sieveData = massRetained{n}(i, :);
% 
%     % % Relative?
%     % binWidth = diff(sieveDiam);
%     % relativeMass = sieveData ./ binWidth;
% 
%     f1(i) = figure('Position', [740, 1873, 560, 420]);
%     histogram('BinEdges',sieveDiam, 'BinCounts',sieveData, ...
%         'Normalization','percentage', 'FaceAlpha',.5, 'FaceColor',cbf.orange)
% 
%     xlabel('Particle diameter (mm)');
%     set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
%     xticks([0.1, 1, 10])
%     ytickformat("percentage")
% 
% end
% 
% clearvars n i sieveDiam sieveData binWidth relativeMass


%% Visualisation (0 v. -0.8)
n = 2; % Select the table index
sieveDiam = [sieveSizes{n}.*1e-3, 16]; % in mm
sieveData = massRetained{n}([1:3, 5:end], :);
f2 = gobjects(height(sieveData)/2, 1);

j = 1;
for i = 1:length(f2)
    f2(i) = figure('Position', [740, 1873, 560, 420]); % Create new figure

    hold on
    histogram('BinEdges',sieveDiam, 'BinCounts',sieveData(j+1, :), ...
    'Normalization','percentage', 'FaceAlpha',1, 'FaceColor',cbf.orange)

    histogram('BinEdges',sieveDiam, 'BinCounts',sieveData(j, :), ...
    'Normalization','percentage', 'FaceAlpha',.6, 'FaceColor',cbf.yellow)
    hold off

    set(gca, 'XScale', 'log'); % Set x-axis to logarithmic scale
    xticks([0.1, 1, 10])
    % ytickformat("percentage")
    ylim([0, 30])

    legend('NAP–0.8m', 'NAP+0.0m')
    % xlabel('particle diameter (mm)')
    % ylabel('class weight (%)')
    box off % Remove the box around the plot
    % axis off % Uncomment if you want to turn off the axis

    j = j + 2;

end

clearvars n i j sieveDiam sieveData

