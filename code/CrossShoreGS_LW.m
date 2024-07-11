%% Initialisation
close all
clear
clc

[~, fontsize, cbf, PHZ, SEDMEX] = sedmex_init;
% fontsize = 30; % ultra-wide screen

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];
sampleIDs = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'};
idxSedi = [2155 2166 2171 2176 2181 2186 2196 2209];


%% Cross-shore profiles
dataPath{6} = [folderPath 'topobathy' filesep 'transects' filesep 'PHZD.nc'];

% Assign variables
z = ncread(dataPath{6}, 'z'); % bed level [m+NAP]
d = ncread(dataPath{6}, 'd'); % cross-shore distance [m]
ID = ncread(dataPath{6}, 'ID'); % profile ID (days since 2020-10-16 00:00:00)

% Timing
startDate = datetime('2020-10-16 00:00:00');
surveyDates = startDate + days(ID);
surveyDates.Format = 'dd MMM yyyy';

% Data selection
track = 2; % L2
survey = [15 18 19 20 21 24 27];
surveyDates = surveyDates(survey);

Z = NaN(length(z), length(survey));
for n = 1:length(survey)
    Z(:,n) = z(:, survey(n), track);
    Z(:,n) = movmean(Z(:,n), 10);
end

xpos = d(idxSedi);

clearvars dataPath startDate n


%% Generate missing data #1

% Dates corresponding to the available data
date_19Sep = datenum('19-Sep-2021'); % idx = 15
date_23Sep = datenum('23-Sep-2021'); % idx = 16

% Date for which we need to interpolate the data
date_20Sep = datenum('20-Sep-2021');

% Available data
Z_19Sep = z(:, 15, track);
Z_23Sep = z(:, 16, track);

% Manually correct outliers through interpolation
Z_23Sep(2334:2341) = interp1([2333, 2342], [Z_23Sep(2333), Z_23Sep(2342)], 2334:2341);

% Smooth data
Z_19Sep = movmean(Z_19Sep, 10);
Z_23Sep = movmean(Z_23Sep, 10);

% Interpolate the data
Z_20Sep = interp1([date_19Sep, date_23Sep], [Z_19Sep'; Z_23Sep'], date_20Sep)';

% Plot the original and interpolated data for visualization
figure;
plot(Z_19Sep, '-o', 'DisplayName', '19 Sep 2021');
hold on;
plot(Z_23Sep, '-o', 'DisplayName', '23 Sep 2021');
plot(Z_20Sep, '-x', 'DisplayName', '20 Sep 2021');
legend;
xlabel('Position along the beach transect');
ylabel('Bed elevation');
title('Interpolated Bed Elevation for 20 Sep 2021');
grid on;

clearvars date_19Sep date_23Sep date_20Sep Z_19Sep Z_23Sep


%% Generate missing data #2

% Dates corresponding to the available data
date_30Sep = datenum('30-Sep-2021'); % idx = 19
date_02Oct = datenum('02-Oct-2021'); % idx = 20

% Date for which we need to interpolate the data
date_01Oct = datenum('01-Oct-2021');

% Available data
Z_30Sep = z(:, 19, track);
Z_02Oct = z(:, 20, track);

% Smooth data
Z_30Sep = movmean(Z_30Sep, 10);
Z_02Oct = movmean(Z_02Oct, 10);

% Interpolate the data
Z_01Oct = interp1([date_30Sep, date_02Oct], [Z_30Sep'; Z_02Oct'], date_01Oct)';

% Plot the original and interpolated data for visualization
figure
plot(Z_30Sep, '-o', 'DisplayName', '30 Sep 2021')
hold on
plot(Z_02Oct, '-o', 'DisplayName', '02 Oct 2021')
plot(Z_01Oct, '-x', 'DisplayName', '01 Oct 2021')
legend
xlabel('Position along the beach transect')
ylabel('Bed elevation')
title('Interpolated Bed Elevation for 01 Oct 2021')
grid on

% Add interpolated data to matrix
Z = [Z(:, 1), Z_20Sep, Z(:, 2:3), Z_01Oct, Z(:, 4:end)];
surveyDates = [surveyDates(1); surveyDates(1)+days(1); surveyDates(2:3); surveyDates(4)-days(1); surveyDates(4:end)];

clearvars z track Z_30Sep Z_20Sep Z_01Oct Z_02Oct date_30Sep date_02Oct date_01Oct


%% Visualisation: L2 profile development
f1 = figure('Position',[740, 1665, 1719, 628]);

hold on
p = gobjects(size(surveyDates));
for n = 1:length(surveyDates)
    p(n) = plot(d, Z(:,n), 'LineWidth',3);
end
hold off

p(2).LineStyle = ":";
p(5).LineStyle = ":";

newcolors = crameri('-roma', length(surveyDates));
colororder(newcolors)

xline(d(idxSedi), '-', sampleIDs, 'FontSize',fontsize, 'LabelHorizontalAlignment','center')

% Tide levels over considered period
yline(.731, '--', 'MHW', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)
yline(.192, '--', 'MSL', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)  
yline(-.535, '--', 'MLW', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)

% Extreme tide levels
yline(0.08, '--', 'storm LW         ', 'FontSize',fontsize, 'Color','red', 'LineWidth',2, 'Alpha',.5) % Oct 01 08:20
yline(1.13, '--', 'storm HW         ', 'FontSize',fontsize, 'Color','red', 'LineWidth',2, 'Alpha',.5) % Sep 30 00:30

xlim([-30, 20])
ylim([-1.7, 1.7])

xlabel('cross-shore distance (m)')
ylabel('bed level (NAP+m)')
legend(string(surveyDates), 'Location','eastoutside')

clearvars n newcolors p d


%% Cross-shore sediment stats
dataPath{1} = [folderPath 'grainsizes' filesep 'GS_20210920.csv'];
dataPath{2} = [folderPath 'grainsizes' filesep 'GS_20210928.csv'];
dataPath{3} = [folderPath 'grainsizes' filesep 'GS_20211001.csv'];
dataPath{4} = [folderPath 'grainsizes' filesep 'GS_20211007.csv'];
dataPath{5} = [folderPath 'grainsizes' filesep 'GS_20211015.csv'];

opts = detectImportOptions(dataPath{1});
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20210920 = readtable(dataPath{1}, opts);
GS_20210920 = [GS_20210920(1:16,:); GS_20210920(1,:); GS_20210920(17:end,:)]; % unchanged L2C2_3
GS_20210920.Sample_Identity(17) = {'L2C2_3'};
GS_20210920.Sample_Number = repelem([1 2 3 4], 8)';
GS_20210920.Sample_Identity = regexprep(GS_20210920.Sample_Identity, '_.*', '');
GS_20210920.Mean_mm = GS_20210920.Mean_mu/1e3;

GS_20210928 = readtable(dataPath{2}, opts);
GS_20210928 = GS_20210928(startsWith(GS_20210928.Sample_Identity, 'L2C'), :);
GS_20210928.Sample_Number = repelem([1 2 3 4], 8)';
GS_20210928.Sample_Identity = regexprep(GS_20210928.Sample_Identity, '_.*', '');
GS_20210928.zNAP_m = repmat(Z(idxSedi, 2), 4, 1);
GS_20210928.Mean_mm = GS_20210928.Mean_mu/1e3;

GS_20211001 = readtable(dataPath{3}, opts);
GS_20211001(1, :) = [];
GS_20211001(endsWith(GS_20211001.Sample_Identity, '1001'), :) = [];
GS_20211001 = [GS_20211001; GS_20211001; GS_20211001; GS_20211001];
GS_20211001{9:end, 6:end} = NaN;
GS_20211001.Sample_Number = repelem([1 2 3 4], 8)';
GS_20211001.zNAP_m = repmat(Z(idxSedi, 3), 4, 1);
GS_20211001.Mean_mm = GS_20211001.Mean_mu/1e3;

GS_20211007 = readtable(dataPath{4}, opts);
GS_20211007.Sample_Number = repelem([1 2 3 4], 8)';
GS_20211007.Sample_Identity = regexprep(GS_20211007.Sample_Identity, '_.*', '');
GS_20211007.zNAP_m = repmat(Z(idxSedi, 4), 4, 1);
GS_20211007.Mean_mm = GS_20211007.Mean_mu/1e3;

GS_20211015 = readtable(dataPath{5}, opts);
GS_20211015 = [GS_20211015(1:7,:); GS_20211015(15,:); GS_20211015(8:end,:)]; % too deep L2C6_1
GS_20211015.Sample_Identity(8) = {'L2C6_1'};
GS_20211015{8, 6:end} = NaN;
GS_20211015 = [GS_20211015; GS_20211015(1:8, :)];
GS_20211015{25:end, 6:end} = NaN;
GS_20211015.Sample_Number = repelem([2 1 3 4], 8)'; % S8 is missing the first time
GS_20211015.Sample_Identity = regexprep(GS_20211015.Sample_Identity, '_.*', '');
GS_20211015.zNAP_m = repmat(Z(idxSedi, 5), 4, 1);
GS_20211015.Mean_mm = GS_20211015.Mean_mu/1e3;

GS_tables = {GS_20210920, GS_20210928, GS_20211001, GS_20211007, GS_20211015};
rowNames = ["GS_20210920"; "GS_20210928"; "GS_20211001"; "GS_20211007"; "GS_20211015"];
columnNames = ["L2C2"; "L2C3"; "L2C3_5"; "L2C4"; "L2C4_5"; "L2C5W"; "L2C5E"; "L2C6"];
samplingDates = datetime(["2021-09-20"; "2021-09-28"; "2021-10-01"; "2021-10-07"; "2021-10-15"], 'Format','dd-MM');

% Replace dots ('.') in sample IDs by underscores ('_')
for i = 1:length(GS_tables)
    % Access the 'name' column
    names = GS_tables{i}.Sample_Identity;

    % Replace '.' with '_' in each name
    newNames = strrep(names, '.', '_');

    % Update the 'name' column in the table
    GS_tables{i}.Sample_Identity = newNames;
end

clearvars folderPath dataPath opts GS_20210920 GS_20210928 GS_20211001 GS_20211007 GS_20211015 i names newNames


%% Isolate the LW-samples
GS_tables_LW = cell(size(GS_tables));

for i = 1:length(GS_tables)

    % Filter rows where Sample_Number is equal to 1
    filteredRows = GS_tables{i}.Sample_Number == 1;
    
    % Create a new table with the filtered rows
    GS_tables_LW{i} = GS_tables{i}(filteredRows, :);

end

clearvars i filteredRows


%% Create new tables
% Number of tables and rows
numTables = numel(GS_tables_LW);
numRows = height(GS_tables_LW{1});

% Initialize the matrix to hold 'column' values
meanMatrix = zeros(numTables, numRows);
sortingMatrix = zeros(numTables, numRows);
skewnessMatrix = zeros(numTables, numRows);
kurtosisMatrix = zeros(numTables, numRows);

% Loop through each table in the cell array
for i = 1:numTables
    % Extract the right column from the table
    mean_values = GS_tables_LW{i}.Mean_mm;
    sorting_values = GS_tables_LW{i}.Sorting;
    skewness_values = GS_tables_LW{i}.Skewness;
    kurtosis_values = GS_tables_LW{i}.Kurtosis;

    % Populate the matrix
    meanMatrix(i, :) = mean_values;
    sortingMatrix(i, :) = sorting_values;
    skewnessMatrix(i, :) = skewness_values;
    kurtosisMatrix(i, :) = kurtosis_values;
end

% Convert the matrix to a table
meanTable = array2table(meanMatrix);
sortingTable = array2table(sortingMatrix);
skewnessTable = array2table(skewnessMatrix);
kurtosisTable = array2table(kurtosisMatrix);

% Assign custom column names
meanTable.Properties.VariableNames = columnNames;
sortingTable.Properties.VariableNames = columnNames;
skewnessTable.Properties.VariableNames = columnNames;
kurtosisTable.Properties.VariableNames = columnNames;

% Assign custom row names
meanTable.Properties.RowNames = rowNames;
sortingTable.Properties.RowNames = rowNames;
skewnessTable.Properties.RowNames = rowNames;
kurtosisTable.Properties.RowNames = rowNames;

clearvars numTables numRows meanMatrix sortingMatrix skewnessMatrix kurtosisMatrix i mean_values sorting_values


%% Visualisation: spatial plots
axs = gobjects(size(GS_tables_LW));
ax_left = gobjects(size(GS_tables_LW));
ax_right = gobjects(size(GS_tables_LW));

f2 = figure(Position=[1222, 957, 1237, 628]);
tiledlayout(3, 2, 'TileSpacing','tight')

for i = 1:length(GS_tables_LW)
    axs(i) = nexttile;
    title(char(samplingDates(i)), FontSize=fontsize*.8)

    yyaxis left
    plot(xpos, meanTable{i,:}, '-ok', 'LineWidth',3, 'MarkerSize',8)
    yline(mean(meanTable{i,:}), '--k', 'LineWidth',1)
    ylim([0 3])
    ax_left(i) = gca;
    ax_left(i).YColor = 'black';
    if i == 3
        ylabel('M_{G} (mm)')
    end
    if i == 2 || i == 4
        yticklabels('')
    end

    yyaxis right
    plot(xpos, sortingTable{i,:}, '-or', 'LineWidth',3, 'MarkerSize',8)
    yline(mean(sortingTable{i,:}), '--r', 'LineWidth',1)
    ylim([1 4])
    ax_right(i) = gca;
    ax_right(i).YColor = 'red';
    if i == 4
        ylabel('\sigma_{G}')
    end
end

xticks(axs, xpos)
xticklabels(axs, sampleIDs)
xticklabels(axs(1:3), '')
yticklabels(axs([1, 3]), '')

xlim(axs, [xpos(1)-1, xpos(end)+1])
grid(axs, 'on')

clearvars i


%% Visualisation: time evolution plots
f3 = figure(Position=[390, 957, 831, 888]);
tiledlayout(4, 2, 'TileSpacing','none')

hold on
for i = 1:length(sampleIDs)
    axs(i) = nexttile;

    yyaxis left
    plot(samplingDates, meanTable{:,i}, '-ok', 'LineWidth',3, 'MarkerSize',8)
    yline(mean(meanTable{:,i}), '--k', 'LineWidth',1)
    ylim([0 3])
    ax_left(i) = gca;
    ax_left(i).YColor = 'black';
    if i == 3
        ylabel('M_{G} (mm)')
    end
    if i == 2 || i == 4 || i == 6 || i == 8
        yticklabels('')
    end

    yyaxis right
    plot(samplingDates, sortingTable{:,i}, '-or', 'LineWidth',3, 'MarkerSize',8)
    yline(mean(sortingTable{:,i}), '--r', 'LineWidth',1)
    ylim([1 4])
    ax_right(i) = gca;
    ax_right(i).YColor = 'red';
    if i == 6
        ylabel('\sigma_{G}')
    end
    if i == 1 || i == 3 || i == 5 || i == 7
        yticklabels('')
    end

    text(samplingDates(end)-days(3.5), 3.5, sampleIDs(i), 'FontSize',fontsize)
end
hold off

xticks(axs, samplingDates(1):days(7):samplingDates(end))
xticklabels(axs(1:6), '')

xlim(axs, [samplingDates(1)-days(1), samplingDates(end)+days(1)])
grid(axs, 'on')

clearvars i ax_right ax_left axs


%% dGSD v.bed-level change 
zMatrix = Z(idxSedi, [2, 3, 5, 8, 9]);
zMatrix = zMatrix';

zTable = array2table(zMatrix);
zTable.Properties.VariableNames = columnNames;
zTable.Properties.RowNames = rowNames;

zDiff = diff(table2array(zTable));
meanDiff = diff(table2array(meanTable));
sortingDiff = diff(table2array(sortingTable));
skewnessDiff = diff(table2array(skewnessTable));
kurtosisDiff = diff(table2array(kurtosisTable));


%% Isolate trends in sorting wrt bed-level change
pos = zDiff>0;

% Sample data
x = zDiff(pos);
y = sortingDiff(pos);

% Manually remove outlier
% x(13) = NaN;
% y(13) = NaN;

% Fit the linear model
mdl = fitlm(x, y);

% Extract the slope and intercept
slope = mdl.Coefficients.Estimate(2);
intercept = mdl.Coefficients.Estimate(1);

% Extract the goodness of fit (R-squared)
R_squared = mdl.Rsquared.Ordinary;

% Display the results
fprintf('Slope: %.4f\n', slope)
fprintf('Intercept: %.4f\n', intercept)
fprintf('R-squared: %.4f\n', R_squared)

% Plot the scatter points
figureRH;
scatter(x, y, 100, 'filled', 'MarkerEdgeColor','k')
yline(0, '--', 'LineWidth',2)
hold on

% Plot the fitted line
x_fit = linspace(min(x), max(x), 100);
y_fit = slope * x_fit + intercept;
plot(x_fit, y_fit, '-r', 'LineWidth', 3)

% Add labels and title
xlabel('\Deltaz (m)')
ylabel('\Delta\sigma_{G}')
title('Linear Fit of Data')
legend('Data points', '', 'Fitted line')

% Add text for the fitted line equation and R-squared
equationText = sprintf('y = %.2fx + %.2f', slope, intercept);
rsquaredText = sprintf('R^2 = %.2f', R_squared);
text(.12, .2, {equationText, rsquaredText}, 'VerticalAlignment','top',...
    'HorizontalAlignment','right', 'FontSize',fontsize*.8)

% Hold off to stop adding to this plot
hold off


%% Visualisation: GSD stats trends wrt bed-level change
f4 = figure(Position=[1, 86, 1470, 754]);
tiledlayout(2,2, 'TileSpacing','loose')

nexttile
scatter(zDiff, meanDiff, 100, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
% xlabel('\Deltaz (m)')
ylabel('\DeltaM_{G} (mm)')

nexttile
scatter(zDiff, sortingDiff, 100, 'filled', 'MarkerEdgeColor','k'); hold on
plot(zDiff(2, 7), sortingDiff(2, 7), 'rx', 'MarkerSize',20, 'LineWidth',2)
plot(x_fit, y_fit, '-r', 'LineWidth', 3); hold off
text(.23, -.2, {equationText, rsquaredText}, 'VerticalAlignment','top',...
    'HorizontalAlignment','right', 'FontSize',fontsize*.6, 'Color','r');
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
% xlabel('\Deltaz (m)')
ylabel('\Delta\sigma_{G}')
legend(sampleIDs, 'Location','eastoutside')

nexttile
scatter(zDiff, skewnessDiff, 100, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
xlabel('\Deltaz (m)')
ylabel('\DeltaSk_{G}')

nexttile
scatter(zDiff, kurtosisDiff, 100, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
xlabel('\Deltaz (m)')
ylabel('\DeltaK_{G}')

newcolors = brewermap(length(meanDiff), 'Greys');
colororder(newcolors)


%% Visualisation: MG and SG trends wrt bed-level change
f5 = figure(Position=[740, 1880, 1719, 413]);
tiledlayout(1,2, 'TileSpacing','compact')

nexttile
scatter(zDiff, meanDiff, 200, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
xlabel('\Deltaz (m)')
ylabel('\DeltaM_{G} (mm)')

nexttile
scatter(zDiff, sortingDiff, 200, 'filled', 'MarkerEdgeColor','k'); hold on
% plot(zDiff(2, 7), sortingDiff(2, 7), 'rx', 'MarkerSize',20, 'LineWidth',2)
% plot(x_fit, y_fit, '-r', 'LineWidth', 3); hold off
% text(.23, -.2, {equationText, rsquaredText}, 'VerticalAlignment','top',...
%     'HorizontalAlignment','right', 'FontSize',fontsize*.6, 'Color','r');
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
xlabel('\Deltaz (m)')
ylabel('\Delta\sigma_{G}')
legend(sampleIDs, 'Location','eastoutside')

newcolors = brewermap(length(meanDiff), 'Greys');
colororder(newcolors)
