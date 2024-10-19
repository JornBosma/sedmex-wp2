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
% survey = [15 18 19 24 27];
survey = [15 18 19 20 21 24 27];
% survey = 15:24;

Z = NaN(length(z), length(survey));
for n = 1:length(survey)
    Z(:,n) = z(:, survey(n), track);
    Z(:,n) = movmean(Z(:,n), 10);
end

xpos = d(idxSedi);

clearvars dataPath z startDate track 


%% L2 profile development
f1 = figure('Position',[740, 1665, 1719, 628]);

hold on
p = nan(length(survey));
% valAlpha = nan(size(survey));
% valCol = nan(length(survey), 3);
% valCol = repmat(linspace(1,.1,length(survey))', 1, 3);
for n = 1:length(survey)
    % valAlpha(n) = (1/length(survey))+(n-1)*(1/length(survey));
    % p(n) = plot(d, Z(:,n), 'LineWidth',3, 'Color',[cbf.orange, valAlpha(n)]);
    % p(n) = plot(d, Z(:,n), 'LineWidth',3, 'Color',valCol(n,:));
    p(n) = plot(d, Z(:,n), 'LineWidth',3);
end
hold off

newcolors = crameri('-roma', length(survey));
% newcolors = brewermap(length(survey), '-RdYlBu');
colororder(newcolors)

xline(d(idxSedi), '-', sampleIDs, 'FontSize',fontsize, 'LabelHorizontalAlignment','center')

% Tide levels over considered period
yline(.731, '--', 'MHW', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)
yline(.192, '--', 'MSL', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)  
yline(-.535, '--', 'MLW', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)

% Extreme tide levels
yline(0.08, '--', 'storm LW         ', 'FontSize',fontsize, 'Color','red', 'LineWidth',2, 'Alpha',.5)  % Oct 01 08:20
yline(1.13, '--', 'storm HW         ', 'FontSize',fontsize, 'Color','red', 'LineWidth',2, 'Alpha',.5) % Sep 30 00:30

xlim([-30, 20])
ylim([-1.7, 1.7])

xlabel('cross-shore distance (m)')
ylabel('bed level (NAP+m)')
legend(string(surveyDates(survey)), 'Location','eastoutside')
% lgd.Position = [0.150 0.172	0.113 0.256];

clearvars p valCol n newcolors


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
GS_20211015.Sample_Number = repelem([1 2 3 4], 8)';
GS_20211015.Sample_Identity = regexprep(GS_20211015.Sample_Identity, '_.*', '');
GS_20211015.zNAP_m = repmat(Z(idxSedi, 5), 4, 1);
GS_20211015.Mean_mm = GS_20211015.Mean_mu/1e3;

GS_tables = {GS_20210920, GS_20210928, GS_20211001, GS_20211007, GS_20211015};
rowNames = ["GS_20210920"; "GS_20210928"; "GS_20211001"; "GS_20211007"; "GS_20211015"];
columnNames = ["L2C2"; "L2C3"; "L2C3_5"; "L2C4"; "L2C4_5"; "L2C5W"; "L2C5E"; "L2C6"];
samplingDates = datetime(["2021-09-20"; "2021-09-28"; "2021-10-01"; "2021-10-07"; "2021-10-15"], 'Format','dd-MM');

clearvars dataPath opts GS_20210920 GS_20210928 GS_20211001 GS_20211007 GS_20211015


%% Replace dots ('.') in sample IDs by underscores ('_')

for i = 1:length(GS_tables)
    % Access the 'name' column
    names = GS_tables{i}.Sample_Identity;

    % Replace '.' with '_' in each name
    newNames = strrep(names, '.', '_');

    % Update the 'name' column in the table
    GS_tables{i}.Sample_Identity = newNames;
end

clearvars i names newNames


%% Calculations

% Initialize tables for mean and standard deviation with sampleIDs as columns
GS_Mean = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Mean_Std = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Sorting = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Sorting_Std = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Count = array2table(zeros(length(GS_tables), length(columnNames)), 'VariableNames', columnNames);

GS_Skewness = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Skewness_Std = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Kurtosis = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);
GS_Kurtosis_Std = array2table(nan(length(GS_tables), length(columnNames)), 'VariableNames',columnNames);

GS_Mean.Properties.RowNames = rowNames;
GS_Mean_Std.Properties.RowNames = rowNames;
GS_Sorting.Properties.RowNames = rowNames;
GS_Sorting_Std.Properties.RowNames = rowNames;
GS_Count.Properties.RowNames = rowNames;

GS_Skewness.Properties.RowNames = rowNames;
GS_Skewness_Std.Properties.RowNames = rowNames;
GS_Kurtosis.Properties.RowNames = rowNames;
GS_Kurtosis_Std.Properties.RowNames = rowNames;

% Loop through each table
for i = 1:length(GS_tables)
    currentTable = GS_tables{i};

    % Loop through each name
    for j = 1:length(columnNames)
        targetID = columnNames{j};

        % Filter the rows where the name matches targetName
        filteredRows = currentTable(strcmp(currentTable.Sample_Identity, targetID), :);
        % filteredRows = filteredRows(1,:);

        % Calculate the mean and standard deviation of the 'Value' column for these rows
        if ~isempty(filteredRows)
            meanValue = mean(filteredRows.Mean_mm, 'omitmissing');
            meanStdValue = std(filteredRows.Mean_mm, 'omitmissing');
            sortingValue = mean(filteredRows.Sorting, 'omitmissing');
            sortingStdValue = std(filteredRows.Sorting, 'omitmissing');
            countValue = sum(~isnan(filteredRows.Mean_mm));

            skewnessValue = mean(filteredRows.Skewness, 'omitmissing');
            skewnessStdValue = std(filteredRows.Skewness, 'omitmissing');
            kurtosisValue = mean(filteredRows.Kurtosis, 'omitmissing');
            kurtosisStdValue = std(filteredRows.Kurtosis, 'omitmissing');
        else
            meanValue = NaN; % Use NaN for no entries
            meanStdValue = NaN;
            sortingValue = NaN;
            sortingStdValue = NaN;
            countValue = 0;  % Zero count for no entries

            skewnessValue = NaN;
            skewnessStdValue = NaN;
            kurtosisValue = NaN;
            kurtosisStdValue = NaN;
        end

        % Assign values to tables
        GS_Mean{i, targetID} = meanValue;
        GS_Mean_Std{i, targetID} = meanStdValue;
        GS_Sorting{i, targetID} = sortingValue;
        GS_Sorting_Std{i, targetID} = sortingStdValue;
        GS_Count{i, targetID} = countValue;

        GS_Skewness{i, targetID} = skewnessValue;
        GS_Skewness_Std{i, targetID} = skewnessStdValue;
        GS_Kurtosis{i, targetID} = kurtosisValue;
        GS_Kurtosis_Std{i, targetID} = kurtosisStdValue;
    end
end

GS_Mean_mean = mean(GS_Mean, 2);
GS_Sorting_mean = mean(GS_Sorting, 2);

GS_Skewness_mean = mean(GS_Skewness, 2);
GS_Kurtosis_mean = mean(GS_Kurtosis, 2);

% dGS_Mean = NaN(height(GS_Mean), width(GS_Mean)-1);
% GS_Mean_grad = NaN(height(GS_Mean), width(GS_Mean)-1);
% dXpos = diff(xpos);
% 
% for i = 1:height(GS_Mean)
%     % Calculate differences
%     dGS_Mean(i,:) = diff(GS_Mean{i,:});
% 
%     % Compute gradient
%     GS_Mean_grad(i,:) = dGS_Mean(i,:) ./ dXpos';
% end
% 
% GS_Mean_grad = GS_Mean_grad*1e3;  % Convert unit from mm to mu
% GS_Mean_grad_mean_x = mean(GS_Mean_grad, 2);
% GS_Mean_gradient = [GS_Mean_grad, GS_Mean_grad_mean_x];
% GS_Mean_grad_mean_t = mean(GS_Mean_gradient, 1);
% GS_Mean_gradient = [GS_Mean_gradient; GS_Mean_grad_mean_t];
% GS_Mean_gradient = round(GS_Mean_gradient);

clearvars filteredRow i j meanValue meanstdValue sortingValue sortingstdValue countValue targetID dGS_Mean GS_Mean_grad GS_Mean_grad_mean_x GS_Mean_gradient GS_Mean_grad_mean_t


%% Visualisation
axs = gobjects(size(GS_tables));
ax_left = gobjects(size(GS_tables));
ax_right = gobjects(size(GS_tables));

f2 = figure(Position=[1222, 957, 1237, 628]);
tiledlayout(3, 2, 'TileSpacing','tight')

for i = 1:length(GS_tables)
    axs(i) = nexttile;
    title(char(samplingDates(i)), FontSize=fontsize*.8)

    yyaxis left
    errorbar(xpos, GS_Mean{i,:}, GS_Mean_Std{i,:}, '-ok', 'LineWidth',3)
    yline(GS_Mean_mean.mean(i), '--k', 'LineWidth',1)
    ylim([0 3])
    ax_left(i) = gca;
    ax_left(i).YColor = 'black';
    if i == 3
        ylabel('M_{G} (mm)')
    end
    if i == 2 || i == 4
        yticklabels('')
    end
    text(xpos(1)+.6, 2.5, ['(n = ', mat2str(GS_Count{i,1}), ')'], 'Fontsize',fontsize*.6)

    yyaxis right
    errorbar(xpos, GS_Sorting{i,:}, GS_Sorting_Std{i,:}, '-or', 'LineWidth',3)
    yline(GS_Sorting_mean.mean(i), '--r', 'LineWidth',1)
    ylim([1 4])
    ax_right(i) = gca;
    ax_right(i).YColor = 'red';
    if i == 4
        ylabel('σ_{G}')
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
tl = tiledlayout(4, 2, 'TileSpacing','tight');

hold on
for i = 1:length(sampleIDs)
    axs(i) = nexttile;

    yyaxis left
    % plot(samplingDates, GS_Mean{:,i}, '-ok', 'LineWidth',3)
    errorbar(samplingDates, GS_Mean{:,i}, GS_Mean_Std{:,i}, '-ok', 'LineWidth',3)
    yline(mean(GS_Mean{:,i}), '--k', 'LineWidth',1)
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
    % plot(samplingDates, GS_Sorting{:,i}, '-or', 'LineWidth',3)
    errorbar(samplingDates, GS_Sorting{:,i}, GS_Sorting_Std{:,i}, '-or', 'LineWidth',3)
    yline(mean(GS_Sorting{:,i}), '--r', 'LineWidth',1)
    ylim([1 4])
    ax_right(i) = gca;
    ax_right(i).YColor = 'red';
    if i == 6
        ylabel('σ_{G}')
    end
    if i == 1 || i == 3 || i == 5 || i == 7
        yticklabels('')
    end

    text(samplingDates(end)-days(3.5), 3.5, sampleIDs(i), 'FontSize',fontsize)
end
hold off

% ylabel(tl, 'M_{G} (mm)', 'FontSize',fontsize)

xticks(axs, samplingDates(1):days(7):samplingDates(end))
xticklabels(axs(1:6), '')

xlim(axs, [samplingDates(1)-days(1), samplingDates(end)+days(1)])
grid(axs, 'on')

clearvars i ax_right ax_left axs


%% dGSD v.bed-level change 
bed_level = Z(idxSedi, [1, 2, 4, 6, 7]);
bed_level = bed_level';
bed_level_change = diff(bed_level, [], 1);

dGS_Mean = diff(table2array(GS_Mean));
dGS_Sorting = diff(table2array(GS_Sorting));
dGS_Skewness = diff(table2array(GS_Skewness));
dGS_Kurtosis = diff(table2array(GS_Kurtosis));


%% Isolate trends in sorting wrt bed-level change
pos = bed_level_change>0;

% Sample data
x = bed_level_change(pos);
y = dGS_Sorting(pos);

% Manually remove outlier
x(13) = NaN;
y(13) = NaN;

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
ylabel('\Deltaσ_{G}')
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
scatter(bed_level_change, dGS_Mean, 100, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
% xlabel('\Deltaz (m)')
ylabel('\DeltaM_{G} (mm)')

nexttile
scatter(bed_level_change, dGS_Sorting, 100, 'filled', 'MarkerEdgeColor','k'); hold on
plot(bed_level_change(2, 7), dGS_Sorting(2, 7), 'rx', 'MarkerSize',20, 'LineWidth',2)
plot(x_fit, y_fit, '-r', 'LineWidth', 3); hold off
text(.29, .35, {equationText, rsquaredText}, 'VerticalAlignment','top',...
    'HorizontalAlignment','right', 'FontSize',fontsize*.6, 'Color','r');
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
% xlabel('\Deltaz (m)')
ylabel('\Deltaσ_{G}')
legend(sampleIDs, 'Location','eastoutside')

nexttile
scatter(bed_level_change, dGS_Skewness, 100, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
xlabel('\Deltaz (m)')
ylabel('\DeltaSk_{G}')

nexttile
scatter(bed_level_change, dGS_Kurtosis, 100, 'filled', 'MarkerEdgeColor','k')
xline(0, '--', 'LineWidth',2)
yline(0, '--', 'LineWidth',2)
xlabel('\Deltaz (m)')
ylabel('\DeltaK_{G}')

newcolors = brewermap(length(dGS_Mean), 'Greys');
colororder(newcolors)
