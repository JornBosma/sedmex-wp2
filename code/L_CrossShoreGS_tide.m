%% Initialisation
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];
sampleIDs = {'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S7', 'S8'};
idxSedi = [2155 2166 2171 2176 2181 2186 2196 2209];
xLocs = {'L2C2','L2C3','L2C3_5','L2C4','L2C4_5','L2C5W','L2C5E','L2C6'};

samplingDates = datetime(["2021-09-20"; "2021-09-28"; "2021-10-01"; "2021-10-07"; "2021-10-15"], 'Format','dd-MMM-yyyy');

scrapeIDs = {'A1', 'A2', 'A3', 'B1', 'B2', 'B3'};
idxScrp = [2111 2160 2195 2100 2137 2156];
zScrp = [0.842, 0.179, -0.738, 1.058, 0.5, 0.168];
MgScrp = [0.7663, 0.8251, 1.1917, 0.6488, 1.1421, 1.2748];
SgScrp = [2.3424, 2.4073, 2.0377, 2.3947, 1.9592, 2.3955];

% L2 sampling times
Sep20 = datetime({'20-Sep-2021 09:15:00'; '20-Sep-2021 12:10:00'; '20-Sep-2021 15:15:00'; '20-Sep-2021 18:40:00'});
Sep28 = datetime({'28-Sep-2021 10:00:00'; '28-Sep-2021 13:00:00'; '28-Sep-2021 16:00:00'; '28-Sep-2021 19:00:00'});
Oct01 = datetime({'01-Oct-2021 09:15:00'});
Oct07 = datetime({'07-Oct-2021 10:00:00'; '07-Oct-2021 13:00:00'; '07-Oct-2021 16:00:00'; '07-Oct-2021 19:00:00'}); % Only first and last times were written down
Oct15 = datetime({'15-Oct-2021 10:20:00'; '15-Oct-2021 13:25:00'; '15-Oct-2021 15:30:00'});

% Combine all sampling times into a single array
allSamplingTimes = [Sep20; Sep28; Oct01; Oct07; Oct15];

clearvars Sep20 Sep28 Oct01 Oct07 Oct15


%% Water level data
instrument = 'L2C6OSSI';
filename = fullfile(folderPath, 'hydrodynamics', 'pressuresensors',...
    instrument, ['tailored_', instrument, '.nc']);

% instrument = 'L2C7ADCP';
% filename = fullfile(folderPath, 'hydrodynamics', 'ADCP',...
%     ['tailored_', instrument, '.nc']);

% instrument = 'L2C10VEC';
% filename = fullfile(folderPath, 'hydrodynamics', 'ADV',...
%     instrument, ['tailored_', instrument, '.nc']);

info = ncinfo(filename);

% Defined starting time
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
t = ncread(filename, 't');      % seconds since 2021-09-01
Time = t0 + seconds(t);
eta = ncread(filename, 'zs');   % water level [NAP+m]
Hm0 = ncread(filename, 'Hm0');  % wave height [m]
L2C6 = timetable(Time, eta, Hm0);

% Retiming the water levels timetable to match the sampling times
alignedWaterLevels = retime(L2C6, allSamplingTimes, 'linear');

% Display the aligned water levels at the sampling times
disp(alignedWaterLevels);

clearvars allSamplingTimes info t t0 Time filename eta Hm0


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

clearvars dataPath startDate n ID


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

% Add interpolated data to matrix
Z = [Z(:, 1), Z_20Sep, Z(:, 2:3), Z_01Oct, Z(:, 4:end)];
surveyDates = [surveyDates(1); surveyDates(1)+days(1); surveyDates(2:3); surveyDates(4)-days(1); surveyDates(4:end)];

% Create colour scheme
surveyColours = crameri('-roma', length(surveyDates));

clearvars z track Z_30Sep Z_20Sep Z_01Oct Z_02Oct date_30Sep date_02Oct date_01Oct


%% Visualisation: L2 profile development
f0a = figure('Position',[740, 1665, 1719, 628]);

hold on
p = gobjects(size(surveyDates));
for n = 1:length(surveyDates)
    p(n) = plot(d, Z(:,n), 'LineWidth',3);
end
scatter(d(idxScrp([1, 4, 5])), zScrp([1, 4, 5]), 200, 'vk', 'LineWidth',3)
text(d(idxScrp([1, 4, 5]))-.4, zScrp([1, 4, 5])-.2, scrapeIDs([1, 4, 5]), 'FontSize',fontsize*.8)
% text(d(idxScrp([2, 3, 6]))-.4, zScrp([2, 3, 6])+.15, scrapeIDs([2, 3, 6]), 'FontSize',fontsize*.8)
hold off

p(2).LineStyle = ":";
p(5).LineStyle = ":";

colororder(surveyColours)

% xline(d(idxScrp), '-', scrapeIDs, 'FontSize',fontsize, 'LabelHorizontalAlignment','center')
xline(d(idxSedi), '-', sampleIDs, 'FontSize',fontsize, 'LabelHorizontalAlignment','center')

% Tide levels over considered period
yline(.731, '--', 'MHW', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)
yline(.192, '--', 'MSL', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)  
yline(-.535, '--', 'MLW', 'FontSize',fontsize, 'Color',cbf.grey, 'LineWidth',2, 'Alpha',.5)

% Extreme tide levels
yline(0.08, '--', 'storm LW         ', 'FontSize',fontsize, 'Color','red', 'LineWidth',2, 'Alpha',.5) % Oct 01 08:20
yline(1.13, '--', 'storm HW         ', 'FontSize',fontsize, 'Color','red', 'LineWidth',2, 'Alpha',.5) % Sep 30 00:30

xlim([-30, 10])
ylim([-1.7, 1.7])

xlabel('cross-shore distance (m)')
ylabel('bed level (NAP+m)')
legend(string(surveyDates, 'dd MMM'), 'Location','eastoutside')

clearvars n p


%% Visualisation: Water levels during sampling
f0b = gobjects(size(samplingDates));
p = gobjects(size(samplingDates));
s1 = gobjects(size(samplingDates));
s2 = gobjects(size(samplingDates));

m = 1;

for n = [2, 3, 5, 8, 9]
    % Visualisation: L2 profile development
    f0b(m) = figure('Position', [740, 1665, 1719, 628]);
    
    hold on
    p(m) = plot(d, Z(:,n), 'LineWidth', 3);
    s1(m) = scatter(d(idxSedi), Z(idxSedi, n), 200, 'vk', 'LineWidth', 3);
    text(d(idxSedi)-.5, Z(idxSedi, n)-.2, sampleIDs, 'FontSize', fontsize * .8)

    if n == 8
        s2(m) = scatter(d(idxScrp(1:3)), zScrp(1:3), 200, 'dk', 'LineWidth',3);
        text(d(idxScrp(1:3))-.4, zScrp(1:3)+.2, scrapeIDs(1:3), 'FontSize',fontsize*.8)
    end
    if n == 9
        s2(m) = scatter(d(idxScrp(4:6)), zScrp(4:6), 200, 'dk', 'LineWidth',3);
        text(d(idxScrp(4:6))-.4, zScrp(4:6)+.2, scrapeIDs(4:6), 'FontSize',fontsize*.8)
    end
    hold off
    
    colororder(surveyColours(n, :))
    
    % Filter water levels for the current survey date
    currentDate = dateshift(surveyDates(n), 'start', 'day');
    nextDate = currentDate + days(1);
    currentWaterLevels = alignedWaterLevels(alignedWaterLevels.Time >= currentDate & alignedWaterLevels.Time < nextDate, :);
    
    % Tide level during sampling
    yl = yline(currentWaterLevels.eta, '--', string(currentWaterLevels.Time, "HH:mm"), 'FontSize', fontsize, 'Color', cbf.skyblue, 'LineWidth', 2, 'Alpha', .5);
    if n == 3
        yl(1).LabelVerticalAlignment = "bottom";
    end
    if n == 9
        yl(2).LabelVerticalAlignment = "bottom";
    end

    xlim([-30, 10])
    ylim([-1.7, 1.7])
    
    xlabel('cross-shore distance (m)')
    ylabel('bed level (NAP+m)')

    try
        legend([p(m), s1(m), s2(m)], {string(surveyDates(n), 'dd MMM'), 'sample', 'scraper'}, 'Location', 'eastoutside')
    catch
        legend([p(m), s1(m)], {string(surveyDates(n), 'dd MMM'), 'sample'}, 'Location', 'eastoutside')
    end

    m = m + 1;
end
p(1).LineStyle = ":";
p(3).LineStyle = ":";

clearvars n m newcolors p d s1 s2 yl year currentDate currentWaterLevels nextDate


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


%% Organisation
% Extract the date part from the 'Time' column and add it as a new variable
alignedWaterLevels.Date = dateshift(alignedWaterLevels.Time, 'start', 'day');

% Convert Date to categorical for easier counting
alignedWaterLevels.Date = categorical(alignedWaterLevels.Date);

% Count occurrences of each date
[counts, dates] = histcounts(alignedWaterLevels.Date);

% Display the results
countT = table(dates', counts', 'VariableNames', {'Date', 'Count'});
disp(countT);

% Extract the bed levels at sample locations
sampleZ = Z(idxSedi, [2, 3, 5, 8, 9]);

clearvars dates counts


%% Visualisation: spatial plots (mean)
axs = gobjects(size(GS_tables));
tideColours = brewermap(4, '-PuOr');

f1b = figure(Position=[1222, 1665, 1237, 628]);
tiledlayout(3, 2, 'TileSpacing','compact')

for i = 1:length(GS_tables)

    % Get the corresponding water level for this sample time
    sampleTimes = alignedWaterLevels.Time(alignedWaterLevels.Date == categorical(samplingDates(i)));
    waterLevels = alignedWaterLevels.eta(alignedWaterLevels.Date == categorical(samplingDates(i)));

    axs(i) = nexttile;
    title(string(samplingDates(i), 'dd MMM'), FontSize=fontsize*.8)
    hold on
    for j = 1:countT.Count(i)
        data = GS_tables{i}(GS_tables{i}.Sample_Number == j, :);
        
        if j > 1
            % Indicate sample locations which have not been submerged since
            % the previous sampling moment
            dry = sampleZ(:, i) > waterLevels(j)+0.1 & sampleZ(:, i) > waterLevels(j-1)+0.1;
            % Bed level has stayed above water level, plot with thick red edge
            scatter(xpos(dry), data.Mean_mm(dry), 300, 'filled', ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'g', ...
                'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
        end

        plot(xpos, data.Mean_mm, '-o', 'LineWidth',3, 'MarkerSize',8, 'Color',tideColours(j, :))

    end

    if i == 3
        ylabel('M_{G} (mm)')
    end
    if i == 2 || i == 4
        yticklabels('')
    end

    yline(mean(GS_tables{i}.Mean_mm, 'all', 'omitmissing'), '--k', 'LineWidth',2)
    hold off

    legend([string(sampleTimes, 'HH:mm'); 'mean'], 'Location','northeastoutside')
end

xlim(axs, [xpos(1)-1, xpos(end)+1])
xticks(axs, xpos)
xticklabels(axs, sampleIDs)
xticklabels(axs(1:3), '')

ylim(axs, [0, 3.5])
yticks(axs, 0:3)

grid(axs, 'on')

clearvars i j dry sampleTimes waterLevels axs


%% Visualisation: spatial plots (sorting)
axs = gobjects(size(GS_tables));
tideColours = brewermap(4, '-PuOr');

f1a = figure(Position=[1222, 1665, 1237, 628]);
tiledlayout(3, 2, 'TileSpacing','compact')

for i = 1:length(GS_tables)

    % Get the corresponding water level for this sample time
    sampleTimes = alignedWaterLevels.Time(alignedWaterLevels.Date == categorical(samplingDates(i)));
    waterLevels = alignedWaterLevels.eta(alignedWaterLevels.Date == categorical(samplingDates(i)));

    axs(i) = nexttile;
    title(string(samplingDates(i), 'dd MMM'), FontSize=fontsize*.8)
    hold on
    for j = 1:countT.Count(i)
        data = GS_tables{i}(GS_tables{i}.Sample_Number == j, :);
        
        if j > 1
            % Indicate sample locations which have not been submerged since
            % the previous sampling moment
            dry = sampleZ(:, i) > waterLevels(j)+0.1 & sampleZ(:, i) > waterLevels(j-1)+0.1;
            % Bed level has stayed above water level, plot with thick red edge
            scatter(xpos(dry), data.Sorting(dry), 300, 'filled', ...
               'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'g', ...
                'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
        end

        plot(xpos, data.Sorting, '-o', 'LineWidth',3, 'MarkerSize',8, 'Color',tideColours(j, :))

    end

    if i == 3
        ylabel('σ_{G}')
    end
    if i == 2 || i == 4
        yticklabels('')
    end

    yline(mean(GS_tables{i}.Sorting, 'all', 'omitmissing'), '--k', 'LineWidth',2)
    hold off

    legend([string(sampleTimes, 'HH:mm'); 'mean'], 'Location','northeastoutside')
end

xlim(axs, [xpos(1)-1, xpos(end)+1])
xticks(axs, xpos)
xticklabels(axs, sampleIDs)
xticklabels(axs(1:3), '')

ylim(axs, [1.5, 4.2])
yticks(axs, 2:4)

grid(axs, 'on')

clearvars i j dry sampleTimes waterLevels axs


%% Visualisation: time evolution plots
axs = gobjects(size(GS_tables));
ax_left = gobjects(size(GS_tables));
ax_right = gobjects(size(GS_tables));
f2 = gobjects(size(GS_tables));
tidalVariation_M = nan(length(GS_tables), length(xpos));
tidalVariation_S = nan(length(GS_tables), length(xpos));

dryThreshold = 0.2;

for i = 1:length(GS_tables)
    f2(i) = figure(Position=[740, 1405, 831, 888]);
    tiledlayout(5,2, 'TileSpacing','compact')

    % Get the corresponding water level for this sample time
    sampleTimes = alignedWaterLevels.Time(alignedWaterLevels.Date == categorical(samplingDates(i)));
    waterLevels = alignedWaterLevels.eta(alignedWaterLevels.Date == categorical(samplingDates(i)));

    % Specify the date of interest
    targetDate = samplingDates(i);

    % Extract the data for the specified date
    dailyData = L2C6(L2C6.Time >= targetDate & L2C6.Time < targetDate + days(1), :);

    axs(1) = nexttile;
    plot(dailyData.Time, dailyData.eta, 'LineWidth',3)
    ylabel('η (NAP+m)')
    ylim([-1.2, 1])

    axs(2) = nexttile;
    yyaxis left
    ax_left(2) = gca;
    ax_left(2).YColor = 'k';
    yticklabels({})
    yyaxis right
    plot(dailyData.Time, dailyData.Hm0, 'LineWidth',3, 'Color',cbf.blue)
    ax_right(2) = gca;
    ax_right(2).YColor = 'k';
    ylabel('H_{m0} (m)')
    ylim([0, .5])
    
    % if i == 3
    %     text(sampleTimes(1)+hours(5.8), .08, instrument, 'FontSize',fontsize*.8)
    % else
    %     text(sampleTimes(1)+hours(5.8), .42, instrument, 'FontSize',fontsize*.8)
    % end

    hold on
    for j = 1:length(xLocs)
        axs(j+2) = nexttile;        
        data = GS_tables{i}(strcmp(GS_tables{i}.Sample_Identity, xLocs{j}), :);
        sampleHeight = sampleZ(j, i); % assume bed level remains constant throughout tidal cycle
        
        yyaxis left
        for k = 2:length(waterLevels)
            % Indicate sample locations which have not been submerged since
            % the previous sampling moment
            dry = sampleZ(j, i) > waterLevels(k)+dryThreshold & sampleHeight > waterLevels(k-1)+dryThreshold;
            % Bed level has stayed above water level, plot with thick red edge
            if dry
                scatter(sampleTimes(k), data.Mean_mm(k), 300, 'filled', ...
                   'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'g', ...
                    'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
            end
        end

        plot(sampleTimes, data.Mean_mm(1:length(sampleTimes)), '-ok', 'LineWidth',3, 'MarkerSize',8)
        yline(mean(data.Mean_mm, 'all', 'omitmissing'), '--k', 'LineWidth',2)
        ylim([0, 3.5])
        yticks(0:3)
        ax_left(j) = gca;
        ax_left(j).YColor = 'black';
        if j == 3
            ylabel('M_{G} (mm)')
        end
        if j == 2 || j == 4 || j == 6 || j == 8
            yticklabels('')
        end
    
        yyaxis right
        for k = 2:length(waterLevels)
            % Indicate sample locations which have not been submerged since
            % the previous sampling moment
            dry = sampleZ(j, i) > waterLevels(k)+dryThreshold & sampleHeight > waterLevels(k-1)+dryThreshold;
            % Bed level has stayed above water level, plot with thick red edge
            if dry
                scatter(sampleTimes(k), data.Sorting(k), 300, 'filled', ...
                   'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'g', ...
                    'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off');
            end
        end

        plot(sampleTimes, data.Sorting(1:length(sampleTimes)), '-or', 'LineWidth',3, 'MarkerSize',8)
        yline(mean(data.Sorting, 'all', 'omitmissing'), '--r', 'LineWidth',2)
        ylim([1.5, 4.2])
        yticks(2:4)
        ax_right(j) = gca;
        ax_right(j).YColor = 'red';
        if j == 6
            ylabel('σ_{G}')
        end
        if j == 1 || j == 3 || j == 5 || j == 7
            yticklabels('')
        end

        % Calculate and store tidal variation of each location
        tidalVariation_M(i, j) = std(data.Mean_mm, 'omitmissing');
        tidalVariation_S(i, j) = std(data.Sorting, 'omitmissing');
    
        text(sampleTimes(1)+hours(8.4), ax_right(j).YLim(2)*.9,...
            sampleIDs(j), 'FontSize',fontsize)
        % if j < 5
        %     annotation('textbox', [0.45, 0.8-j*.15, 0.1, 0.1], 'String',sampleIDs(j), 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
        % else
        %     annotation('textbox', [0.85, 1.4-j*.15, 0.1, 0.1], 'String',sampleIDs(j), 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
        % end
    end
    hold off

    % xticks(axs, sampleTimes)    
    
    xlim(axs, [sampleTimes(1)-hours(1), sampleTimes(1)+hours(10)])
    grid(axs, 'on')
    
    xticklabels(axs(1:8), '')

end

% Only 1 sample on 1 Oct
tidalVariation_M(3, :) = NaN;
tidalVariation_S(3, :) = NaN;

clearvars i j k ax_right ax_left axs dry sampleTimes waterLevels sampleHeight targetDate dailyData


%% Check for co-variation
% If a p-value is below 0.05, the corresponding correlation coefficient is statistically significant.
% If above 0.05, the result is not statistically significant, and the observed correlation could be due to random chance.

n = [1, 2, 4, 5];
tideVar_M = tidalVariation_M(n, :);
tideVar_S = tidalVariation_S(n, :);

% Preallocate for storing results
numLocations = size(tideVar_M, 2); % Number of columns
correlationCoeffs = zeros(1, numLocations);
pValues = zeros(1, numLocations);

% Compute correlation and p-values for each column pair
for i = 1:numLocations
    [correlationCoeffs(i), pValues(i)] = corr(tideVar_M(:, i), tideVar_S(:, i), 'Rows', 'complete');
end

% Display results
disp('Correlation coefficients for each location:');
disp(correlationCoeffs);
disp('P-values for each location:');
disp(pValues);

clearvars numLocations


%% Visualisation: tidal variation per location
f3 = figure(Position=[740, 1838, 1458, 455]);
tiledlayout(2, 5, 'TileSpacing','compact')

nexttile(1, [1, 3])
hold on
plot(xpos, tideVar_M', 'o-', 'LineWidth',2)
colororder(surveyColours([2, 3, 8, 9], :))
plot(xpos, mean(tideVar_M, 'omitmissing'), 'r', 'LineWidth',4)
% plot(xpos, median(tideVar_M, 'omitmissing'), ':r', 'LineWidth',4)
hold off

ylabel('σ_{tide,M_G} (mm)')

xlim([xpos(1)-1, xpos(end)+1])
xticks(xpos)
xticklabels('')
grid on

legend([cellstr(string(samplingDates(n), 'dd-MMM')); 'mean'], 'Location', 'northeastoutside')

nexttile(6, [1, 3])
hold on
plot(xpos, tideVar_S', 'o-', 'LineWidth',2)
plot(xpos, mean(tideVar_S, 'omitmissing'), 'r', 'LineWidth',4)
% plot(xpos, median(tideVar_S, 'omitmissing'), ':r', 'LineWidth',4)
hold off

ylabel('σ_{tide,σ_{G}}')

xlim([xpos(1)-1, xpos(end)+1])
xticks(xpos)
xticklabels(sampleIDs)
grid on

annotation('textbox', [0.49, 0.815, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.49, 0.375, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.895, 0.815, 0.1, 0.1], 'String','(c)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');

nexttile(4, [2,2])
bar(correlationCoeffs)
% xticks(xpos)
xticklabels(sampleIDs)
ylabel('r')
grid on


% clearvars tideVar_M tideVar_S n


%% Visualisation: tidal variation per location
n = [1, 2, 4, 5];
tideVar_M = tidalVariation_M(n, :);
tideVar_S = tidalVariation_S(n, :);

f3b = figure(Position=[740, 1813, 936, 480]);
tiledlayout(2, 1, 'TileSpacing','compact')

nexttile
hold on
plot(xpos, tideVar_M', 'o-', 'LineWidth',2)
colororder(surveyColours([2, 3, 8, 9], :))
plot(xpos, mean(tideVar_M, 'omitmissing'), 'r', 'LineWidth',4)
% plot(xpos, median(tideVar_M, 'omitmissing'), ':r', 'LineWidth',4)
hold off

ylabel('σ_{tide,M_G} (mm)')

xlim([xpos(1)-1, xpos(end)+1])
xticks(xpos)
xticklabels('')
grid on

legend([cellstr(string(samplingDates(n), 'dd-MMM')); 'mean'], 'Location', 'northeastoutside')

nexttile
hold on
plot(xpos, tideVar_S', 'o-', 'LineWidth',2)
plot(xpos, mean(tideVar_S, 'omitmissing'), 'r', 'LineWidth',4)
% plot(xpos, median(tideVar_S, 'omitmissing'), ':r', 'LineWidth',4)
hold off

ylabel('σ_{tide,σ_{G}}')

xlim([xpos(1)-1, xpos(end)+1])
xticks(xpos)
xticklabels(sampleIDs)
grid on

annotation('textbox', [0.78, 0.815, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.78, 0.375, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');

% clearvars tideVar_M tideVar_S n


%% Horizontal sampling error

% Emerged samples: 20 cm above water level
M_S1_20Sep = [0.945037; 1.22433];
S_S1_20Sep = [2.67655; 3.3442];

M_S1_28Sep = [0.759304; 0.999734];
S_S1_28Sep = [2.40133; 2.60146];

M_S1_07Oct = [0.89689; 1.48307];
S_S1_07Oct = [2.40958; 3.23303];

M_S2_07Oct = [0.646713; 0.731363];
S_S2_07Oct = [2.06916; 2.43918];

% Calculate mean stats
Mg_mean = mean([M_S1_20Sep, M_S1_28Sep, M_S1_07Oct, M_S2_07Oct]);
Mg_std = std([M_S1_20Sep, M_S1_28Sep, M_S1_07Oct, M_S2_07Oct]);

Sg_mean = mean([S_S1_20Sep, S_S1_28Sep, S_S1_07Oct, S_S2_07Oct]);
Sg_std = std([S_S1_20Sep, S_S1_28Sep, S_S1_07Oct, S_S2_07Oct]);

% Display bulk stats
% mean(Mg_mean([1,2,4]))
% mean(Mg_std([1,2,4]))
% 
% mean(Sg_mean([1,2,4]))
% mean(Sg_std([1,2,4]))