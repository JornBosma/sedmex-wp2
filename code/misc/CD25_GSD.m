%% Initialisation
close all
clear
clc

[~, ~, cbf, ~, ~] = sedmex_init;

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];


%% GS_20201202
dataPath = [folderPath 'grainsizes' filesep 'GS_20201202.csv'];

opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20201202 = readtable(dataPath, opts);
sieveSizes_20201202 = fliplr(extractSieveNumbers(GS_20201202));
massRetained_20201202 = fliplr(GS_20201202{:, 17:end});

clearvars dataPath opts


%% GS_20211008/9
dataPath1 = [folderPath 'grainsizes' filesep 'GS_20211008.csv'];
dataPath2 = [folderPath 'grainsizes' filesep 'GS_20211009.csv'];

opts = detectImportOptions(dataPath1);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20211008 = readtable(dataPath1, opts);
GS_20211009 = readtable(dataPath2, opts);
GS_20211008 = [GS_20211009; GS_20211008];

sieveSizes_20211008 = fliplr(extractSieveNumbers(GS_20211008));
massRetained_20211008 = fliplr(GS_20211008{:, 17:end});

clearvars GS_20211009 dataPath1 dataPath2 opts


%% GS_20221026
dataPath = [folderPath 'grainsizes' filesep 'GS_20221026.csv'];

opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20221026 = readtable(dataPath, opts);
sieveSizes_20221026 = fliplr(extractSieveNumbers(GS_20221026));
massRetained_20221026 = fliplr(GS_20221026{:, 17:end});

clearvars dataPath opts


%% Merge data
tables = {GS_20201202, GS_20211008, GS_20221026}; % Collection of tables
sieveSizes = {sieveSizes_20201202, sieveSizes_20211008, sieveSizes_20221026}; % Corresponding sieve sizes
massRetained = {massRetained_20201202, massRetained_20211008, massRetained_20221026};


%% Visualisation
n = 2; % Select the table index
for i = 22%1:height(tables{n}) % Loop through each row of the selected table

    % Computations
    sieve_sizes = [sieveSizes{n}, 16000];
    bin_widths = diff(sieve_sizes);
    mass = massRetained{n}(i, :); % Extract mass retained values
    relativeMass = mass ./ bin_widths;
    cumulativeMass = cumsum(relativeMass); % Compute cumulative mass
    totalMass = cumulativeMass(end); % Total mass is the last cumulative value
    frequency = [relativeMass / totalMass * 100, 0]; % Compute frequency distribution

    % Visualisation
    f(i) = figure('Position', [740, 1873, 560, 420]); % Create new figure
    ax = axes; % Create new axes

    boundaries = sieve_sizes; % Define boundaries from sieve sizes
    h = stairs(boundaries, frequency, 'LineStyle', 'none'); hold on
    x = repelem(h.XData(2:end), 2);
    y = repelem(h.YData(1:end-1), 2);
    x(end) = []; y(1) = [];
    fill([x, fliplr(x)], [y, min(h.YData)*ones(size(y))], cbf.orange, 'LineWidth', 1)
    
    for j = 1:length(h.XData)
        line([h.XData(j), h.XData(j)], [0, h.YData(j)], 'Color', 'k', 'LineWidth', 1)
    end
    
    ax.XScale = 'log'; % Set x-axis to log scale
    box off % Remove the box around the plot
    axis off % Uncomment if you want to turn off the axis

end

clearvars n i j cumulativeMass totalMass frequency ax boundaries h x y sieve_sizes


%% Visualisation (combined)

idx = [49, 39, 26, 14, 4];
idx2 = [38, 27, 18, 14, 4];

f = gobjects(size(idx));
for i = 1:length(idx)
    f(i) = figure('Position', [740, 1873, 560, 420]); % Create new figure
    ax = axes; % Create new axes

    % Computations
    sieve_sizes = [sieveSizes{2}, 16000];
    bin_widths = diff(sieve_sizes);
    mass = massRetained{2}(idx(i), :); % Extract mass retained values
    relativeMass = mass ./ bin_widths;
    cumulativeMass = cumsum(relativeMass); % Compute cumulative mass
    totalMass = cumulativeMass(end); % Total mass is the last cumulative value
    frequency = [relativeMass / totalMass * 100, 0]; % Compute frequency distribution

    boundaries = sieve_sizes; % Define boundaries from sieve sizes
    h = stairs(boundaries, frequency, 'LineStyle', 'none', 'HandleVisibility','off'); hold on
    x = repelem(h.XData(2:end), 2);
    y = repelem(h.YData(1:end-1), 2);
    x(end) = []; y(1) = [];
    fill([x, fliplr(x)], [y, min(h.YData)*ones(size(y))], cbf.orange, 'LineWidth', 1)
    for j = 1:length(h.XData)
        line([h.XData(j), h.XData(j)], [0, h.YData(j)], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility','off')
    end


    % Computations
    sieve_sizes = [sieveSizes{3}, 16000];
    bin_widths = diff(sieve_sizes);
    mass = massRetained{3}(idx2(i), :); % Extract mass retained values
    relativeMass = mass ./ bin_widths;
    cumulativeMass = cumsum(relativeMass); % Compute cumulative mass
    totalMass = cumulativeMass(end); % Total mass is the last cumulative value
    frequency = [relativeMass / totalMass * 100, 0]; % Compute frequency distribution

    boundaries = sieve_sizes; % Define boundaries from sieve sizes
    h = stairs(boundaries, frequency, 'LineStyle', 'none', 'HandleVisibility','off'); hold on
    x = repelem(h.XData(2:end), 2);
    y = repelem(h.YData(1:end-1), 2);
    x(end) = []; y(1) = [];
    fill([x, fliplr(x)], [y, min(h.YData)*ones(size(y))], cbf.yellow, 'LineWidth', 1, 'FaceAlpha',.6)
    for j = 1:length(h.XData)
        line([h.XData(j), h.XData(j)], [0, h.YData(j)], 'Color', 'k', 'LineWidth', 1, 'HandleVisibility','off')
    end

    legend('2021-10-08', '2022-10-26')
    xlabel('Particle Diameter (μm)')
    ylabel('Class Weight (%)')
    ylim([0, 30])
    ax.XScale = 'log'; % Set x-axis to log scale
    box off % Remove the box around the plot
    % axis off % Uncomment if you want to turn off the axis

end

clearvars n i j cumulativeMass totalMass frequency ax boundaries h x y sieve_sizes


%% Save figures
% for i = 1:length(f)
%     exportgraphics(f(i), ['/Users/jorn/Library/CloudStorage/OneDrive-UniversiteitUtrecht/Events/CD/2025/figures/GSDs/', mat2str(i), '.png'], 'Resolution',300)
% end