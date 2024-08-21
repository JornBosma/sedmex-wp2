% Data sieve analyses conversion to cumulative sieve curves and percentiles
% Maarten Kleinhans August 2007

%% Initialisation
close all
clear
clc

[~, ~, cbf, ~, ~] = sedmex_init;

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];

%% Load input files
% ZeefData = load('lapeyne.txt');
% ZeefDiam = load('diamlapeyne.txt') ./ 1000; % Convert to meters

dataPath = [folderPath 'grainsizes' filesep 'GS_20201202.csv'];

opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20201202 = readtable(dataPath, opts);
sieveSizes_20201202 = [extractSieveNumbers(GS_20201202)]';
massRetained_20201202 = [GS_20201202{:, 17:end}]';

ZeefData = massRetained_20201202(1:end-1, :); % NOTE: 45 mu sieve too fine (rounds off to zero)
ZeefDiam = [16000; sieveSizes_20201202(1:end-1)] ./ 1e6; % Convert to meters
sieveDiam = [flip([16000; sieveSizes_20201202] ./ 1e6)]'; % Convert to meters
sieveData = flip(massRetained_20201202);

% % Exclude largest sieve size
% ZeefDiam = ZeefDiam(2:end);
% ZeefData = ZeefData(2:end,:);
% sieveDiam = sieveDiam(1:end-1);
% sieveData = sieveData(1:end-1,:);

% bin_widths = diff(sieveDiam);
% relSieveData = sieveData ./ bin_widths';

% Get dimensions of the data
[M, N] = size(ZeefData);

% Calculate fractional diameters
FracDiam = exp((log(ZeefDiam(1:end-1)) + log(ZeefDiam(2:end))) / 2);

% Correct negative values and normalize data
ZeefData(ZeefData < 0) = 0;
ZeefData = ZeefData ./ sum(ZeefData);

% Calculate cumulative probabilities (ZeefPerc)
ZeefPerc = flipud(cumsum(flipud(ZeefData)));

% Percentiles to calculate
Perc = [0.1, 0.16, 0.5, 0.84, 0.9, 0.95]';
Percen = NaN(length(Perc), N);

% Calculate percentiles
for dat = 1:N
    for per = 1:length(Perc)
        i = find(ZeefPerc(:, dat) > Perc(per), 1, 'last');
        if isempty(i) || i == M
            continue;
        end
        top = i;
        bot = i + 1;
        Percen(per, dat) = exp(log(ZeefDiam(top)) - ...
            (ZeefPerc(top, dat) - Perc(per)) / (ZeefPerc(top, dat) - ZeefPerc(bot, dat)) * ...
            (log(ZeefDiam(top)) - log(ZeefDiam(bot))));
    end
end

% Calculate geometric mean diameter and standard deviation
GeoMeanStd = NaN(2, N);
PsiMeanStd = NaN(2, N);
psi = log(FracDiam) / log(2);

for dat = 1:N
    PsiMeanStd(1, dat) = nansum(psi .* ZeefData(:, dat));
    GeoMeanStd(1, dat) = 2 ^ PsiMeanStd(1, dat);
    PsiMeanStd(2, dat) = sqrt(nansum(((psi - PsiMeanStd(1, dat)) .^ 2) .* ZeefData(:, dat)));
    GeoMeanStd(2, dat) = 2 ^ PsiMeanStd(2, dat);
end

% Plot figures
f = figure('Position',[740, 1665, 1719, 628]);
tiledlayout(1, 2)

nexttile
for i = 15%1:length(ZeefPerc)
    semilogx(ZeefDiam(1:end-1), ZeefPerc(:, i), '-o', 'LineWidth',3);
    hold on
    plot(Percen(:, i), Perc, '.', 'LineWidth',3);
    plot(GeoMeanStd(1, i), 0.5, 'og', ...
        2 ^ (PsiMeanStd(1, i) - PsiMeanStd(2, i)), 0.5, '>g', ...
        2 ^ (PsiMeanStd(1, i) + PsiMeanStd(2, i)), 0.5, '<g', 'LineWidth',3);
end
hold off

yline(0.1, '--')
yline(0.5, '--')
yline(0.9, '--')
xline(Percen(1,i), ':', 'LineWidth',3)
xline(Percen(3,i), '-', 'LineWidth',3)
xline(Percen(5,i), ':', 'LineWidth',3)

xlabel('Particle diameter (mm)');
ylabel('Probability');
title('CDF');

nexttile
for i = 15%1:length(ZeefPerc)
    yyaxis left
    semilogx(FracDiam, ZeefData(:, i), 'o-', 'LineWidth',2);
    hold on
    plot(GeoMeanStd(1, i), 0.01, 'og', ...
        2 ^ (PsiMeanStd(1, i) - PsiMeanStd(2, i)), 0.01, '>g', ...
        2 ^ (PsiMeanStd(1, i) + PsiMeanStd(2, i)), 0.01, '<g', 'LineWidth',3);
    ylabel('Fraction')

    yyaxis right
    % Plotting the histogram
    histogram('BinEdges',sieveDiam, 'BinCounts',sieveData(:, i), 'FaceAlpha',.1, 'LineStyle','none');
    yticks([])
end
hold off

xline(Percen(1,i), ':', 'LineWidth',3)
xline(Percen(3,i), '-', 'LineWidth',3)
xline(Percen(5,i), ':', 'LineWidth',3)

xlabel('Particle diameter (m)');
title('PDF');


%% Store figure
exportgraphics(f, '/Users/jorn/Downloads/2020-12-02_T05S02.png', 'Resolution',300)