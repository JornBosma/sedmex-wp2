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

ZeefData = massRetained_20201202(1:end-1, :); % 45 mu sieve too fine (rounds off to zero)
ZeefDiam = [16000; sieveSizes_20201202(1:end-1)] ./ 1e6; % Convert to meters

sieveDiam = [flip([16000; sieveSizes_20201202] ./ 1e6)]'; % Convert to meters
sieveData = flip(massRetained_20201202);

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
f = figure('Position',[740, 957, 1719, 1336]);
tiledlayout(2, 1)

nexttile
for i = 1%:length(ZeefPerc)
    semilogx(ZeefDiam(1:end-1), ZeefPerc(:, i));
    hold on
    plot(Percen(:, i), Perc, '.');
    plot(GeoMeanStd(1, i), 0.5, 'o', ...
        2 ^ (PsiMeanStd(1, i) - PsiMeanStd(2, i)), 0.5, '>', ...
        2 ^ (PsiMeanStd(1, i) + PsiMeanStd(2, i)), 0.5, '<');
end
hold off
xlabel('Particle diameter (mm)');
ylabel('Probability');
title('Cumulative Sieve Curve with Percentiles and Mean Diameter');

nexttile
for i = 1%:length(ZeefPerc)
    yyaxis left
    semilogx(FracDiam, ZeefData(:, i));
    hold on
    plot(GeoMeanStd(1, i), 0.01, 'o', ...
        2 ^ (PsiMeanStd(1, i) - PsiMeanStd(2, i)), 0.01, '>', ...
        2 ^ (PsiMeanStd(1, i) + PsiMeanStd(2, i)), 0.01, '<');
    ylabel('Fraction')

    yyaxis right
    % Plotting the histogram
    histogram('BinEdges',sieveDiam, 'BinCounts',sieveData(:, i), 'FaceAlpha',.1);
    yticks([])
end
hold off
xlabel('Particle diameter (mm)');
title('Non-Cumulative Probability Distribution');
