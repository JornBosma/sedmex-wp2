%% Initialisation
[basePath, ~, ~, ~, ~] = sedmex_init;

locResults = [basePath 'results' filesep];
locOverleaf = '/Users/jorn/Library/CloudStorage/Dropbox/Apps/Overleaf/SEDMEX_WP2/figures/';

%% studySite
% exportgraphics(f1, [locResults 'studySite.png'], 'Resolution',300)

%% exSample
% exportgraphics(f1, [locResults 'exSample.png'], 'Resolution',300)

%% beachfaceGS
% exportgraphics(f1, [locResults 'BF_GS_mean.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'BF_GS_sorting.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'BF_GS_skewness.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'BF_GS_kurtosis.png'], 'Resolution',300)
% exportgraphics(f8, [locResults 'BF_GS_locations.png'], 'Resolution',300)

%% boundaryConditions
% exportgraphics(f1, [locOverleaf 'boundCond.png'], 'Resolution',300)

%% CrossShoreGS
% exportgraphics(f1, [locResults 'CS_profile.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'CS_GS_prof.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'CS_GS_time.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'CS_dGSdZ.png'], 'Resolution',300)

%% CrossShoreGS_LW
exportgraphics(f1, [locResults 'CS_LW_GS_prof.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'CS_LW_GS_time.png'], 'Resolution',300)
% exportgraphics(f5, [locOverleaf 'CS_LW_dMSdz.png'], 'Resolution',300)

% exportgraphics(f4, [locResults 'CS_LW_dGSdZ.png'], 'Resolution',300)

%% CrossShoreGS_tide
% exportgraphics(f0a, [locResults 'CS_profile.png'], 'Resolution',300)

% exportgraphics(f1a, [locResults 'CS_GS_mean.png'], 'Resolution',300)
% exportgraphics(f1b, [locResults 'CS_GS_sort.png'], 'Resolution',300)

% for i = 1:5
%     exportgraphics(f2(i), [locOverleaf 'CS_GS_time_' mat2str(i) '.png'], 'Resolution',300)
% end

% exportgraphics(f3, [locOverleaf 'CS_tidal_variation.png'], 'Resolution',300)

%% DEM
% exportgraphics(f1, [locResults 'DEM_2022_Q1.png'], 'Resolution',300)

%% flowRoses
% exportgraphics(f1, [locResults 'roseCurrent.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'roseWave.png'], 'Resolution',300)

%% LongShoreGS_2
% exportgraphics(f6, [locResults 'LS_GSD.png'], 'Resolution',300)

% exportgraphics(f1, [locResults 'LS_GS_mean_0m.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'LS_GS_sorting_0m.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'LS_GS_mean_080m.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'LS_GS_sorting_080m.png'], 'Resolution',300)

%% LongShoreGS_3
% for i = 1:length(f2)
% for i = 6%[1, 3, 5, 6] % L6, L4, L2, L0
%     exportgraphics(f2(i), [locResults 'LS_GSD_' mat2str(i) '.png'], 'Resolution',300)
% end

%% LongShoreMobility
% exportgraphics(f1, [locResults 'BSS_components.png'], 'Resolution',300)
% exportgraphics(f2b, [locResults 'Shields_WvC.png'], 'Resolution',300)
% exportgraphics(f6d, [locOverleaf 'PDFs_MG.png'], 'Resolution',300)

% exportgraphics(f3a, [locResults 'TS_tau.png'], 'Resolution',300)
% exportgraphics(f4a, [locResults 'TS_Einstein.png'], 'Resolution',300)
% exportgraphics(f5a, [locResults 'TS_BL_rate.png'], 'Resolution',300)
% exportgraphics(f5c, [locResults 'TS_BL_rate_zeros.png'], 'Resolution',300)
% exportgraphics(f6a, [locResults 'PDFs_10.png'], 'Resolution',300)
% exportgraphics(f6b, [locResults 'PDFs_50.png'], 'Resolution',300)
% exportgraphics(f6c, [locResults 'PDFs_90.png'], 'Resolution',300)

%% scraper
% exportgraphics(f, [locOverleaf 'scraper.png'], 'Resolution',300)

% exportgraphics(f1, [locResults 'scraper_mean.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'scraper_sort.png'], 'Resolution',300)