%% Initialisation
[basePath, ~, ~, ~, ~] = sedmex_init;

locResults = [basePath 'results' filesep];
locOverleaf = '/Users/jorn/Library/CloudStorage/Dropbox/Apps/Overleaf/SEDMEX_WP2/figures/';

%% beachfaceGS
% exportgraphics(f1, [locResults 'BF_GS_mean.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'BF_GS_sorting.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'BF_GS_skewness.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'BF_GS_kurtosis.png'], 'Resolution',300)
% exportgraphics(f8, [locResults 'BF_GS_locations.png'], 'Resolution',300)

%% boundaryConditions
% exportgraphics(f1, [locResults 'conditions.png'], 'Resolution',300)

%% CrossShoreGS
% exportgraphics(f1, [locResults 'CS_profile.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'CS_GS_prof.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'CS_GS_time.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'CS_dGSdZ.png'], 'Resolution',300)

%% CrossShoreGS LW only
% exportgraphics(f1, [locResults 'CS_LW_profile.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'CS_LW_GS_prof.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'CS_LW_GS_time.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'CS_LW_dGSdZ.png'], 'Resolution',300)
% exportgraphics(f5, [locResults 'CS_LW_dMSdZ.png'], 'Resolution',300)

%% DEM
% exportgraphics(f1, [locResults 'DEM_locations.png'], 'Resolution',300)

%% flowRoses
% exportgraphics(f1, [locResults 'roseCurrent.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'roseWave.png'], 'Resolution',300)

%% LongShoreGS_2
% exportgraphics(f1, [locResults 'LS_GS_mean_0m.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'LS_GS_sorting_0m.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'LS_GS_mean_080m.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'LS_GS_sorting_080m.png'], 'Resolution',300)
exportgraphics(f5, [locResults 'LSGSD.png'], 'Resolution',300)

%% LongShoreMobility
% exportgraphics(f1, [locResults 'BSS_components.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'Shields_WvC.png'], 'Resolution',300)
% exportgraphics(f3, [locResults 'TS_tau.png'], 'Resolution',300)
% exportgraphics(f4, [locResults 'TS_Einstein.png'], 'Resolution',300)
% exportgraphics(f5, [locResults 'TS_BL_rate.png'], 'Resolution',300)
% exportgraphics(f6, [locResults 'TS_BL_rate_zeros.png'], 'Resolution',300)
% exportgraphics(f6, [locResults 'PDFs.png'], 'Resolution',300)

%% scraper
% exportgraphics(f1, [locResults 'scraper_mean.png'], 'Resolution',300)
% exportgraphics(f2, [locResults 'scraper_sorting.png'], 'Resolution',300)
