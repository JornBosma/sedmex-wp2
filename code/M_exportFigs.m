%% Initialisation
[basePath, ~, ~, ~, ~] = sedmex_init;

locResults = [basePath 'results' filesep];
locOverleaf = '/Users/jorn/Library/CloudStorage/Dropbox/Apps/Overleaf/SEDMEX_WP2/figures/';

%% A_studySite
exportgraphics(f1, [locOverleaf 'studySite.png'], 'Resolution',300) % Fig. 1

%% B_exSample
exportgraphics(f1, [locResults 'exSample.png'], 'Resolution',300) % Fig. 2

%% C_DEM
exportgraphics(f1, [locResults 'DEM_2022_Q1.png'], 'Resolution',300) % Fig. 3a

%% D_flowRoses
exportgraphics(f1, [locResults 'roseCurrent.png'], 'Resolution',300) % Fig. 3b
exportgraphics(f2, [locResults 'roseWave.png'], 'Resolution',300) % Fig. 3c

%% E_noncumGSD
exportgraphics(f2, [locResults 'noncumGSD.png'], 'Resolution',300) % Fig. 3d

%% F_boundaryConditions
exportgraphics(f1, [locOverleaf 'boundCond.png'], 'Resolution',300) % Fig. 4

%% G_scraper
exportgraphics(f, [locResults 'scraper.png'], 'Resolution',300) % Fig. 5

%% H_LongShoreGS
exportgraphics(f6, [locResults 'LS_GSD.png'], 'Resolution',300) % Fig. 6

%% I_ShieldsDiagram
exportgraphics(f1, [locOverleaf 'BSS_timeseries.png'], 'Resolution',300) % Fig. 7
exportgraphics(f2b, [locOverleaf 'Shields_WvC.png'], 'Resolution',300) % Fig. 8

%% J_bedloadTransport
exportgraphics(f6f, [locOverleaf 'PDFs_FC.png'], 'Resolution',300) % Fig. 9

%% K_CrossShoreGS_LW
exportgraphics(f1, [locResults 'CS_profile.png'], 'Resolution',300) % Fig. 10a
exportgraphics(f3, [locResults 'CS_LW_GS_time.png'], 'Resolution',300) % Fig. 10b
exportgraphics(f2, [locResults 'CS_LW_GS_prof.png'], 'Resolution',300) % Fig. 10c
exportgraphics(f5, [locOverleaf 'CS_LW_dMSdz.png'], 'Resolution',300) % Fig. 12

%% L_CrossShoreGS_tide
exportgraphics(f3b, [locOverleaf 'CS_tidal_variation.png'], 'Resolution',300) % Fig. 11

% Fig. A1
for i = 1:5
    exportgraphics(f2(i), [locOverleaf 'CS_GS_time_' mat2str(i) '.png'], 'Resolution',300)
end


