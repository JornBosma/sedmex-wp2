function [basePath, fontsize, cbf, PHZ, SEDMEX] = sedmex_init

% define global variables
% global basePath

% set base path. All subdirectories are branches from the base path.
basePath = strrep(which('sedmex_init'), 'sedmex_init.m', '');

% add paths
addpath(basePath);
addpath(genpath([basePath 'code']))
addpath(genpath([basePath 'dataRaw']))
addpath(genpath([basePath 'data']))
addpath(genpath([basePath 'results']))

% set default settings
fontsize = 30;
set(groot, 'DefaultAxesFontSize', fontsize)
% set(groot, 'DefaultTextInterpreter', 'latex')
% set(groot, 'DefaultAxesTickLabelInterpreter', 'latex')
% set(groot, 'DefaultLegendInterpreter', 'latex')
% set(groot, 'DefaultLegendLocation', 'northwest')
set(groot, 'DefaultLegendBox', 'on')
set(groot, 'DefaultAxesBox', 'on')

set(groot, 'defaultUicontrolFontName', 'Arial')
set(groot, 'defaultUitableFontName', 'Arial')
set(groot, 'defaultAxesFontName', 'Arial')
set(groot, 'defaultUipanelFontName', 'Arial')
set(groot, 'defaultTextFontName', 'Arial')

% colourblind-friendly colour palette
cbf.orange = [230/255, 159/255, 0];
cbf.skyblue = [86/255, 180/255, 233/255];
cbf.bluegreen = [0, 158/255, 115/255];
cbf.yellow = [240/255, 228/255, 66/255];
cbf.blue = [0/255, 114/255, 178/255];
cbf.vermilion = [213/255, 94/255, 0/255];
cbf.redpurp = [204/255, 121/255, 167/255];
% not colourblind-friendly
cbf.qual12 = ["#a6cee3"; ...
        "#1f78b4"; ...
        "#b2df8a"; ...
        "#33a02c"; ...
        "#fb9a99"; ...
        "#e31a1c"; ...
        "#fdbf6f"; ...
        "#ff7f00"; ...
        "#cab2d6"; ...
        "#6a3d9a"; ...
        "#ffff99"; ...
        "#b15928"];
cbf.custom12 = ...
   [1.000000000000000   0.500000000000000   0.800000000000000; ...
    0.992156862745098	0.749019607843137	0.435294117647059; ...
    0.000000000000000   1.000000000000000   1.000000000000000; ...
    0.792156862745098	0.698039215686275	0.839215686274510; ...
    1.000000000000000   0.498039215686275	0.000000000000000; ...
    0.698039215686275	0.874509803921569	0.541176470588235; ...
    0.900000000000000	0.900000000000000	0.000000000000000; ...
    0.890196078431373	0.101960784313725	0.109803921568627; ...
    0.415686274509804	0.239215686274510	0.603921568627451; ...
    0.694117647058824	0.349019607843137	0.156862745098039; ...
    0.121568627450980	0.470588235294118	0.705882352941177; ...
    0.200000000000000	0.627450980392157	0.172549019607843];
cbf.grey = [0.5, 0.5, 0.5];
% cmap = crameri('batlowK');
% % cmap = brewermap(256, 'RdYlBu');
% numColors = size(cmap, 1);
% indices = round(linspace(1, numColors, 6));
% cbf.six = cmap(indices, :);
% cbf.six = brewermap(9, 'Set1');
% cbf.six = [
%     0.0000, 0.4470, 0.7410;  % Blue
%     0.8500, 0.3250, 0.0980;  % Orange
%     0.9290, 0.6940, 0.1250;  % Yellow
%     0.4940, 0.1840, 0.5560;  % Purple
%     0.4660, 0.6740, 0.1880;  % Green
%     0.3010, 0.7450, 0.9330;  % Cyan
% ];
cbf.six = [
    0.0000, 0.4470, 0.7410;  % Blue
    0.8500, 0.3250, 0.0980;  % Orange
    1.0000, 0.8276, 0.0000;  % Yellow
    0.4940, 0.1840, 0.5560;  % Purple
    0.4660, 0.6740, 0.1880;  % Green
    0.1000, 0.9000, 1.0000;  % Cyan
];

% tidal datums wrt NAP @ PHZ (waterinfo)
% PHZ.NAP = 0; % local reference datum [NAP+m]
PHZ.HWS = 0.81; % typical high water spring level [NAP+m] (visually determined from waterinfo)
PHZ.LWS = -1.07; % typical low water spring level [NAP+m] (visually determined from waterinfo)
PHZ.HWlim = 1.95; % upper limit normal water-level range [NAP+m]
PHZ.LWlim = -1.17;  % upper limit normal water-level range [NAP+m]

% tidal datums during monitoring period (2019 - 2022)
PHZ.MeanHW = 0.664; % mean high water [NAP+m]
PHZ.MeanLW = -0.703; % mean low water [NAP+m]
PHZ.MeanSL = 0.0833; % mean water level [NAP+m]
PHZ.MaxWL = 2.41; % maximum water level [NAP+m]
PHZ.MinWL = -1.81; % minimum water level [NAP+m]
PHZ.MeanTR = 1.37; % mean tidal range [NAP+m]
PHZ.MaxTR = 2.49; % maximum tidal range [NAP+m]
PHZ.MinTR = 0.52; % minimum tidal range [NAP+m]

% tidal datums (slotgemiddelden 2011): HHNK & Witteveen+Bos, 2016
PHZ.DHW = 2.95; % decennial (1/10y) high water level [NAP+m]
PHZ.BHW = 2.4; % biennial (1/2y) high water level [NAP+m]
PHZ.AHW = 2.25; % annual (1/1y) high water level [NAP+m]
% PHZ.MHW = 0.64; % mean high water level [NAP+m]
% PHZ.MSL = 0.04; % mean sea level [NAP+m]
% PHZ.MLW = -0.69; % mean low water level [NAP+m]
PHZ.LAT = -1.17; % lowest astronomical tide [NAP+m]

% axis limits and ticks
PHZ.xLim = [114800, 118050]; % RD-2008
PHZ.yLim = [557750, 560650];
PHZ.xLimWGS = [53.005205, 53.031496]; % WGS-84
PHZ.yLimWGS = [4.788338, 4.836422];
PHZ.xTick = 114000:1e3:118000;
PHZ.yTick = 558000:1e3:561000;
PHZ.xLimHook = [116950, 117650]; % spit hook zoom-in
PHZ.yLimHook = [559760, 560400];

% tidal datums during SEDMEX: Woerdman et al., 2022
% SEDMEX.MHW = 0.68; % mean high water [NAP+m]
% SEDMEX.MLW = -0.57; % mean low water [NAP+m]
% SEDMEX.MWL = 0.16; % mean water level [NAP+m]
% SEDMEX.MaxWL = 1.34; % maximum water level [NAP+m]
% SEDMEX.MinWL = -1.07; % minimum water level [NAP+m]
% SEDMEX.MSTR = 1.77; % mean spring tidal range [NAP+m]
% SEDMEX.MNTR = 0.92; % mean neap tidal range [NAP+m]
% SEDMEX.MTR = 1.25; % mean tidal range [NAP+m]
% SEDMEX.MaxTR = 1.82; % maximum tidal range [NAP+m]
% SEDMEX.MinTR = 0.55; % minimum tidal range [NAP+m]

% tidal datums during SEDMEX
SEDMEX.MeanHW = 0.684; % mean high water [NAP+m]
SEDMEX.MeanLW = -0.609; % mean low water [NAP+m]
SEDMEX.MeanSL = 0.139; % mean water level [NAP+m]
SEDMEX.MaxWL = 1.35; % maximum water level [NAP+m]
SEDMEX.MinWL = -1.06; % minimum water level [NAP+m]
SEDMEX.MeanTR = 1.29; % mean tidal range [NAP+m]
SEDMEX.MaxTR = 1.85; % maximum tidal range [NAP+m]
SEDMEX.MinTR = 0.53; % minimum tidal range [NAP+m]

% SEDMEX instrument coordinates
SEDMEX.L2C1 = [117158, 559855];  % L2C1KELLER
SEDMEX.L2C2 = [117193, 559822];  % L2C2VEC
SEDMEX.L2C4 = [117197, 559818];  % L2C4VEC
SEDMEX.L2C10 = [117235, 559781]; % L2C10VEC

SEDMEX.L1C1 = [117421, 560054];  % L1C1VEC
SEDMEX.L2C5 = [117198, 559815];  % L2C5SONTEK1
SEDMEX.L3C1 = [116839, 559536];  % L3C1VEC
SEDMEX.L4C1 = [116103, 558945];  % L4C1VEC
SEDMEX.L5C1 = [115670, 558604];  % L5C1VEC
SEDMEX.L6C1 = [115402, 558225];  % L6C1VEC

SEDMEX.L1C2 = [117445, 560045];  % L1C2OSSI
SEDMEX.L2C9 = [117222, 559793];  % L2C9OSSI
SEDMEX.L4C3 = [116125, 558917];  % L4C3OSSI
SEDMEX.L5C2 = [115716, 558560];  % L5C2OSSI
SEDMEX.L6C2 = [115470, 558176];  % L6C2OSSI

% ready
return
