%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, SEDMEX] = sedmex_init;
% fontsize = 22; % ultra-wide screen

% Load DEMs
DEMpath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'PHZ' filesep];
A = load([DEMpath, 'PHZ_2022_Q2.mat']);

% Load polygons
pgns = getPgons;

% Instrument locations
% OSSI = [SEDMEX.L6C2; SEDMEX.L5C2; SEDMEX.L4C3; SEDMEX.L2C9; SEDMEX.L1C2];
ADV = [SEDMEX.L6C1; SEDMEX.L5C1; SEDMEX.L4C1; SEDMEX.L3C1; SEDMEX.L2C5; SEDMEX.L1C1];
trackNames = {"L6", "L5", "L4", "L3", "L2", "L1"};
% sampleNames = {'L6', 'Tmb', 'L4', 'L3.5', 'L2', 'L0.5', 'Lgn'};
sampleNames = {'L6', 'Tmb', 'L4', 'L3.5', 'L2', 'L0', 'Lgn'};

clearvars DEMpath 


%% Longshore sampling locations
folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep];

% Load longshore sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210921.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_20210921 = readtable(dataPath, opts);
GS_20210921 = GS_20210921(:, 1:15);

% Load L2 sediment data
dataPath = [folderPath 'grainsizes' filesep 'GS_20210920.csv'];
opts = detectImportOptions(dataPath);
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');
GS_L2 = readtable(dataPath, opts);
GS_L2 = GS_L2(startsWith(GS_L2.Sample_Identity, ["L2C3","L2C5"]), :);
GS_L2 = GS_L2(:, 1:15);

% Merge tables
GS_LS = [GS_20210921; GS_L2];

% Create logical index of rows where 'zNAP_m' is greater than or equal to -0.2
logicalIndex = GS_LS.zNAP_m >= -0.2;

% Filter the table based on the logical index
GS_LS = GS_LS(logicalIndex, :);

% Keep only one timestamp of L2
GS_LS = GS_LS(1:7, :);

% Fix the longshore order
GS_LS = [GS_LS(1:4, :); GS_LS(7, :); GS_LS(5:6, :)];

clearvars GS_20210921 GS_L2


%% Computations #1
mask1 = inpolygon(A.DEM.X, A.DEM.Y, pgns.site(:, 1), pgns.site(:, 2));

A.DEM.X(~mask1) = NaN;
A.DEM.Y(~mask1) = NaN;
A.DEM.Z(~mask1) = NaN;

mask2 = inpolygon(A.DEM.X, A.DEM.Y, pgns.harbour(:, 1), pgns.harbour(:, 2));

B = A;
B.DEM.X(mask2) = NaN;
B.DEM.Y(mask2) = NaN;
[B.DEM.Z, winsize] = smoothdata(B.DEM.Z, 'gaussian', 100);

clearvars mask1 mask2


%% Visualisation: DEM
f1 = figure('Position',[740, 957, 1719, 1336]);
surf(A.DEM.X, A.DEM.Y, A.DEM.Z); hold on

shading flat
ax = gca; ax.SortMethod = 'childorder';
axis off equal
view(50.7, 90)

cb = colorbar;
cb.Location = 'northoutside';
set(cb, 'position', [.24 .64 .46 .015])
cb.Label.String = 'bed level (NAP+m)';
cb.FontSize = fontsize;
clim([-8, 8])
colormap(sandyToMarineBlue(256, true));

% Add black background
ax.SortMethod = 'childorder';
patch(ax, pgns.site(:,1), pgns.site(:,2), 'k', 'EdgeColor','k', 'LineWidth',3)
h = get(ax,'Children');
set(ax,'Children',[h(2) h(1)])

% Add relief shadow
light

% Add north arrow
Narrow(fontsize)

% Add relevant contours
[~, c1] = contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [SEDMEX.MeanLW, SEDMEX.MeanHW], '-k', 'ShowText','off', 'LineWidth',1);
[~, c2] = contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [SEDMEX.MinWL, SEDMEX.MaxWL], ':r', 'ShowText','off', 'LineWidth',2);

% Instrument locations
s1 = scatter(ADV(:,1), ADV(:,2), 150, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor','r');
% text(ADV([1:2 5],1)+40, ADV([1:2 5],2)+20, trackNames([1:2 5]), 'FontSize',fontsize*.8)
% text(ADV(3:4,1)+70, ADV(3:4,2), trackNames(3:4), 'FontSize',fontsize*.8)
text(ADV(:,1)+60, ADV(:,2)-50, trackNames, 'FontSize',fontsize*.8, 'FontWeight','bold')

% scatter(SEDMEX.L3C1(1), SEDMEX.L3C1(2), 150, 'filled', 'MarkerEdgeColor','k', 'MarkerFaceColor','m')
% text(SEDMEX.L3C1(1)+70, SEDMEX.L3C1(2)-80, 'L3', 'FontSize',fontsize*.8)

% Sampling locations
s2 = scatter(GS_LS.xRD_m, GS_LS.yRD_m, 100, 'filled', 'MarkerEdgeColor','k',...
    'MarkerFaceColor','y', 'Marker','v');
text(GS_LS.xRD_m-120, GS_LS.yRD_m-60, sampleNames, 'FontSize',fontsize*.8)

% Legend
% lgnd = legend([c1, s1, s2], {'beach face', 'instruments', 'samples'}, 'Position',[0.1044, 0.4664, 0.1093, 0.1321], 'FontSize',fontsize);
lgnd = legend([c1, c2, s1, s2], {'beach face', 'max tidal range', 'instruments', 'samples'}, 'Position',[0.105, 0.4664, 0.1093, 0.1321], 'FontSize',fontsize);

% Add relevant contours
% [~, c1] = contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [SEDMEX.MeanHW, SEDMEX.MeanHW], '-g', 'ShowText','off', 'LineWidth',2);
% [~, c2] = contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [SEDMEX.MeanLW, SEDMEX.MeanLW], '-k', 'ShowText','off', 'LineWidth',2);
% [~, c3] = contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [-0.75, -0.75], '-m', 'ShowText','off', 'LineWidth',2);
% [~, c4] = contour(B.DEM.X, B.DEM.Y, B.DEM.Z, [-1.4, -1.4], '-r', 'ShowText','off', 'LineWidth',2);

% Legend
% lgnd = legend([c1,c2,c3,c4], {num2str(SEDMEX.MeanHW, 2), num2str(SEDMEX.MeanLW, 2), '-0.75', '-1.40'},...
%     'Position',[0.1313, 0.4876, 0.0779, 0.1247]);
