%% Initialisation
close all
clear
clc

[~, fontsize, ~, PHZ, SEDMEX] = sedmex_init;
% fontsize = 22; % ultra-wide screen

% Load DEMs
A = load('PHZ_2019_Q2','-mat'); % first survey
B = load('PHZ_2022_Q3','-mat');

% Load polygons
pgns = getPgons;


%% Computations #1
mask = inpolygon(A.DEM.X, A.DEM.Y, pgns.site(:, 1), pgns.site(:, 2));

A.DEM.X(~mask) = NaN;
A.DEM.Y(~mask) = NaN;
A.DEM.Z(~mask) = NaN;
B.DEM.X(~mask) = NaN;
B.DEM.Y(~mask) = NaN;
B.DEM.Z(~mask) = NaN;

% calculate DEM of Difference (DoD)
BminA = B.DEM.Z-A.DEM.Z;

mask = inpolygon(A.DEM.X, A.DEM.Y, pgns.harbour(:, 1), pgns.harbour(:, 2));
BminA(mask) = NaN;

C = A;
C.DEM.X(mask) = NaN;
C.DEM.Y(mask) = NaN;

D = B;
D.DEM.X(mask) = NaN;
D.DEM.Y(mask) = NaN;


%% Visualisation
% Original colormap data
originalColormap = colormap(brewermap(256, '-RdBu'));

% Extend the white part in the middle
numColors = size(originalColormap, 1);
extendedWhitePart = repmat([1, 1, 1], round(numColors * 0.2), 1); % Increase the white part

% Combine the original parts and the extended white part
customColormap = [
    originalColormap(1:round(numColors * 0.31), :); % First part of the original colormap
    extendedWhitePart; % Extended white part
    originalColormap(round(numColors * 0.7):end, :) % Last part of the original colormap
]; close all

f1 = figure('Position',[740, 957, 1719, 1336]);
surf(C.DEM.X, C.DEM.Y, BminA); hold on

shading('flat')
axis off equal
view(46.4, 90)

cb = colorbar;
cb.Location = 'northoutside';
set(cb, 'position', [.24 .64 .46 .015])
cb.Label.String = 'bed-level change (m)';
cb.TickLabels = {'-2', '-1.5', '-1', '-0.5', '0', '+0.5', '+1', '+1.5', '+2'};
% cb.FontSize = 26;

% Set color limits and colormap
clim([-2, 2])
colormap(customColormap)

% % North arrow
% Narrow(fontsize) % make arrow white

% % MHW & MLW contours
% contour(D.DEM.X, D.DEM.Y, D.DEM.Z, [PHZ.MeanLW, PHZ.MeanHW], '-k', 'ShowText', 'off', 'LineWidth', 1)
