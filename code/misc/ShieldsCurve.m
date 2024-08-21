%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;


%% Parameters
g = 9.81;      % gravitational acceleration [m/s^2]
nu = 1.3e-6;   % kinematic viscosity of water [m^2/s] (psu = 35 ppt, T = 12°C)
rhoW = 1025;   % density of seawater [kg/m3]
rhoS = 2650;   % density of sediment [kg/m3]

% d10 = 272e-6;  % 10th percentile grain diameter [m]
% d50 = 557e-6;  % median grain diameter [m]
% d90 = 1939e-6; % 90th percentile grain diameter [m]

% d10 = 451e-6;  % 10th percentile grain diameter [m]
% d50 = 1334e-6; % median grain diameter [m]
% d90 = 4105e-6; % 90th percentile grain diameter [m]

d10 = 277e-6;  % 10th percentile grain diameter [m]
d50 = 672e-6; % median grain diameter [m]
d90 = 2427e-6; % 90th percentile grain diameter [m]

dn = [d10, d50, d90];

dnVec = 1.4e-4:1e-4:1e-1;


%% Computations
[thetaCr, tauCr, dStar] = computeCBSS(dn, d50, rhoS, rhoW);
[thetaCrVec, tauCrVec, dStarVec] = computeCBSS(dnVec, d50, rhoS, rhoW);


%% Visualisation
figureRH;

yyaxis left
p1 = loglog(dStarVec, thetaCrVec, 'LineWidth',5);
yline(thetaCr, '--', {'d_{10}','d_{50}','d_{90}'},...
    'LineWidth',2, 'FontSize',fontsize, 'Color',p1.Color,...
    'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','right')
ylim([.01, 1])
ylabel('\theta_{cr}')

yyaxis right
p2 = loglog(dStarVec, tauCrVec, 'LineWidth',5);
yline(tauCr, '--', {'d_{10}','d_{50}','d_{90}'},...
    'LineWidth',2, 'FontSize',fontsize, 'Color',p2.Color,...
    'LabelVerticalAlignment','middle', 'LabelHorizontalAlignment','left')
ylabel('\tau_{cr} (Pa)')

xline(dStar, '--', {'d_{10}','d_{50}','d_{90}'},...
    'LineWidth',2, 'FontSize',fontsize, 'LabelHorizontalAlignment','center')
xlabel('D^*')

