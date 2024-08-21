%% Initialization
close all
clear
clc

[~, ~, ~, ~, ~] = sedmex_init;

dataPath = fullfile(filesep, 'Volumes', 'T7 Shield', 'DataDescriptor', 'hydrodynamics', 'ADV', filesep);

% Constants
g = 9.81;
rho_s = 2650;
rho_w = 1025;
D90 = 2138e-6;

% Defined starting time
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');

% Load ADV data
filename = [dataPath 'L4C1VEC' filesep 'tailored_L4C1VEC.nc'];
info = ncinfo(filename);

t = ncread(filename, 't');          % seconds since 2021-09-01
tStart = t0 + seconds(t(1));
tEnd = t0 + seconds(t(end));

u_z = ncread(filename, 'umag');     % flow velocity [m/s] at depth z 
z = ncread(filename, 'h');          % height above the bed [m]
H = ncread(filename, 'Hm0');        % significant wave height [m]
k = ncread(filename, 'k');          % wave number [m^-1]
% T = ncread(filename, 'Tmm10');      % wave period [s]
T = ncread(filename, 'Tp');         % peak wave period [s]
h = ncread(filename, 'd');          % water depth [m]
phi_c = ncread(filename, 'uang');   % current direction [°]
phi_w = ncread(filename, 'puvdir'); % wave propagation direction [°]
Urms = ncread(filename, 'u_ssm');   % rms total (u_ss+v_ss) orbital velocity [m/s] at depth z
Urms = sqrt(2) .* Urms;             % rms orbital velocity amplitude at depth z
% The sqrt(2) factor accounts for the conversion from RMS values of the
% individual velocity components to the combined amplitude of the
% orbital velocity. It ensures that the relationship between the peak
% values and the RMS values of the sinusoidal components is correctly
% represented when combining the two orthogonal components into a
% single RMS amplitude.

% Estimate the depth-averaged current velocity
k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)
[u_c, ~] = compute_DAV(u_z, z, k_sc, h);

% Compute the shear stress components
[tau_c1, tau_w1, tau_cw1] = compute_BSS_orbital(u_c, h, Urms, T, rho_w, phi_c, phi_w, D90, g);
[tau_c2, tau_w2, tau_cw2] = compute_BSS_linearTheory(u_c, H, k, T, h, rho_w, phi_c, phi_w, D90, g);

% Visualise results
figure
hold on
scatter(tau_w2, tau_w1)
plot([0, 4], [0, 4], 'r-')
hold off

xlabel('\tau_{w,linearTheory} (Pa)')
ylabel('\tau_{w,Urms} (Pa)')

xlim([0, 4])
ylim([0, 4])
axis square