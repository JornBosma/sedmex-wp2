%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];

% shoreNormal = 165; % shore-normal direction (cartesian convention) L1C1
shoreNormal = 135; % shore-normal direction (cartesian convention) L2Cx,L5C1
% shoreNormal = 122; % shore-normal direction (cartesian convention) L3C1,L4C1
% shoreNormal = 142; % shore-normal direction (cartesian convention) L6C1
shoreLine = shoreNormal - 90; % direction of the shoreline

g = 9.81; % gravitational acceleration


%% Correct measured height of control volume above bed
hab_measured = load('hab_C10.mat');
hab_measured = hab_measured.hab./100; % [cm] to [m]

% Iterate over each variable in the timetable
for varName = hab_measured.Properties.VariableNames
    % Get the variable data
    data = hab_measured.(varName{1});
    
    % Apply the condition: if value < 0.1, set it to 0.1 (to prevent
    % extremely high depth-average velocities)
    data(data < 0.12) = 0.12;
    
    % Update the variable in the timetable
    hab_measured.(varName{1}) = data;
end

clearvars varName data


%% Load ADV data
adv = 'L2C10VEC';
% adv = 'L2C5SONTEK1';
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)

% ADVpath = [dataPath 'ADV' filesep adv filesep 'tailored_' adv '_enc.nc'];
ADVpath = [dataPath 'ADV' filesep adv filesep 'tailored_' adv '.nc'];

t = ncread(ADVpath, 't'); % seconds since 2021-09-01 00:00:00
info = ncinfo(ADVpath);
% z = ncread(ADVpath, 'h');
h = ncread(ADVpath, 'd');
Uz = ncread(ADVpath, 'umag');
Ulong = ncread(ADVpath, 'ulm');
Ucros = ncread(ADVpath, 'ucm');
Udir = ncread(ADVpath, 'uang');
Hm0 = ncread(ADVpath, 'Hm0');
Hdir = ncread(ADVpath, 'puvdir');

% Interpolate the measurements to the new time vector
time = t0 + seconds(t);
habCVnew = retime(hab_measured, time, 'pchip');
z = habCVnew.(adv);

% Estimate the depth-averaged current velocity
[Ud, ~] = compute_DAV(Uz, z, k_sc, h);
[Udl, ~] = compute_DAV(Ulong, z, k_sc, h);
[Udc, ~] = compute_DAV(Ucros, z, k_sc, h);

TT = timetable(time, z, Uz, Ud, Udl, Udc, Udir, Hm0, Hdir);
TT = TT(90:end, :); % Udc shows outliers at the start


%% Visualisation
figure
plot(TT.time, TT.Udl, 'k'); hold on
plot(TT.time, TT.Udc, 'r')

figure
plot(t, Ulong, 'k'); hold on
plot(t, Ucros, 'r')


%% Computations
Hs_br = Hm0; % breaker height
theta_br = Hdir - shoreLine; % breaker wave angle relative to the shoreline

% Based on analysis of measured longshore current velocities at the Duck
% site (USA) and at the Egmond site (Netherlands) (Van Rijn, 2001)
K = 0.3; % empirical coefficient (ChatGPT: about 0.39 for sandy beaches)

% Longshore current velocity in mid of surf zone due to breaking waves
Vwave = K * sqrt(g * Hs_br) .* sind(2*theta_br);

% Model verification results based on field data the CERC formula yields
% results, which are slightly too large (factor 2) compared with measured
% values for storms but much too large (factor 5) for low wave conditions.


%% Visualisation
figure('Position',[743, 1668, 1708, 617]);
tiledlayout(2,1)

nexttile
plot(time, Vwave, 'LineWidth',2); hold on
plot(time, Udl, 'LineWidth',2)
ylabel('v_{long} (m s^{-1})')
legend('CERC','observed')

nexttile
yyaxis left
plot(time, Hs_br, 'LineWidth',2)
ylabel('H_{m0} (m)')

yyaxis right
plot(time, theta_br, '.', 'LineWidth',2)
ylabel('\theta (°)')


%% Check correlation
% Assuming Vwave and Udl are both time series of the same length

% Ensure both Vwave and Udl are column vectors
Vwave = Vwave(:);
Udl = Udl(:);

% Remove NaNs from both series
validIndices = ~isnan(Vwave) & ~isnan(Udl);

% Filter out NaN values
Vwave_clean = Vwave(validIndices);
Udl_clean = Udl(validIndices);

% Check if there are still valid data points after removing NaNs
if isempty(Vwave_clean) || isempty(Udl_clean)
    error('After removing NaNs, there are no valid data points to correlate.');
else
    % Calculate Pearson correlation coefficient
    correlation = corr(Vwave_clean, Udl_clean);
    fprintf('The Pearson correlation coefficient between Vwave and Udl (after handling NaNs) is: %.4f\n', correlation);
    
    % Calculate R² (Coefficient of Determination)
    % R² is the square of the Pearson correlation coefficient
    R2 = correlation^2;
    fprintf('The R² value between Vwave and Udl is: %.4f\n', R2);
    
    % Calculate degree of coherence using mscohere function
    % Set parameters for mscohere: window, noverlap, and nfft
    window = hamming(1024);  % Window length of 1024
    noverlap = 512;          % 50% overlap
    nfft = 1024;             % FFT points set to window size
    
    % Calculate degree of coherence using mscohere function
    [Cxy, f] = mscohere(Udl_clean, Vwave_clean, window, noverlap, nfft);
    
    % Plot the degree of coherence
    figure('Position',[743, 965, 1708, 617]);
    plot(f, Cxy);
    title('Coherence between Vwave and Udl');
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    
    % Optionally, calculate average coherence over all frequencies
    avg_coherence = mean(Cxy);
    fprintf('The average coherence between Vwave and Udl is: %.4f\n', avg_coherence);
end


%% Interpretation

% A Pearson correlation coefficient of -0.4207 falls in the moderate
% negative correlation range, indicating that there is a clear inverse
% relationship between the two variables, but it is not very strong.

% While your Pearson correlation coefficient of -0.4207 suggests a moderate
% negative correlation, the R² value of 0.18 reveals that the linear model
% is not particularly strong in explaining the relationship between the two
% series. In other words, while there is some linear relationship, the two
% time series have a lot of unexplained variance that might come from
% noise, nonlinearity, or other external influences.

% An average coherence of 0.23 suggests low to moderate
% correlation between Vwave and Udl across the frequency spectrum. This
% means that, on average, only about 23% of the power in Udl can be
% explained by Vwave at corresponding frequencies.