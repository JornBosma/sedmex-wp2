%% Initialization
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

dataPath = fullfile(filesep, 'Volumes', 'T7 Shield', 'DataDescriptor', 'hydrodynamics', 'ADV', filesep);

instruments = ["L1C1VEC", "L2C10VEC", "L3C1VEC", "L4C1VEC", "L5C1VEC", "L6C1VEC"];
sampleLocs = {'L0', 'L2', 'L3.5', 'L4', 'Tmb', 'L6'};
instruLocs = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'};
GS_fractions = {'Mg', 'd10', 'd50', 'd90', 'Fine', 'Coarse'};
GS_headers = {'M_G', 'D_{10}', 'D_{50}', 'D_{90}', '0.25', '2.00'};

% newcolors = [cbf.vermilion; cbf.blue; cbf.bluegreen; cbf.yellow; cbf.redpurp; cbf.skyblue; cbf.orange];
newcolors = cbf.six([2, 1, 3:5, 7:end], :);
newsymbols = {'o', 's', '^', 'd', 'v', 'p'};

g = 9.81;     
rho_s = 2650;  
rho_w = 1025; 

stormy_period = [datetime('27-Sep-2021 12:00'), datetime('05-Oct-2021 18:00')];
storm_period = [datetime('30-Sep-2021 16:00'), datetime('03-Oct-2021 03:00')];


%% Computations: Critical Bed-Shear Stress & Critical Shields

% Sample data [L0, L2, L3.5, L4, Tmb, L6] averaged over 21/09 and 28/09

% NAP +0.0m [mm]
Mg_a = [0.7266; 0.8357; 0.7682; 0.6992; 0.6538; 0.6608];
d10_a = [0.2776; 0.2802; 0.3102; 0.2799; 0.3031; 0.2256];
d50_a = [0.6714; 0.8127; 0.7671; 0.5428; 0.6009; 0.5739];
d90_a = [2.3978; 2.9071; 2.1240; 2.6873; 1.7590; 3.0086];
fg_a = [0.1203; 0.1828; 0.1066; 0.1546; 0.0630; 0.1383];

% NAP -0.8m [mm]
Mg_b = [0.4290; 0.9000; 1.3128; 1.4467; 0.3249; 0.2641];
d10_b = [0.1580; 0.3006; 0.5094; 0.4683; 0.2097; 0.1817];
d50_b = [0.3083; 0.8962; 1.2409; 1.3588; 0.3305; 0.2670];
d90_b = [3.4175; 3.2338; 3.8910; 5.0204; 0.4951; 0.4003];
fg_b = [0.1209; 0.2079; 0.3062; 0.3336; 0.0129; 0.0084];

% Take the average over both isobaths
Mg = mean([Mg_a, Mg_b], 2);
d10 = mean([d10_a, d10_b], 2);
d50 = mean([d50_a, d50_b], 2);
d90 = mean([d90_a, d90_b], 2);
fg = mean([fg_a, fg_b], 2);

% Convert mm to m
Mg = Mg*1e-3;
d10 = d10*1e-3;
d50 = d50*1e-3;
d90 = d90*1e-3;
Fine = repmat(0.25*1e-3, length(Mg), 1);
Coarse = repmat(2*1e-3, length(Mg), 1);

GS_stats = table(Mg, d10, d50, d90, Fine, Coarse, fg, 'RowNames', flipud(sampleLocs));

% Initialize arrays for critical Shields and bed-shear stress computations
crit_methods = {'Soulsby', 'Egiazaroff', 'McCarron'};
num_methods = length(crit_methods);
num_percentiles = width(GS_stats)-1;
[theta_cr, tau_cr] = deal(NaN(height(GS_stats), num_percentiles, num_methods));

% Perform computations using critical bed-shear stress methods
for i = 1:num_percentiles
    for j = 1:num_methods
        method = crit_methods{j};
        switch method
            case 'Soulsby'
                [theta_cr(:, i, j), tau_cr(:, i, j)] = crit_Soulsby(GS_stats.(i), GS_stats.d50, rho_s, rho_w, g);
            case 'Egiazaroff'
                [theta_cr(:, i, j), tau_cr(:, i, j)] = crit_Egiazaroff(GS_stats.(i), GS_stats.d50, rho_s, rho_w, g);
            case 'McCarron'
                [theta_cr(:, i, j), tau_cr(:, i, j)] = crit_McCarron(GS_stats.(i), GS_stats.d50, GS_stats.fg, rho_s, rho_w, g);
        end
    end
end

% Arrays of field names and corresponding variables
prefixes = {'theta_cr', 'tau_cr'};
fractions = {'Mg', '10', '50', '90', 'Fine', 'Coarse'};

% Loop over prefixes, methods, and percentiles to assign values
for p = 1:length(prefixes)
    for j = 1:num_methods
        method = crit_methods{j};
        for perc = 1:length(fractions)
            fieldName = sprintf('%s%s_%s', prefixes{p}, fractions{perc}, method);
            GS_stats.(fieldName) = eval(sprintf('%s(:, %d, %d)', prefixes{p}, perc, j));
        end
    end
end

clearvars d10 d50 d90 fg theta_cr tau_cr crit_methods prefixes percentiles i j p perc num_methods num_percentiles method fieldName Mg Mg_a Mg_b d10_a d10_b d50_a d50_b d90_a d90_b fg_a fg_b Fine Coarse


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


%% Computations: bed-shear stress
% GS_stats(2,:) = [];  % Exclude L2

% Defined starting time
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');

% Create a datetime array for the time axis
dateStart = datetime(2021, 9, 10);
dateEnd = datetime(2021, 10, 19);
dt = minutes(10);
timeAxis = dateStart:dt:dateEnd;
timeAxis = timeAxis';

% Interpolate the measurements to the new time vector
habCVnew = retime(hab_measured, timeAxis, 'pchip');

% Preallocate the data matrix with NaNs
n = length(timeAxis);
emptyData = nan(n, length(instruments));

% List of variable suffixes
suffixes = {'tau_c', 'tau_w', 'tau_cw', ...
    'theta_c10', 'theta_w10', 'theta_cw10', 'phi_b10', 'q_b10', 'qL_b10_net', 'qC_b10_net', ...
    'theta_c50', 'theta_w50', 'theta_cw50', 'phi_b50', 'q_b50', 'qL_b50_net', 'qC_b50_net', ...
    'theta_c90', 'theta_w90', 'theta_cw90', 'phi_b90', 'q_b90', 'qL_b90_net', 'qC_b90_net', ...
    'theta_cMg', 'theta_wMg', 'theta_cwMg', 'phi_bMg', 'q_bMg', 'qL_bMg_net', 'qC_bMg_net', ...
    'theta_cFine', 'theta_wFine', 'theta_cwFine', 'phi_bFine', 'q_bFine', 'qL_bFine_net', 'qC_bFine_net', ...
    'theta_cCoarse', 'theta_wCoarse', 'theta_cwCoarse', 'phi_bCoarse', 'q_bCoarse', 'qL_bCoarse_net', 'qC_bCoarse_net', ...
    'long_sign', 'long_frac','cross_sign', 'cross_frac'};

% Create a timetable for each suffix and assign variable names
for i = 1:length(suffixes)
    suffix = suffixes{i};
    timetableName = ['TT_', suffix];
    eval([timetableName ' = array2timetable(emptyData, ''RowTimes'', timeAxis);']);
    eval([timetableName '.Properties.VariableNames = instruments;']);
end

for i = 1:length(instruments)
filename = [dataPath char(instruments(i)) filesep 'tailored_' char(instruments(i)) '.nc'];
    
    % Load ADV data
    info = ncinfo(filename);
    
    t = ncread(filename, 't');          % seconds since 2021-09-01
    time = t0 + seconds(t);
    tStart = time(1);                   % instrument-specific starting time
    tEnd = time(end);                   % instrument-specific end time

    % Interpolate the measurements to the new time vector
    habCVnew = retime(hab_measured, time, 'pchip');
    z_new = habCVnew.(i);

    uLong_z = ncread(filename, 'ulm');  % longshore flow velocity [m/s] at depth z
    uLong_sign = sign(uLong_z);         % direction sign of longshore current
    % Positive onshore direction at this cross-section is 135◦
    % counterclockwise from east, and positive alongshore direction is
    % 225◦ counterclockwise from east.
    uLong_sign = -uLong_sign;           % reverse the sign so that NE (flood) is positive
    
    uCross_z = ncread(filename, 'ucm');  % cross-shore flow velocity [m/s] at depth z
    uCross_sign = sign(uCross_z);        % direction sign of cross-shore current

    u_z = ncread(filename, 'umag');     % flow velocity [m/s] at depth z
    % z = ncread(filename, 'h');          % height above the bed [m]
    H = ncread(filename, 'Hm0');        % significant wave height [m]
    k = ncread(filename, 'k');          % wave number [m^-1]
    % T = ncread(filename, 'Tmm10');      % wave period [s]
    T = ncread(filename, 'Tp');         % peak wave period [s]
    h = ncread(filename, 'd');          % water depth [m]
    phi_c = ncread(filename, 'uang');   % current direction [°]
    phi_w = ncread(filename, 'puvdir'); % wave propagation direction [°]
    Urms = ncread(filename, 'u_ssm');   % rms total (u_ss+v_ss) orbital velocity [m/s] at depth z
    Urms = sqrt(2) .* Urms;             % rms peak orbital velocity amplitude at depth z

    % Calculate the longshore and cross-shore fractions of the total velocity vector
    uLong_frac = abs(uLong_z).^2 ./ (abs(uLong_z).^2 + abs(uCross_z).^2);
    uCross_frac = abs(uCross_z).^2 ./ (abs(uLong_z).^2 + abs(uCross_z).^2);

    % Estimate the depth-averaged current velocity
    k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)
    [u_c, ~] = compute_DAV(abs(u_z), z_new, k_sc, h);

    % Compute the shear stress components
    [tau_c, tau_w, tau_cw] = compute_BSS_orbital(u_c, h, Urms, T, rho_w, phi_c, phi_w, GS_stats.d50(i), g);
    % [tau_c, tau_w, tau_cw] = compute_BSS_linearTheory(u_c, H, k, T, h, rho_w, phi_c, phi_w, GS_stats.d50(i), g);

    % Append to timetable
    firstRow = find(TT_tau_c.Time == tStart);
    lastRow = find(TT_tau_c.Time == tEnd);

    TT_long_sign.(i)(firstRow:lastRow) = uLong_sign;
    TT_long_frac.(i)(firstRow:lastRow) = uLong_frac;
    TT_cross_sign.(i)(firstRow:lastRow) = uCross_sign;
    TT_cross_frac.(i)(firstRow:lastRow) = uCross_frac;

    TT_tau_c.(i)(firstRow:lastRow) = tau_c;
    TT_tau_w.(i)(firstRow:lastRow) = tau_w;
    TT_tau_cw.(i)(firstRow:lastRow) = tau_cw;

    %%%%%% Quality check %%%%%%
    % Find locations of BSS spikes
    locations_c = find(TT_tau_c.(i) > 3.0);
    locations_w = find(TT_tau_w.(i) > 3.0);

    % Change the values to NaN
    TT_tau_c.(i)(locations_c) = NaN;
    TT_tau_w.(i)(locations_w) = NaN;
    TT_tau_cw.(i)(locations_c) = NaN;
    TT_tau_cw.(i)(locations_w) = NaN;

    % Manual corrections
    TT_tau_w.L1C1VEC(3936) = NaN;
    TT_tau_cw.L1C1VEC(3936) = NaN;
    % TT_tau_w.L2C4VEC([1978, 1979]) = NaN;
    % TT_tau_cw.L2C4VEC([1978, 1979]) = NaN;
    %%%%%% Quality check %%%%%%
    
    % Fraction-specific Shields numbers
    [theta_cMg, theta_wMg, theta_cwMg] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.Mg(i), g);
    [theta_c10, theta_w10, theta_cw10] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d10(i), g);
    [theta_c50, theta_w50, theta_cw50] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d50(i), g);
    [theta_c90, theta_w90, theta_cw90] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d90(i), g);
    [theta_cFine, theta_wFine, theta_cwFine] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.Fine(i), g);
    [theta_cCoarse, theta_wCoarse, theta_cwCoarse] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.Coarse(i), g);

    TT_theta_cMg.(i)(firstRow:lastRow) = theta_cMg;
    TT_theta_wMg.(i)(firstRow:lastRow) = theta_wMg;
    TT_theta_cwMg.(i)(firstRow:lastRow) = theta_cwMg;

    TT_theta_c10.(i)(firstRow:lastRow) = theta_c10;
    TT_theta_w10.(i)(firstRow:lastRow) = theta_w10;
    TT_theta_cw10.(i)(firstRow:lastRow) = theta_cw10;

    TT_theta_c50.(i)(firstRow:lastRow) = theta_c50;
    TT_theta_w50.(i)(firstRow:lastRow) = theta_w50;
    TT_theta_cw50.(i)(firstRow:lastRow) = theta_cw50;
    
    TT_theta_c90.(i)(firstRow:lastRow) = theta_c90;
    TT_theta_w90.(i)(firstRow:lastRow) = theta_w90;
    TT_theta_cw90.(i)(firstRow:lastRow) = theta_cw90;

    TT_theta_cFine.(i)(firstRow:lastRow) = theta_cFine;
    TT_theta_wFine.(i)(firstRow:lastRow) = theta_wFine;
    TT_theta_cwFine.(i)(firstRow:lastRow) = theta_cwFine;

    TT_theta_cCoarse.(i)(firstRow:lastRow) = theta_cCoarse;
    TT_theta_wCoarse.(i)(firstRow:lastRow) = theta_wCoarse;
    TT_theta_cwCoarse.(i)(firstRow:lastRow) = theta_cwCoarse;

    % Nondimensional bedload predictors (gross)
    alpha = 11;  % calibration coefficient of Ribberink (1998)
    beta = 1.65; % calibration exponent

    phi_bMg = compute_Einstein_parameter(TT_theta_cwMg.(i), GS_stats.theta_crMg_McCarron(i), alpha, beta);
    phi_b10 = compute_Einstein_parameter(TT_theta_cw10.(i), GS_stats.theta_cr10_McCarron(i), alpha, beta);
    phi_b50 = compute_Einstein_parameter(TT_theta_cw50.(i), GS_stats.theta_cr50_McCarron(i), alpha, beta);
    phi_b90 = compute_Einstein_parameter(TT_theta_cw90.(i), GS_stats.theta_cr90_McCarron(i), alpha, beta);
    phi_bFine = compute_Einstein_parameter(TT_theta_cwFine.(i), GS_stats.theta_crFine_McCarron(i), alpha, beta);
    phi_bCoarse = compute_Einstein_parameter(TT_theta_cwCoarse.(i), GS_stats.theta_crCoarse_McCarron(i), alpha, beta);

    TT_phi_bMg.(i) = phi_bMg;
    TT_phi_b10.(i) = phi_b10;
    TT_phi_b50.(i) = phi_b50;
    TT_phi_b90.(i) = phi_b90;
    TT_phi_bFine.(i) = phi_bFine;
    TT_phi_bCoarse.(i) = phi_bCoarse;

    % Dimensional bedload transport rate (gross)
    q_bMg = compute_transport_rate(phi_bMg, rho_w, rho_s, GS_stats.Mg(i), g);
    q_b10 = compute_transport_rate(phi_b10, rho_w, rho_s, GS_stats.d10(i), g);
    q_b50 = compute_transport_rate(phi_b50, rho_w, rho_s, GS_stats.d50(i), g);
    q_b90 = compute_transport_rate(phi_b90, rho_w, rho_s, GS_stats.d90(i), g);
    q_bFine = compute_transport_rate(phi_bFine, rho_w, rho_s, GS_stats.Fine(i), g);
    q_bCoarse = compute_transport_rate(phi_bCoarse, rho_w, rho_s, GS_stats.Coarse(i), g);

    TT_q_bMg.(i) = q_bMg;
    TT_q_b10.(i) = q_b10;
    TT_q_b50.(i) = q_b50;
    TT_q_b90.(i) = q_b90;
    TT_q_bFine.(i) = q_bFine;
    TT_q_bCoarse.(i) = q_bCoarse;

    % Dimensional bedload transport rate (net: longshore direction)
    qL_bMg_net = TT_long_sign.(i) .* q_bMg .* TT_long_frac.(i);
    qL_b10_net = TT_long_sign.(i) .* q_b10 .* TT_long_frac.(i);
    qL_b50_net = TT_long_sign.(i) .* q_b50 .* TT_long_frac.(i);
    qL_b90_net = TT_long_sign.(i) .* q_b90 .* TT_long_frac.(i);
    qL_bFine_net = TT_long_sign.(i) .* q_bFine .* TT_long_frac.(i);
    qL_bCoarse_net = TT_long_sign.(i) .* q_bCoarse .* TT_long_frac.(i);

    TT_qL_bMg_net.(i) = qL_bMg_net;
    TT_qL_b10_net.(i) = qL_b10_net;
    TT_qL_b50_net.(i) = qL_b50_net;
    TT_qL_b90_net.(i) = qL_b90_net;
    TT_qL_bFine_net.(i) = qL_bFine_net;
    TT_qL_bCoarse_net.(i) = qL_bCoarse_net;

    % Dimensional bedload transport rate (net: cross-shore direction)
    qC_bMg_net = TT_cross_sign.(i) .* q_bMg .* TT_cross_frac.(i);
    qC_b10_net = TT_cross_sign.(i) .* q_b10 .* TT_cross_frac.(i);
    qC_b50_net = TT_cross_sign.(i) .* q_b50 .* TT_cross_frac.(i);
    qC_b90_net = TT_cross_sign.(i) .* q_b90 .* TT_cross_frac.(i);
    qC_bFine_net = TT_cross_sign.(i) .* q_bFine .* TT_cross_frac.(i);
    qC_bCoarse_net = TT_cross_sign.(i) .* q_bCoarse .* TT_cross_frac.(i);

    TT_qC_bMg_net.(i) = qC_bMg_net;
    TT_qC_b10_net.(i) = qC_b10_net;
    TT_qC_b50_net.(i) = qC_b50_net;
    TT_qC_b90_net.(i) = qC_b90_net;
    TT_qC_bFine_net.(i) = qC_bFine_net;
    TT_qC_bCoarse_net.(i) = qC_bCoarse_net;

end

% Manual corrections
TT_cross_frac.L4C1VEC(1:97) = NaN;
TT_cross_sign.L4C1VEC(1:97) = NaN;

% clearvars dataPath dt emptyData filename firstRow h H i info k lastRow n phi_c phi_w t T tau_c tau_w tau_cw theta_c10 theta_c50 theta_c90 theta_cFine theta_cCoarse theta_w10 theta_w50 theta_w90 theta_wFine theta_wCoarse theta_cw10 theta_cw50 theta_cw90 theta_cwFine theta_cwCoarse timeAxis rho_w u_z u_c h Urms emptyData firstRow lastRow z z_new z_measured nonNaNValues nonNaNIndices interpolatedValues phi_b10 phi_b50 phi_b90 phi_bFine phi_bCoarse q_b10 q_b50 q_b90 q_bFine q_bCoarse uLong_z uLong_sign uLong_frac uCross_z uCross_sign uCross_frac q_b10_net q_b50_net q_b90_net q_bFine_net q_bCoarse_net suffix suffixes timetableName tStart tEnd time t0 habCVnew theta_cMg theta_wMg theta_cwMg phi_bMg q_bMg q_bMg_net locations_c locations_w qL_bMg_net qL_b10_net qL_b50_net qL_b90_net qL_bFine_net qL_bCoarse_net qC_bMg_net qC_b10_net qC_b50_net qC_b90_net qC_bFine_net qC_bCoarse_net


%% Longshore vs Cross-shore flow
% Low percentage -> flood + offshore & ebb + onshore
% High percentage -> flood + onshore & ebb + offshore

for i = 1:6
    series1 = TT_long_sign{:,i};
    series2 = TT_cross_sign{:,i};
    
    % Find matching '1's
    both_ones = (series1 == 1) & (series2 == 1);
    
    % Find matching '-1's
    both_minus_ones = (series1 == -1) & (series2 == -1);
    
    % Count occurrences where both are '1'
    count_ones = sum(both_ones);
    
    % Count occurrences where both are '-1'
    count_minus_ones = sum(both_minus_ones);
    
    total_counts = sum(~isnan(series1) & ~isnan(series2));
    
    % Total matching occurrences
    total_matches = count_ones + count_minus_ones;
    match_perc = total_matches ./ total_counts .* 100;
    
    % Display the result
    fprintf('Total matching occurrences: %d (%.f%%)\n', total_matches, match_perc);
end


%% Only consider time steps without NaNs
rowsWithoutNaNs = ~any(ismissing(TT_tau_cw), 2);

TT_theta_cwMg_NaN = TT_theta_cwMg(rowsWithoutNaNs, :);
TT_theta_cw10_NaN = TT_theta_cw10(rowsWithoutNaNs, :);
TT_theta_cw50_NaN = TT_theta_cw50(rowsWithoutNaNs, :);
TT_theta_cw90_NaN = TT_theta_cw90(rowsWithoutNaNs, :);
TT_theta_cwFine_NaN = TT_theta_cwFine(rowsWithoutNaNs, :);
TT_theta_cwCoarse_NaN = TT_theta_cwCoarse(rowsWithoutNaNs, :);


%% Only consider stormy period

% Find indices within the date range
idx = (TT_theta_cMg.Time >= stormy_period(1)) & (TT_theta_cMg.Time <= stormy_period(2));
% idx = (TT_theta_cMg.Time >= storm_period(1)) & (TT_theta_cMg.Time <= storm_period(2));

% Get the first and last index
firstIdx = find(idx, 1, 'first');
lastIdx = find(idx, 1, 'last');

nonStormIdx = [1:firstIdx-1, lastIdx+1:height(TT_theta_cMg)];
stormIdx = firstIdx:lastIdx;

% Define suffixes to iterate over
suffixes = {'Coarse', 'Fine', 'Mg', '90', '50', '10'};

for i = 1:length(suffixes)
    suffix = suffixes{i};
    
    % Dynamically generate variable names for storm data
    eval(['TT_theta_cw' suffix '_NaN_storm = TT_theta_cw' suffix '(stormIdx, :);']);
    eval(['TT_theta_cw' suffix '_NaN_storm = TT_theta_cw' suffix '_NaN_storm(rowsWithoutNaNs(stormIdx), :);']);
    
    % Dynamically generate variable names for non-storm data
    eval(['TT_theta_cw' suffix '_NaN_noStorm = TT_theta_cw' suffix '(nonStormIdx, :);']);
    eval(['TT_theta_cw' suffix '_NaN_noStorm = TT_theta_cw' suffix '_NaN_noStorm(rowsWithoutNaNs(nonStormIdx), :);']);
end


%% Time periods and frequency
% Get the frequency of the time axis
dt = TT_theta_cw50.Properties.TimeStep;
totalDuration = dateEnd-dateStart;

% Find rows without any NaNs
numRowsWithoutNaNs = sum(rowsWithoutNaNs);
totalDuration_noNaN = numRowsWithoutNaNs * dt;
totalDurationFraction = totalDuration_noNaN / totalDuration;

numRowsWithoutNaNs_noStorm = sum(rowsWithoutNaNs(nonStormIdx));
totalDuration_noNaN_noStorm = numRowsWithoutNaNs_noStorm * dt;
totalDurationFraction_noStorm = totalDuration_noNaN_noStorm / totalDuration_noNaN;

numRowsWithoutNaNs_storm = sum(rowsWithoutNaNs(stormIdx));
totalDuration_noNaN_storm = numRowsWithoutNaNs_storm * dt;
totalDurationFraction_storm = totalDuration_noNaN_storm / totalDuration_noNaN;


%% Computations: mobilisation duration (Shields)

% Create an empty table
percentExceed_Soulsby = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);
percentExceed_Egiazaroff = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);
percentExceed_McCarron = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);
percentExceed_McCarron_noStorm = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);
percentExceed_McCarron_storm = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);

% Calculate the time percentages of exceedance
disp('McCarron:')
for i = 1:width(TT_theta_cwMg_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_Soulsby = sum(TT_theta_cwMg_NaN.(i) > GS_stats.theta_crMg_Soulsby(i));
    exceededTimes_Egiazaroff = sum(TT_theta_cwMg_NaN.(i) > GS_stats.theta_crMg_Soulsby(i));
    exceededTimes_McCarron = sum(TT_theta_cwMg_NaN.(i) > GS_stats.theta_crMg_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
    percent_Mg_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
    exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
    percent_Mg_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
    exceededDuration_McCarron = exceededTimes_McCarron*dt;
    percent_Mg_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;

    % Append the percentage value to the table
    percentExceed_Soulsby.Mg(i) = percent_Mg_Soulsby;
    percentExceed_Egiazaroff.Mg(i) = percent_Mg_Egiazaroff;
    percentExceed_McCarron.Mg(i) = percent_Mg_McCarron;

    % Display the total duration
    disp(['Threshold Mg at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
        ' h (', num2str(percent_Mg_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);

end

for i = 1:width(TT_theta_cw10_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_Soulsby = sum(TT_theta_cw10_NaN.(i) > GS_stats.theta_cr10_Soulsby(i));
    exceededTimes_Egiazaroff = sum(TT_theta_cw10_NaN.(i) > GS_stats.theta_cr10_Soulsby(i));
    exceededTimes_McCarron = sum(TT_theta_cw10_NaN.(i) > GS_stats.theta_cr10_McCarron(i));
        
    % Find the duration when the threshold is exceeded
    exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
    percent_d10_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
    exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
    percent_d10_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
    exceededDuration_McCarron = exceededTimes_McCarron*dt;
    percent_d10_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;

    % Append the percentage value to the table
    percentExceed_Soulsby.d10(i) = percent_d10_Soulsby;
    percentExceed_Egiazaroff.d10(i) = percent_d10_Egiazaroff;
    percentExceed_McCarron.d10(i) = percent_d10_McCarron;
    
    % Display the total duration
    disp(['Threshold d10 at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
        ' h (', num2str(percent_d10_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);

end

for i = 1:width(TT_theta_cw50_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_Soulsby = sum(TT_theta_cw50_NaN.(i) > GS_stats.theta_cr50_Soulsby(i));
    exceededTimes_Egiazaroff = sum(TT_theta_cw50_NaN.(i) > GS_stats.theta_cr50_Egiazaroff(i));
    exceededTimes_McCarron = sum(TT_theta_cw50_NaN.(i) > GS_stats.theta_cr50_McCarron(i));
        
    % Find the duration when the threshold is exceeded
    exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
    percent_d50_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
    exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
    percent_d50_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
    exceededDuration_McCarron = exceededTimes_McCarron*dt;
    percent_d50_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;

    % Append the percentage value to the table
    percentExceed_Soulsby.d50(i) = percent_d50_Soulsby;
    percentExceed_Egiazaroff.d50(i) = percent_d50_Egiazaroff;
    percentExceed_McCarron.d50(i) = percent_d50_McCarron;
    
    % Display the total duration
    disp(['Threshold d50 at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
        ' h (', num2str(percent_d50_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);

end

for i = 1:width(TT_theta_cw90_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_Soulsby = sum(TT_theta_cw90_NaN.(i) > GS_stats.theta_cr90_Soulsby(i));
    exceededTimes_Egiazaroff = sum(TT_theta_cw90_NaN.(i) > GS_stats.theta_cr90_Egiazaroff(i));
    exceededTimes_McCarron = sum(TT_theta_cw90_NaN.(i) > GS_stats.theta_cr90_McCarron(i));
        
    % Find the duration when the threshold is exceeded
    exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
    percent_d90_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
    exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
    percent_d90_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
    exceededDuration_McCarron = exceededTimes_McCarron*dt;
    percent_d90_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;

    % Append the percentage value to the table
    percentExceed_Soulsby.d90(i) = percent_d90_Soulsby;
    percentExceed_Egiazaroff.d90(i) = percent_d90_Egiazaroff;
    percentExceed_McCarron.d90(i) = percent_d90_McCarron;
    
    % Display the total duration
    disp(['Threshold d90 at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
        ' h (', num2str(percent_d90_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);

end

for i = 1:width(TT_theta_cwFine_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_Soulsby = sum(TT_theta_cwFine_NaN.(i) > GS_stats.theta_crFine_Soulsby(i));
    exceededTimes_Egiazaroff = sum(TT_theta_cwFine_NaN.(i) > GS_stats.theta_crFine_Egiazaroff(i));
    exceededTimes_McCarron = sum(TT_theta_cwFine_NaN.(i) > GS_stats.theta_crFine_McCarron(i));
        
    % Find the duration when the threshold is exceeded
    exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
    percent_Fine_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
    exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
    percent_Fine_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
    exceededDuration_McCarron = exceededTimes_McCarron*dt;
    percent_Fine_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;

    % Append the percentage value to the table
    percentExceed_Soulsby.Fine(i) = percent_Fine_Soulsby;
    percentExceed_Egiazaroff.Fine(i) = percent_Fine_Egiazaroff;
    percentExceed_McCarron.Fine(i) = percent_Fine_McCarron;
    
    % Display the total duration
    disp(['Threshold Fine at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
        ' h (', num2str(percent_Fine_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);

end

for i = 1:width(TT_theta_cwCoarse_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_Soulsby = sum(TT_theta_cwCoarse_NaN.(i) > GS_stats.theta_crCoarse_Soulsby(i));
    exceededTimes_Egiazaroff = sum(TT_theta_cwCoarse_NaN.(i) > GS_stats.theta_crCoarse_Egiazaroff(i));
    exceededTimes_McCarron = sum(TT_theta_cwCoarse_NaN.(i) > GS_stats.theta_crCoarse_McCarron(i));
        
    % Find the duration when the threshold is exceeded
    exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
    percent_Coarse_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
    exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
    percent_Coarse_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
    exceededDuration_McCarron = exceededTimes_McCarron*dt;
    percent_Coarse_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;

    % Append the percentage value to the table
    percentExceed_Soulsby.Coarse(i) = percent_Coarse_Soulsby;
    percentExceed_Egiazaroff.Coarse(i) = percent_Coarse_Egiazaroff;
    percentExceed_McCarron.Coarse(i) = percent_Coarse_McCarron;
    
    % Display the total duration
    disp(['Threshold Coarse at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
        ' h (', num2str(percent_Coarse_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);

end

for i = 1:width(TT_theta_cwMg_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_McCarron_storm = sum(TT_theta_cwMg_NaN_storm.(i) > GS_stats.theta_crMg_McCarron(i));
    exceededTimes_McCarron_noStorm = sum(TT_theta_cwMg_NaN_noStorm.(i) > GS_stats.theta_crMg_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_McCarron_storm = exceededTimes_McCarron_storm*dt;
    percent_Mg_McCarron_storm = exceededDuration_McCarron_storm/totalDuration_noNaN_storm*100;
    exceededDuration_McCarron_noStorm = exceededTimes_McCarron_noStorm*dt;
    percent_Mg_McCarron_noStorm = exceededDuration_McCarron_noStorm/totalDuration_noNaN_noStorm*100;

    % Append the percentage value to the table
    percentExceed_McCarron_storm.Mg(i) = percent_Mg_McCarron_storm;
    percentExceed_McCarron_noStorm.Mg(i) = percent_Mg_McCarron_noStorm;

end

for i = 1:width(TT_theta_cw10_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_McCarron_storm = sum(TT_theta_cw10_NaN_storm.(i) > GS_stats.theta_cr10_McCarron(i));
    exceededTimes_McCarron_noStorm = sum(TT_theta_cw10_NaN_noStorm.(i) > GS_stats.theta_cr10_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_McCarron_storm = exceededTimes_McCarron_storm*dt;
    percent_10_McCarron_storm = exceededDuration_McCarron_storm/totalDuration_noNaN_storm*100;
    exceededDuration_McCarron_noStorm = exceededTimes_McCarron_noStorm*dt;
    percent_10_McCarron_noStorm = exceededDuration_McCarron_noStorm/totalDuration_noNaN_noStorm*100;

    % Append the percentage value to the table
    percentExceed_McCarron_storm.d10(i) = percent_10_McCarron_storm;
    percentExceed_McCarron_noStorm.d10(i) = percent_10_McCarron_noStorm;

end

for i = 1:width(TT_theta_cw50_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_McCarron_storm = sum(TT_theta_cw50_NaN_storm.(i) > GS_stats.theta_cr50_McCarron(i));
    exceededTimes_McCarron_noStorm = sum(TT_theta_cw50_NaN_noStorm.(i) > GS_stats.theta_cr50_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_McCarron_storm = exceededTimes_McCarron_storm*dt;
    percent_50_McCarron_storm = exceededDuration_McCarron_storm/totalDuration_noNaN_storm*100;
    exceededDuration_McCarron_noStorm = exceededTimes_McCarron_noStorm*dt;
    percent_50_McCarron_noStorm = exceededDuration_McCarron_noStorm/totalDuration_noNaN_noStorm*100;

    % Append the percentage value to the table
    percentExceed_McCarron_storm.d50(i) = percent_50_McCarron_storm;
    percentExceed_McCarron_noStorm.d50(i) = percent_50_McCarron_noStorm;

end

for i = 1:width(TT_theta_cw90_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_McCarron_storm = sum(TT_theta_cw90_NaN_storm.(i) > GS_stats.theta_cr90_McCarron(i));
    exceededTimes_McCarron_noStorm = sum(TT_theta_cw90_NaN_noStorm.(i) > GS_stats.theta_cr90_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_McCarron_storm = exceededTimes_McCarron_storm*dt;
    percent_90_McCarron_storm = exceededDuration_McCarron_storm/totalDuration_noNaN_storm*100;
    exceededDuration_McCarron_noStorm = exceededTimes_McCarron_noStorm*dt;
    percent_90_McCarron_noStorm = exceededDuration_McCarron_noStorm/totalDuration_noNaN_noStorm*100;

    % Append the percentage value to the table
    percentExceed_McCarron_storm.d90(i) = percent_90_McCarron_storm;
    percentExceed_McCarron_noStorm.d90(i) = percent_90_McCarron_noStorm;

end

for i = 1:width(TT_theta_cwFine_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_McCarron_storm = sum(TT_theta_cwFine_NaN_storm.(i) > GS_stats.theta_crFine_McCarron(i));
    exceededTimes_McCarron_noStorm = sum(TT_theta_cwFine_NaN_noStorm.(i) > GS_stats.theta_crFine_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_McCarron_storm = exceededTimes_McCarron_storm*dt;
    percent_Fine_McCarron_storm = exceededDuration_McCarron_storm/totalDuration_noNaN_storm*100;
    exceededDuration_McCarron_noStorm = exceededTimes_McCarron_noStorm*dt;
    percent_Fine_McCarron_noStorm = exceededDuration_McCarron_noStorm/totalDuration_noNaN_noStorm*100;

    % Append the percentage value to the table
    percentExceed_McCarron_storm.Fine(i) = percent_Fine_McCarron_storm;
    percentExceed_McCarron_noStorm.Fine(i) = percent_Fine_McCarron_noStorm;

end

for i = 1:width(TT_theta_cwCoarse_NaN)

    % Find the moments when the threshold is exceeded
    exceededTimes_McCarron_storm = sum(TT_theta_cwCoarse_NaN_storm.(i) > GS_stats.theta_crCoarse_McCarron(i));
    exceededTimes_McCarron_noStorm = sum(TT_theta_cwCoarse_NaN_noStorm.(i) > GS_stats.theta_crCoarse_McCarron(i));

    % Find the duration when the threshold is exceeded
    exceededDuration_McCarron_storm = exceededTimes_McCarron_storm*dt;
    percent_Coarse_McCarron_storm = exceededDuration_McCarron_storm/totalDuration_noNaN_storm*100;
    exceededDuration_McCarron_noStorm = exceededTimes_McCarron_noStorm*dt;
    percent_Coarse_McCarron_noStorm = exceededDuration_McCarron_noStorm/totalDuration_noNaN_noStorm*100;

    % Append the percentage value to the table
    percentExceed_McCarron_storm.Coarse(i) = percent_Coarse_McCarron_storm;
    percentExceed_McCarron_noStorm.Coarse(i) = percent_Coarse_McCarron_noStorm;

end

% Calculate mean row and column
meanRow_Soulsby = mean(percentExceed_Soulsby.Variables, 1);
meanCol_Soulsby = mean(percentExceed_Soulsby.Variables, 2);

meanRow_Egiazaroff = mean(percentExceed_Egiazaroff.Variables, 1);
meanCol_Egiazaroff = mean(percentExceed_Egiazaroff.Variables, 2);

meanRow_McCarron = mean(percentExceed_McCarron.Variables, 1);
meanCol_McCarron = mean(percentExceed_McCarron.Variables, 2);

meanRow_McCarron_storm = mean(percentExceed_McCarron_storm.Variables, 1);
meanCol_McCarron_storm = mean(percentExceed_McCarron_storm.Variables, 2);

meanRow_McCarron_noStorm = mean(percentExceed_McCarron_noStorm.Variables, 1);
meanCol_McCarron_noStorm = mean(percentExceed_McCarron_noStorm.Variables, 2);

% Add mean row and column to the table
percentExceed_Soulsby{'Mean', :} = meanRow_Soulsby;
percentExceed_Soulsby.Mean = [meanCol_Soulsby; mean(meanRow_Soulsby)];

percentExceed_Egiazaroff{'Mean', :} = meanRow_Egiazaroff;
percentExceed_Egiazaroff.Mean = [meanCol_Egiazaroff; mean(meanRow_Egiazaroff)];

percentExceed_McCarron{'Mean', :} = meanRow_McCarron;
percentExceed_McCarron.Mean = [meanCol_McCarron; mean(meanRow_McCarron)];

percentExceed_McCarron_storm{'Mean', :} = meanRow_McCarron_storm;
percentExceed_McCarron_storm.Mean = [meanCol_McCarron_storm; mean(meanRow_McCarron_storm)];

percentExceed_McCarron_noStorm{'Mean', :} = meanRow_McCarron_noStorm;
percentExceed_McCarron_noStorm.Mean = [meanCol_McCarron_noStorm; mean(meanRow_McCarron_noStorm)];


clearvars dt exceededDuration_Soulsby exceededDuration_Egiazaroff exceededDuration_McCarron exceededTimes_Soulsby exceededTimes_Egiazaroff exceededTimes_McCarron i percent_d10_Soulsby percent_d10_Egiazaroff percent_d10_McCarron percent_d50_Soulsby percent_d50_Egiazaroff percent_d50_McCarron percent_d90_Soulsby percent_d90_Egiazaroff percent_d90_McCarron meanCol_Soulsby meanCol_Egiazaroff meanCol_McCarron meanRow_Soulsby meanRow_Egiazaroff meanRow_McCarron totalDuration totalDurationFraction numRowsWithoutNaNs rowsWithoutNaNs


%% Visualisation: Bed Shear Stress (timeseries)
ax = gobjects(1, length(instruLocs));

f3a = figure('Position',[740, 957, 1719, 1336]);
tl = tiledlayout(length(instruLocs), 1, 'TileSpacing','tight');

for i = 1:length(instruLocs)
    above_d00 = TT_tau_cw.(i)<=GS_stats.tau_cr50_McCarron(i);
    above_d10 = TT_tau_cw.(i)>GS_stats.tau_cr10_McCarron(i) & TT_tau_cw.(i)<=GS_stats.tau_cr90_McCarron(i);
    above_d50 = TT_tau_cw.(i)>GS_stats.tau_cr50_McCarron(i) & TT_tau_cw.(i)<=GS_stats.tau_cr10_McCarron(i);
    above_d90 = TT_tau_cw.(i)>GS_stats.tau_cr90_McCarron(i);

    tau_d10 = TT_tau_cw.(i);
    tau_d50 = TT_tau_cw.(i);
    tau_d90 = TT_tau_cw.(i);
    tau_d00 = TT_tau_cw.(i);

    tau_d10(~above_d10) = NaN;
    tau_d50(~above_d50) = NaN;
    tau_d90(~above_d90) = NaN;
    tau_d00(~above_d00) = NaN;

    ax(i) = nexttile;
    % plot(TT_tau_cw.Time, TT_tau_cw.(i), 'LineWidth',2, 'Color',newcolors(i,:))
    plot(TT_tau_cw.Time, tau_d00, '-k', 'LineWidth',2); hold on
    plot(TT_tau_cw.Time, tau_d50, '-b', 'LineWidth',2)
    plot(TT_tau_cw.Time, tau_d10, '-r', 'LineWidth',2)
    plot(TT_tau_cw.Time, tau_d90, '-m', 'LineWidth',2)

    % scatter(TT_tau_cw.Time(above_d00), TT_tau_cw.(i)(above_d00), 30, 'filled', 'MarkerFaceColor','k'); hold on
    % scatter(TT_tau_cw.Time(above_d90), TT_tau_cw.(i)(above_d90), 30, 'filled', 'MarkerFaceColor','g'); hold on
    % scatter(TT_tau_cw.Time(above_d50), TT_tau_cw.(i)(above_d50), 30, 'filled', 'MarkerFaceColor','b'); hold on
    % scatter(TT_tau_cw.Time(above_d10), TT_tau_cw.(i)(above_d10), 30, 'filled', 'MarkerFaceColor','r'); hold on

    yline(GS_stats.tau_cr10_McCarron(i), '-.', 'LineWidth',2)
    yline(GS_stats.tau_cr50_McCarron(i), ':', 'LineWidth',2)
    yline(GS_stats.tau_cr90_McCarron(i), '--', 'LineWidth',2)

    % text(tEnd+hours(5), GS_stats.tau_cr10_McCarron(i), '\tau_{cr,10}', 'FontSize',fontsize*.8)
    % text(tEnd+hours(5), GS_stats.tau_cr50_McCarron(i), '\tau_{cr,50}', 'FontSize',fontsize*.8)
    % text(tEnd+hours(5), GS_stats.tau_cr90_McCarron(i), '\tau_{cr,90}', 'FontSize',fontsize*.8)

    xticks(dateStart+1:days(2):dateEnd)
    xlim([dateStart, dateEnd])
    ylim([0, 3])

    text(dateEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

    grid on
    box off
end
xticklabels(ax(1:end-1), [])
ylabel(tl, '\tau_{cw} (Pa)', 'FontSize',fontsize)
legend(ax(1), '', '', '', '', '\tau_{cr,10}', '\tau_{cr,50}', '\tau_{cr,90}', 'Location','northoutside', 'NumColumns',3)
zoom xon
linkaxes

clearvars above_d00 above_d10 above_d50 above_d90 tau_d00 tau_d10 tau_d50 tau_d90 tl i ax


%% Visualisation: Bed Shear Stress (pdf)
% f3b = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_tau_cw.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Create the KDE plot
%     [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
%         'Support','nonnegative');
%     plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
% 
% end
% hold off
% xlim([0, 4])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Bed Shear Stress')
% xlabel('\tau_{cw} (Pa)')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i


%% Visualisation: Nondimensional Bedload Transport Rate (timeseries)
ax = gobjects(1, length(instruLocs));

f4a = figure('Position',[740, 957, 1719, 1336]);
tl = tiledlayout(length(instruLocs), 1, 'TileSpacing','tight');

for i = 1:length(instruLocs)

    ax(i) = nexttile;
    plot(TT_phi_b10.Time, TT_phi_b10.(i), '-k', 'LineWidth',2); hold on
    plot(TT_phi_b50.Time, TT_phi_b50.(i), '-r', 'LineWidth',2)
    plot(TT_phi_b90.Time, TT_phi_b90.(i), '-b', 'LineWidth',2)

    xticks(dateStart+1:days(2):dateEnd)
    xlim([dateStart, dateEnd])
    ylim([0, 3])

    text(dateEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

    grid on
    box off
    
end
xticklabels(ax(1:end-1), [])
ylabel(tl, '\phi_{b}', 'FontSize',fontsize)
legend(ax(1), 'D_{10}', 'D_{50}', 'D_{90}', 'Location','northoutside', 'NumColumns',3)
zoom xon
linkaxes

clearvars tl ax i


%% Visualisation: Nondimensional Bedload Transport Rate (pdf)
% f4b = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_phi_b50.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Create the KDE plot
%     [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
%         'Support','nonnegative');
%     plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
% 
% end
% hold off
% xlim([0, 1])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Nondimensional Bedload Transport Rate')
% xlabel('\phi_{b}')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i


%% Visualisation: Dimensional Bedload Transport Rate (timeseries)
ax = gobjects(1, length(instruLocs));

f5a = figure('Position',[740, 957, 1719, 1336]);
tl = tiledlayout(length(instruLocs), 1, 'TileSpacing','tight');

for i = 1:length(instruLocs)

    ax(i) = nexttile;
    plot(TT_q_b10.Time, TT_q_b10.(i), '-k', 'LineWidth',2); hold on
    plot(TT_q_b50.Time, TT_q_b50.(i), '-r', 'LineWidth',2)
    plot(TT_q_b90.Time, TT_q_b90.(i), '-b', 'LineWidth',2)

    xticks(dateStart+1:days(2):dateEnd)
    xlim([dateStart, dateEnd])
    ylim([0, 8e-5])

    text(dateEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

    grid on
    box off
    
end
xticklabels(ax(1:end-1), [])
ylabel(tl, 'q_{b} (m^2 s^{-1})', 'FontSize',fontsize)
legend(ax(1), 'D_{10}', 'D_{50}', 'D_{90}', 'Location','northoutside', 'NumColumns',3)
zoom xon
linkaxes

clearvars tl ax i


%% Visualisation: Dimensional Bedload Transport Rate (pdf)
% f5b = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_q_b50.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Create the KDE plot
%     [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
%         'Support','nonnegative');
%     plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
% end
% hold off
% xlim([0, 1e-4])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Bedload Transport Rate')
% xlabel('q_{b} (m^2 s^{-1})')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i


%% Find rows with any NaNs and replace with zeros

% Find rows with one or more NaNs
nanRows = any(ismissing(TT_tau_cw), 2);

% Get indices of rows with NaNs
nanRowIndices = find(nanRows);

% Initiate new tables
TT_tau_c_NaN = TT_tau_c;
TT_tau_w_NaN = TT_tau_w;
TT_tau_cw_NaN = TT_tau_cw;

TT_theta_cMg_NaN = TT_theta_cMg;
TT_theta_c10_NaN = TT_theta_c10;
TT_theta_c50_NaN = TT_theta_c50;
TT_theta_c90_NaN = TT_theta_c90;
TT_theta_cFine_NaN = TT_theta_cFine;
TT_theta_cCoarse_NaN = TT_theta_cCoarse;

TT_theta_wMg_NaN = TT_theta_wMg;
TT_theta_w10_NaN = TT_theta_w10;
TT_theta_w50_NaN = TT_theta_w50;
TT_theta_w90_NaN = TT_theta_w90;
TT_theta_wFine_NaN = TT_theta_wFine;
TT_theta_wCoarse_NaN = TT_theta_wCoarse;

TT_theta_cwMg_NaN = TT_theta_cwMg;
TT_theta_cw10_NaN = TT_theta_cw10;
TT_theta_cw50_NaN = TT_theta_cw50;
TT_theta_cw90_NaN = TT_theta_cw90;
TT_theta_cwFine_NaN = TT_theta_cwFine;
TT_theta_cwCoarse_NaN = TT_theta_cwCoarse;

TT_phi_bMg_NaN = TT_phi_bMg;
TT_phi_b10_NaN = TT_phi_b10;
TT_phi_b50_NaN = TT_phi_b50;
TT_phi_b90_NaN = TT_phi_b90;
TT_phi_bFine_NaN = TT_phi_bFine;
TT_phi_bCoarse_NaN = TT_phi_bCoarse;

TT_q_bMg_NaN = TT_q_bMg;
TT_q_b10_NaN = TT_q_b10;
TT_q_b50_NaN = TT_q_b50;
TT_q_b90_NaN = TT_q_b90;
TT_q_bFine_NaN = TT_q_bFine;
TT_q_bCoarse_NaN = TT_q_bCoarse;

TT_q_bMg_zeros = TT_q_bMg;
TT_q_b10_zeros = TT_q_b10;
TT_q_b50_zeros = TT_q_b50;
TT_q_b90_zeros = TT_q_b90;
TT_q_bFine_zeros = TT_q_bFine;
TT_q_bCoarse_zeros = TT_q_bCoarse;

TT_qL_bMg_net_zeros = TT_qL_bMg_net;
TT_qL_b10_net_zeros = TT_qL_b10_net;
TT_qL_b50_net_zeros = TT_qL_b50_net;
TT_qL_b90_net_zeros = TT_qL_b90_net;
TT_qL_bFine_net_zeros = TT_qL_bFine_net;
TT_qL_bCoarse_net_zeros = TT_qL_bCoarse_net;

TT_qC_bMg_net_zeros = TT_qC_bMg_net;
TT_qC_b10_net_zeros = TT_qC_b10_net;
TT_qC_b50_net_zeros = TT_qC_b50_net;
TT_qC_b90_net_zeros = TT_qC_b90_net;
TT_qC_bFine_net_zeros = TT_qC_bFine_net;
TT_qC_bCoarse_net_zeros = TT_qC_bCoarse_net;

% Assign NaNs to the entire row
TT_tau_c_NaN{nanRowIndices, :} = NaN;
TT_tau_w_NaN{nanRowIndices, :} = NaN;
TT_tau_cw_NaN{nanRowIndices, :} = NaN;

TT_theta_cMg_NaN{nanRowIndices, :} = NaN;
TT_theta_c10_NaN{nanRowIndices, :} = NaN;
TT_theta_c50_NaN{nanRowIndices, :} = NaN;
TT_theta_c90_NaN{nanRowIndices, :} = NaN;
TT_theta_cFine_NaN{nanRowIndices, :} = NaN;
TT_theta_c90_NaN{nanRowIndices, :} = NaN;

TT_theta_wMg_NaN{nanRowIndices, :} = NaN;
TT_theta_w10_NaN{nanRowIndices, :} = NaN;
TT_theta_w50_NaN{nanRowIndices, :} = NaN;
TT_theta_w90_NaN{nanRowIndices, :} = NaN;
TT_theta_wFine_NaN{nanRowIndices, :} = NaN;
TT_theta_wCoarse_NaN{nanRowIndices, :} = NaN;

TT_theta_cwMg_NaN{nanRowIndices, :} = NaN;
TT_theta_cw10_NaN{nanRowIndices, :} = NaN;
TT_theta_cw50_NaN{nanRowIndices, :} = NaN;
TT_theta_cw90_NaN{nanRowIndices, :} = NaN;
TT_theta_cwFine_NaN{nanRowIndices, :} = NaN;
TT_theta_cwCoarse_NaN{nanRowIndices, :} = NaN;

TT_phi_bMg_NaN{nanRowIndices, :} = NaN;
TT_phi_b10_NaN{nanRowIndices, :} = NaN;
TT_phi_b50_NaN{nanRowIndices, :} = NaN;
TT_phi_b90_NaN{nanRowIndices, :} = NaN;
TT_phi_bFine_NaN{nanRowIndices, :} = NaN;
TT_phi_bCoarse_NaN{nanRowIndices, :} = NaN;

TT_q_bMg_NaN{nanRowIndices, :} = NaN;
TT_q_b10_NaN{nanRowIndices, :} = NaN;
TT_q_b50_NaN{nanRowIndices, :} = NaN;
TT_q_b90_NaN{nanRowIndices, :} = NaN;
TT_q_bFine_NaN{nanRowIndices, :} = NaN;
TT_q_bCoarse_NaN{nanRowIndices, :} = NaN;

% Assign zeros to the entire row
TT_q_bMg_zeros{nanRowIndices, :} = 0;
TT_q_b10_zeros{nanRowIndices, :} = 0;
TT_q_b50_zeros{nanRowIndices, :} = 0;
TT_q_b90_zeros{nanRowIndices, :} = 0;
TT_q_bFine_zeros{nanRowIndices, :} = 0;
TT_q_bCoarse_zeros{nanRowIndices, :} = 0;

TT_qL_bMg_net_zeros{nanRowIndices, :} = 0;
TT_qL_b10_net_zeros{nanRowIndices, :} = 0;
TT_qL_b50_net_zeros{nanRowIndices, :} = 0;
TT_qL_b90_net_zeros{nanRowIndices, :} = 0;
TT_qL_bFine_net_zeros{nanRowIndices, :} = 0;
TT_qL_bCoarse_net_zeros{nanRowIndices, :} = 0;

TT_qC_bMg_net_zeros{nanRowIndices, :} = 0;
TT_qC_b10_net_zeros{nanRowIndices, :} = 0;
TT_qC_b50_net_zeros{nanRowIndices, :} = 0;
TT_qC_b90_net_zeros{nanRowIndices, :} = 0;
TT_qC_bFine_net_zeros{nanRowIndices, :} = 0;
TT_qC_bCoarse_net_zeros{nanRowIndices, :} = 0;

clearvars nanRows nanRowIndices


%% Calculation: Total Bedload Transport Volume (Gross)
Mg = NaN(length(instruLocs), 1);
d10 = NaN(length(instruLocs), 1);
d50 = NaN(length(instruLocs), 1);
d90 = NaN(length(instruLocs), 1);
Fine = NaN(length(instruLocs), 1);
Coarse = NaN(length(instruLocs), 1);

% Convert datetime to duration and then to seconds
time_seconds = seconds(TT_q_b10_zeros.Time - TT_q_b10_zeros.Time(1));

for i = 1:length(instruLocs)
    
    % Integrate using the trapezoidal rule
    Mg(i) = trapz(time_seconds, TT_q_bMg_zeros.(i));
    d10(i) = trapz(time_seconds, TT_q_b10_zeros.(i));
    d50(i) = trapz(time_seconds, TT_q_b50_zeros.(i));
    d90(i) = trapz(time_seconds, TT_q_b90_zeros.(i));
    Fine(i) = trapz(time_seconds, TT_q_bFine_zeros.(i));
    Coarse(i) = trapz(time_seconds, TT_q_bCoarse_zeros.(i));

    % Display the result
    % fprintf('Total amount of bedload transport: %.2f m^3/m width\n', totalBedload);
end

totalBedload_gross = table(Mg, d10, d50, d90, Fine, Coarse, 'RowNames', flipud(instruLocs));

% Calculate mean row and column
meanRow = mean(totalBedload_gross.Variables, 1);
meanCol = mean(totalBedload_gross.Variables, 2);

% Add mean row and column to the table
totalBedload_gross{'Mean', :} = meanRow;
totalBedload_gross.Mean = [meanCol; mean(meanRow)];

clearvars Mg d10 d50 d90 Fine Coarse i meanCol meanRow time_seconds


%% Calculation: Total Bedload Transport Volume (Net Longshore)
Mg = NaN(length(instruLocs), 1);
d10 = NaN(length(instruLocs), 1);
d50 = NaN(length(instruLocs), 1);
d90 = NaN(length(instruLocs), 1);
Fine = NaN(length(instruLocs), 1);
Coarse = NaN(length(instruLocs), 1);

Mg_neg = NaN(length(instruLocs), 1);
d10_neg = NaN(length(instruLocs), 1);
d50_neg = NaN(length(instruLocs), 1);
d90_neg = NaN(length(instruLocs), 1);
Fine_neg = NaN(length(instruLocs), 1);
Coarse_neg = NaN(length(instruLocs), 1);

Mg_pos = NaN(length(instruLocs), 1);
d10_pos = NaN(length(instruLocs), 1);
d50_pos = NaN(length(instruLocs), 1);
d90_pos = NaN(length(instruLocs), 1);
Fine_pos = NaN(length(instruLocs), 1);
Coarse_pos = NaN(length(instruLocs), 1);

% Convert datetime to duration and then to seconds
time_seconds = seconds(TT_qL_b10_net_zeros.Time - TT_qL_b10_net_zeros.Time(1));

% Zero out positive values
neg_tt_bMg = TT_qL_bMg_net_zeros;
neg_tt_b10 = TT_qL_b10_net_zeros;
neg_tt_b50 = TT_qL_b50_net_zeros;
neg_tt_b90 = TT_qL_b90_net_zeros;
neg_tt_bFine = TT_qL_bFine_net_zeros;
neg_tt_bCoarse = TT_qL_bCoarse_net_zeros;

% Loop through each variable in the timetable
for varIdx = 1:width(neg_tt_bMg)
    % Get the current variable's data
    data_Mg = neg_tt_bMg{:, varIdx};
    data_10 = neg_tt_b10{:, varIdx};
    data_50 = neg_tt_b50{:, varIdx};
    data_90 = neg_tt_b90{:, varIdx};
    data_Fine = neg_tt_bFine{:, varIdx};
    data_Coarse = neg_tt_bCoarse{:, varIdx};

    % Replace positive or zero values with NaN
    data_Mg(data_Mg >= 0) = 0;
    data_10(data_10 >= 0) = 0;
    data_50(data_50 >= 0) = 0;
    data_90(data_90 >= 0) = 0;
    data_Fine(data_Fine >= 0) = 0;
    data_Coarse(data_Coarse >= 0) = 0;

    % Assign the modified data to the new timetable
    neg_tt_bMg{:, varIdx} = data_Mg;
    neg_tt_b10{:, varIdx} = data_10;
    neg_tt_b50{:, varIdx} = data_50;
    neg_tt_b90{:, varIdx} = data_90;
    neg_tt_bFine{:, varIdx} = data_Fine;
    neg_tt_bCoarse{:, varIdx} = data_Coarse;
end

% Zero out negative values
pos_tt_bMg = TT_qL_bMg_net_zeros;
pos_tt_b10 = TT_qL_b10_net_zeros;
pos_tt_b50 = TT_qL_b50_net_zeros;
pos_tt_b90 = TT_qL_b90_net_zeros;
pos_tt_bFine = TT_qL_bFine_net_zeros;
pos_tt_bCoarse = TT_qL_bCoarse_net_zeros;

% Loop through each variable in the timetable
for varIdx = 1:width(pos_tt_bMg)
    % Get the current variable's data
    data_Mg = pos_tt_bMg{:, varIdx};
    data_10 = pos_tt_b10{:, varIdx};
    data_50 = pos_tt_b50{:, varIdx};
    data_90 = pos_tt_b90{:, varIdx};
    data_Fine = pos_tt_bFine{:, varIdx};
    data_Coarse = pos_tt_bCoarse{:, varIdx};

    % Replace positive or zero values with NaN
    data_Mg(data_Mg < 0) = 0;
    data_10(data_10 < 0) = 0;
    data_50(data_50 < 0) = 0;
    data_90(data_90 < 0) = 0;
    data_Fine(data_Fine < 0) = 0;
    data_Coarse(data_Coarse < 0) = 0;

    % Assign the modified data to the new timetable
    pos_tt_bMg{:, varIdx} = data_Mg;
    pos_tt_b10{:, varIdx} = data_10;
    pos_tt_b50{:, varIdx} = data_50;
    pos_tt_b90{:, varIdx} = data_90;
    pos_tt_bFine{:, varIdx} = data_Fine;
    pos_tt_bCoarse{:, varIdx} = data_Coarse;
end

for i = 1:length(instruLocs)
    
    % Integrate using the trapezoidal rule
    Mg(i) = trapz(time_seconds, TT_qL_bMg_net_zeros.(i));
    d10(i) = trapz(time_seconds, TT_qL_b10_net_zeros.(i));
    d50(i) = trapz(time_seconds, TT_qL_b50_net_zeros.(i));
    d90(i) = trapz(time_seconds, TT_qL_b90_net_zeros.(i));
    Fine(i) = trapz(time_seconds, TT_qL_bFine_net_zeros.(i));
    Coarse(i) = trapz(time_seconds, TT_qL_bCoarse_net_zeros.(i));

    % Display the result
    % fprintf('Total amount of bedload transport: %.2f m^3/m width\n', totalBedload);

    % Perform the integration over negative values
    Mg_neg(i) = trapz(time_seconds, neg_tt_bMg.(i));
    d10_neg(i) = trapz(time_seconds, neg_tt_b10.(i));
    d50_neg(i) = trapz(time_seconds, neg_tt_b50.(i));
    d90_neg(i) = trapz(time_seconds, neg_tt_b90.(i));
    Fine_neg(i) = trapz(time_seconds, neg_tt_bFine.(i));
    Coarse_neg(i) = trapz(time_seconds, neg_tt_bCoarse.(i));

    % Perform the integration over negative values
    Mg_pos(i) = trapz(time_seconds, pos_tt_bMg.(i));
    d10_pos(i) = trapz(time_seconds, pos_tt_b10.(i));
    d50_pos(i) = trapz(time_seconds, pos_tt_b50.(i));
    d90_pos(i) = trapz(time_seconds, pos_tt_b90.(i));
    Fine_pos(i) = trapz(time_seconds, pos_tt_bFine.(i));
    Coarse_pos(i) = trapz(time_seconds, pos_tt_bCoarse.(i));
end

totalBedload_netL = table(Mg, d10, d50, d90, Fine, Coarse, 'RowNames', flipud(instruLocs));
totalBedload_L_neg = table(Mg_neg, d10_neg, d50_neg, d90_neg, Fine_neg, Coarse_neg, 'RowNames', flipud(instruLocs));
totalBedload_L_pos = table(Mg_pos, d10_pos, d50_pos, d90_pos, Fine_pos, Coarse_pos, 'RowNames', flipud(instruLocs));

% Calculate mean row and column
meanRow = mean(totalBedload_netL.Variables, 1);
meanCol = mean(totalBedload_netL.Variables, 2);
meanRow_neg = mean(totalBedload_L_neg.Variables, 1);
meanCol_neg = mean(totalBedload_L_neg.Variables, 2);
meanRow_pos = mean(totalBedload_L_pos.Variables, 1);
meanCol_pos = mean(totalBedload_L_pos.Variables, 2);

% Add mean row and column to the table
totalBedload_netL{'Mean', :} = meanRow;
totalBedload_netL.Mean = [meanCol; mean(meanRow)];
totalBedload_L_neg{'Mean', :} = meanRow_neg;
totalBedload_L_neg.Mean = [meanCol_neg; mean(meanRow_neg)];
totalBedload_L_pos{'Mean', :} = meanRow_pos;
totalBedload_L_pos.Mean = [meanCol_pos; mean(meanRow_pos)];

clearvars Mg d10 d50 d90 Fine Coarse i meanCol meanRow time_seconds Mg_neg d10_neg d50_neg d90_neg Fine_neg Coarse_neg Mg_pos d10_pos d50_pos d90_pos Fine_pos Coarse_pos neg_tt_bMg neg_tt_b10 neg_tt_b50 neg_tt_b90 neg_tt_bFine neg_tt_bCoarse pos_tt_bMg pos_tt_b10 pos_tt_b50 pos_tt_b90 pos_tt_bFine pos_tt_bCoarse data_Mg data_10 data_50 data_90 data_Fine data_Coarse


%% Calculation: Total Bedload Transport Volume (Net Cross-shore)
Mg = NaN(length(instruLocs), 1);
d10 = NaN(length(instruLocs), 1);
d50 = NaN(length(instruLocs), 1);
d90 = NaN(length(instruLocs), 1);
Fine = NaN(length(instruLocs), 1);
Coarse = NaN(length(instruLocs), 1);

Mg_neg = NaN(length(instruLocs), 1);
d10_neg = NaN(length(instruLocs), 1);
d50_neg = NaN(length(instruLocs), 1);
d90_neg = NaN(length(instruLocs), 1);
Fine_neg = NaN(length(instruLocs), 1);
Coarse_neg = NaN(length(instruLocs), 1);

Mg_pos = NaN(length(instruLocs), 1);
d10_pos = NaN(length(instruLocs), 1);
d50_pos = NaN(length(instruLocs), 1);
d90_pos = NaN(length(instruLocs), 1);
Fine_pos = NaN(length(instruLocs), 1);
Coarse_pos = NaN(length(instruLocs), 1);

% Convert datetime to duration and then to seconds
time_seconds = seconds(TT_qC_b10_net_zeros.Time - TT_qC_b10_net_zeros.Time(1));

% Zero out positive values
neg_tt_bMg = TT_qC_bMg_net_zeros;
neg_tt_b10 = TT_qC_b10_net_zeros;
neg_tt_b50 = TT_qC_b50_net_zeros;
neg_tt_b90 = TT_qC_b90_net_zeros;
neg_tt_bFine = TT_qC_bFine_net_zeros;
neg_tt_bCoarse = TT_qC_bCoarse_net_zeros;

% Loop through each variable in the timetable
for varIdx = 1:width(neg_tt_bMg)
    % Get the current variable's data
    data_Mg = neg_tt_bMg{:, varIdx};
    data_10 = neg_tt_b10{:, varIdx};
    data_50 = neg_tt_b50{:, varIdx};
    data_90 = neg_tt_b90{:, varIdx};
    data_Fine = neg_tt_bFine{:, varIdx};
    data_Coarse = neg_tt_bCoarse{:, varIdx};

    % Replace positive or zero values with NaN
    data_Mg(data_Mg >= 0) = 0;
    data_10(data_10 >= 0) = 0;
    data_50(data_50 >= 0) = 0;
    data_90(data_90 >= 0) = 0;
    data_Fine(data_Fine >= 0) = 0;
    data_Coarse(data_Coarse >= 0) = 0;

    % Assign the modified data to the new timetable
    neg_tt_bMg{:, varIdx} = data_Mg;
    neg_tt_b10{:, varIdx} = data_10;
    neg_tt_b50{:, varIdx} = data_50;
    neg_tt_b90{:, varIdx} = data_90;
    neg_tt_bFine{:, varIdx} = data_Fine;
    neg_tt_bCoarse{:, varIdx} = data_Coarse;
end

% Zero out negative values
pos_tt_bMg = TT_qC_bMg_net_zeros;
pos_tt_b10 = TT_qC_b10_net_zeros;
pos_tt_b50 = TT_qC_b50_net_zeros;
pos_tt_b90 = TT_qC_b90_net_zeros;
pos_tt_bFine = TT_qC_bFine_net_zeros;
pos_tt_bCoarse = TT_qC_bCoarse_net_zeros;

% Loop through each variable in the timetable
for varIdx = 1:width(pos_tt_bMg)
    % Get the current variable's data
    data_Mg = pos_tt_bMg{:, varIdx};
    data_10 = pos_tt_b10{:, varIdx};
    data_50 = pos_tt_b50{:, varIdx};
    data_90 = pos_tt_b90{:, varIdx};
    data_Fine = pos_tt_bFine{:, varIdx};
    data_Coarse = pos_tt_bCoarse{:, varIdx};

    % Replace positive or zero values with NaN
    data_Mg(data_Mg < 0) = 0;
    data_10(data_10 < 0) = 0;
    data_50(data_50 < 0) = 0;
    data_90(data_90 < 0) = 0;
    data_Fine(data_Fine < 0) = 0;
    data_Coarse(data_Coarse < 0) = 0;

    % Assign the modified data to the new timetable
    pos_tt_bMg{:, varIdx} = data_Mg;
    pos_tt_b10{:, varIdx} = data_10;
    pos_tt_b50{:, varIdx} = data_50;
    pos_tt_b90{:, varIdx} = data_90;
    pos_tt_bFine{:, varIdx} = data_Fine;
    pos_tt_bCoarse{:, varIdx} = data_Coarse;
end

for i = 1:length(instruLocs)
    
    % Integrate using the trapezoidal rule
    Mg(i) = trapz(time_seconds, TT_qC_bMg_net_zeros.(i));
    d10(i) = trapz(time_seconds, TT_qC_b10_net_zeros.(i));
    d50(i) = trapz(time_seconds, TT_qC_b50_net_zeros.(i));
    d90(i) = trapz(time_seconds, TT_qC_b90_net_zeros.(i));
    Fine(i) = trapz(time_seconds, TT_qC_bFine_net_zeros.(i));
    Coarse(i) = trapz(time_seconds, TT_qC_bCoarse_net_zeros.(i));

    % Display the result
    % fprintf('Total amount of bedload transport: %.2f m^3/m width\n', totalBedload);

    % Perform the integration over negative values
    Mg_neg(i) = trapz(time_seconds, neg_tt_bMg.(i));
    d10_neg(i) = trapz(time_seconds, neg_tt_b10.(i));
    d50_neg(i) = trapz(time_seconds, neg_tt_b50.(i));
    d90_neg(i) = trapz(time_seconds, neg_tt_b90.(i));
    Fine_neg(i) = trapz(time_seconds, neg_tt_bFine.(i));
    Coarse_neg(i) = trapz(time_seconds, neg_tt_bCoarse.(i));

    % Perform the integration over negative values
    Mg_pos(i) = trapz(time_seconds, pos_tt_bMg.(i));
    d10_pos(i) = trapz(time_seconds, pos_tt_b10.(i));
    d50_pos(i) = trapz(time_seconds, pos_tt_b50.(i));
    d90_pos(i) = trapz(time_seconds, pos_tt_b90.(i));
    Fine_pos(i) = trapz(time_seconds, pos_tt_bFine.(i));
    Coarse_pos(i) = trapz(time_seconds, pos_tt_bCoarse.(i));
end

totalBedload_netC = table(Mg, d10, d50, d90, Fine, Coarse, 'RowNames', flipud(instruLocs));
totalBedload_C_neg = table(Mg_neg, d10_neg, d50_neg, d90_neg, Fine_neg, Coarse_neg, 'RowNames', flipud(instruLocs));
totalBedload_C_pos = table(Mg_pos, d10_pos, d50_pos, d90_pos, Fine_pos, Coarse_pos, 'RowNames', flipud(instruLocs));

% Calculate mean row and column
meanRow = mean(totalBedload_netC.Variables, 1);
meanCol = mean(totalBedload_netC.Variables, 2);
meanRow_neg = mean(totalBedload_C_neg.Variables, 1);
meanCol_neg = mean(totalBedload_C_neg.Variables, 2);
meanRow_pos = mean(totalBedload_C_pos.Variables, 1);
meanCol_pos = mean(totalBedload_C_pos.Variables, 2);

% Add mean row and column to the table
totalBedload_netC{'Mean', :} = meanRow;
totalBedload_netC.Mean = [meanCol; mean(meanRow)];
totalBedload_C_neg{'Mean', :} = meanRow_neg;
totalBedload_C_neg.Mean = [meanCol_neg; mean(meanRow_neg)];
totalBedload_C_pos{'Mean', :} = meanRow_pos;
totalBedload_C_pos.Mean = [meanCol_pos; mean(meanRow_pos)];

clearvars varIdx Mg d10 d50 d90 Fine Coarse i meanCol meanRow time_seconds Mg_neg d10_neg d50_neg d90_neg Fine_neg Coarse_neg Mg_pos d10_pos d50_pos d90_pos Fine_pos Coarse_pos neg_tt_bMg neg_tt_b10 neg_tt_b50 neg_tt_b90 neg_tt_bFine neg_tt_bCoarse pos_tt_bMg pos_tt_b10 pos_tt_b50 pos_tt_b90 pos_tt_bFine pos_tt_bCoarse


%% Visualisation: Bed Shear Stress (timeseries) - all measuring
ax = gobjects(1, length(instruLocs));

f3c = figure('Position',[740, 957, 1719, 1336]);
tl = tiledlayout(length(instruLocs), 1, 'TileSpacing','tight');

for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    shearStressData = TT_tau_cw_NaN.(i);

    above_d00 = shearStressData<=GS_stats.tau_cr50_McCarron(i);
    above_d10 = shearStressData>GS_stats.tau_cr10_McCarron(i) & shearStressData<=GS_stats.tau_cr90_McCarron(i);
    above_d50 = shearStressData>GS_stats.tau_cr50_McCarron(i) & shearStressData<=GS_stats.tau_cr10_McCarron(i);
    above_d90 = shearStressData>GS_stats.tau_cr90_McCarron(i);

    tau_d10 = shearStressData;
    tau_d50 = shearStressData;
    tau_d90 = shearStressData;
    tau_d00 = shearStressData;

    tau_d10(~above_d10) = NaN;
    tau_d50(~above_d50) = NaN;
    tau_d90(~above_d90) = NaN;
    tau_d00(~above_d00) = NaN;

    ax(i) = nexttile;
    plot(TT_tau_cw_NaN.Time, tau_d00, '-k', 'LineWidth',2); hold on
    plot(TT_tau_cw_NaN.Time, tau_d50, '-b', 'LineWidth',2)
    plot(TT_tau_cw_NaN.Time, tau_d10, '-r', 'LineWidth',2)
    plot(TT_tau_cw_NaN.Time, tau_d90, '-m', 'LineWidth',2)

    yline(GS_stats.tau_cr10_McCarron(i), '-.', 'LineWidth',2)
    yline(GS_stats.tau_cr50_McCarron(i), ':', 'LineWidth',2)
    yline(GS_stats.tau_cr90_McCarron(i), '--', 'LineWidth',2)

    xticks(dateStart+1:days(2):dateEnd)
    xlim([dateStart, dateEnd])
    ylim([0, 3])

    text(dateEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

    grid on
    box off
end
xticklabels(ax(1:end-1), [])
ylabel(tl, '\tau_{cw} (Pa)', 'FontSize',fontsize)
legend(ax(1), '', '', '', '', '\tau_{cr,10}', '\tau_{cr,50}', '\tau_{cr,90}', 'Location','northoutside', 'NumColumns',3)
zoom xon
linkaxes

clearvars above_d00 above_d10 above_d50 above_d90 tau_d00 tau_d10 tau_d50 tau_d90 tl i ax shearStressData


%% Visualisation: Bed Shear Stress (pdf) - all measuring
% f3d = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_tau_cw_NaN.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Create the KDE plot
%     [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
%         'Support','nonnegative');
%     plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
% end
% hold off
% xlim([0, 3])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Bed Shear Stress')
% xlabel('\tau_{cw} (Pa)')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i


%% Visualisation: Nondimensional Bedload Transport Rate (timeseries) - all measuring
ax = gobjects(1, length(instruLocs));

f4c = figure('Position',[740, 957, 1719, 1336]);
tl = tiledlayout(length(instruLocs), 1, 'TileSpacing','tight');

for i = 1:length(instruLocs)

    ax(i) = nexttile;
    plot(TT_phi_b10_NaN.Time, TT_phi_b10_NaN.(i), '-k', 'LineWidth',2); hold on
    plot(TT_phi_b50_NaN.Time, TT_phi_b50_NaN.(i), '-r', 'LineWidth',2)
    plot(TT_phi_b90_NaN.Time, TT_phi_b90_NaN.(i), '-b', 'LineWidth',2)

    xticks(dateStart+1:days(2):dateEnd)
    xlim([dateStart, dateEnd])
    ylim([0, 3])

    text(dateEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

    grid on
    box off
    
end
xticklabels(ax(1:end-1), [])
ylabel(tl, '\phi_{b}', 'FontSize',fontsize)
legend(ax(1), 'D_{10}', 'D_{50}', 'D_{90}', 'Location','northoutside', 'NumColumns',3)
zoom xon
linkaxes

clearvars ax tl i


%% Visualisation: Nondimensional Bedload Transport Rate (pdf) - all measuring
% f4d = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_phi_b50_NaN.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Create the KDE plot
%     [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
%         'Support','nonnegative');
%     plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
% 
% end
% hold off
% xlim([0, .2])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Nondimensional Bedload Transport Rate')
% xlabel('\phi_{b}')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i


%% Visualisation: Gross Bedload Transport Rate (timeseries) - all measuring
ax = gobjects(1, length(instruLocs));

f5c = figure('Position',[740, 957, 1719, 1336]);
tl = tiledlayout(length(instruLocs), 1, 'TileSpacing','tight');

for i = 1:length(instruLocs)

    ax(i) = nexttile;
    plot(TT_q_b10_zeros.Time, TT_q_b10_zeros.(i), '-k', 'LineWidth',2); hold on
    plot(TT_q_b50_zeros.Time, TT_q_b50_zeros.(i), '-r', 'LineWidth',2)
    plot(TT_q_b90_zeros.Time, TT_q_b90_zeros.(i), '-b', 'LineWidth',2)

    xticks(dateStart+1:days(2):dateEnd)
    xlim([dateStart, dateEnd])
    ylim([0, 8e-5])

    text(dateEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

    grid on
    box off
    
end
xticklabels(ax(1:end-1), [])
ylabel(tl, 'q_{b} (m^2 s^{-1})', 'FontSize',fontsize)
legend(ax(1), 'D_{10}', 'D_{50}', 'D_{90}', 'Location','northoutside', 'NumColumns',3)
zoom xon
linkaxes

clearvars ax tl i


%% Visualisation: Gross Bedload Transport Rate (pdf) - all measuring
% f5d = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_q_b50_zeros.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Create the KDE plot
%     [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
%         'Support','nonnegative');
%     plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
% 
% end
% hold off
% xlim([0, 6e-5])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Bedload Transport')
% xlabel('q_{b} (m^2 s^{-1})')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i


%% Visualisation: Net Bedload Transport Rate (pdf) - all measuring
% p = gobjects(size(instruLocs));
% 
% f5e = figure('Position',[740, 957, 1719, 1336]);
% 
% hold on
% % Loop through each location to create KDE plots
% for i = 1:length(instruLocs)
% 
%     % Extract the variable data for the current location
%     data = TT_q_b90_net_zeros.(i);
% 
%      % Remove NaN values
%     data = data(~isnan(data));
% 
%     % Number of points for higher resolution
%     numPoints = 2e3; % Increase this value for higher resolution
% 
%     % Perform kernel density estimation with higher resolution
%     [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
%     'Support','unbounded', 'NumPoints',numPoints);
% 
%     % Apply Gaussian smoothing
%     gaussianWindowSize = 10; % Adjust window size for more or less smoothing
%     smoothed_f = smoothdata(f, 'gaussian', gaussianWindowSize);
% 
%     % Find the location of the peak frequency (mode) and 
%     [max_f, max_index] = max(smoothed_f);
%     peak_location = xi(max_index);
%     meanValue = mean(data);
% 
%     % Create the KDE plot
%     plot(xi, smoothed_f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
%     xline(peak_location, 'LineStyle','-', 'LineWidth',2, 'Color',newcolors(i,:), 'HandleVisibility','off')
%     xline(meanValue, 'LineStyle','--', 'LineWidth',2, 'Color',newcolors(i,:), 'HandleVisibility','off')
% 
% end
% hold off
% xlim([-.6e-4, .6e-4])
% 
% % Add title, labels and legend
% title('Kernel Density Estimates of Bedload Transport')
% xlabel('< NE        q_{b,50} (m^2 s^{-1})        SW >')
% ylabel('Density')
% legend('show')
% 
% clearvars data xi f i max_f max_index peak_location smoothed_f gaussianWindowSize numPoints meanValue p


%% Visualisation: Bed Shear Stress / Gross & Net BL Transport Rate (pdf)
f6d = figure('Position',[740, 1665, 1719, 628]);
tiledlayout(1,3, "TileSpacing","tight")

nexttile; hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_tau_cw_NaN.(i);

     % Remove NaN values
    data = data(~isnan(data));
    
    % Create the KDE plot
    [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',4, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
axis tight
xlim([0, 2])

% Add title, labels and legend
xlabel('\tau_{cw} (Pa)')
ylabel('Density')
yticklabels({})

nexttile; hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_bMg_zeros.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Create the KDE plot
    [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',4, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
axis tight
xlim([0, 5e-5])

% Add title, labels and legend
xlabel('q_{b,M_G,gross} (m^2 s^{-1})')
legend('show', 'Location','east')
yticklabels({})

nexttile; hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_qL_bMg_net_zeros.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Number of points for higher resolution
    numPoints = 2e3; % Increase this value for higher resolution
    
    % Perform kernel density estimation with higher resolution
    [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
    'Support','unbounded', 'NumPoints',numPoints);

    % Apply Gaussian smoothing
    gaussianWindowSize = 10; % Adjust window size for more or less smoothing
    smoothed_f = smoothdata(f, 'gaussian', gaussianWindowSize);

    % Find the location of the peak frequency (mode)
    [max_f, max_index] = max(smoothed_f);
    peak_location = xi(max_index);
    meanValue = mean(data);

    % Create the KDE plot
    plot(xi, smoothed_f, 'LineWidth',4, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
    % xline(peak_location, 'LineStyle','-', 'LineWidth',4, 'Color',newcolors(i,:), 'HandleVisibility','off')
    % xline(meanValue, 'LineStyle','--', 'LineWidth',4, 'Color',newcolors(i,:), 'HandleVisibility','off')
    % xline(peak_location, 'LineStyle',':', 'LineWidth',2, 'Color','k', 'HandleVisibility','off')

end
hold off
axis tight
xlim([-.4e-4, .4e-4])

% Add title, labels and legend
xlabel('q_{b,M_G,net} (m^2 s^{-1})')
yticklabels({})

% Add figure annotations
annotation('textbox', [0.305, 0.8, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.605, 0.8, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.905, 0.8, 0.1, 0.1], 'String','(c)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.69, 0.6, 0.1, 0.1], 'String',{'SW (ebb)'; 'direction'}, 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');
annotation('textbox', [0.84, 0.6, 0.1, 0.1], 'String',{'NE (flood)'; 'direction'}, 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');

clearvars data xi f i numPoints bw gaussianWindowSize smoothed_f max_f max_index peak_location meanValue


%% Reorganisation
varNames = {'onshore', 'offshore', 'easterly', 'westerly'};
BLV_Mg = table(totalBedload_C_pos.Mg_pos(1:end-1), totalBedload_C_neg.Mg_neg(1:end-1), totalBedload_L_pos.Mg_pos(1:end-1), totalBedload_L_neg.Mg_neg(1:end-1), 'VariableNames',varNames, 'RowNames',instruLocs);
BLV_Mg = flipud(abs(BLV_Mg));

BLV_d10 = table(totalBedload_C_pos.d10_pos(1:end-1), totalBedload_C_neg.d10_neg(1:end-1), totalBedload_L_pos.d10_pos(1:end-1), totalBedload_L_neg.d10_neg(1:end-1), 'VariableNames',varNames, 'RowNames',instruLocs);
BLV_d10 = flipud(abs(BLV_d10));

BLV_d50 = table(totalBedload_C_pos.d50_pos(1:end-1), totalBedload_C_neg.d50_neg(1:end-1), totalBedload_L_pos.d50_pos(1:end-1), totalBedload_L_neg.d50_neg(1:end-1), 'VariableNames',varNames, 'RowNames',instruLocs);
BLV_d50 = flipud(abs(BLV_d50));

BLV_d90 = table(totalBedload_C_pos.d90_pos(1:end-1), totalBedload_C_neg.d90_neg(1:end-1), totalBedload_L_pos.d90_pos(1:end-1), totalBedload_L_neg.d90_neg(1:end-1), 'VariableNames',varNames, 'RowNames',instruLocs);
BLV_d90 = flipud(abs(BLV_d90));

BLV_Fine = table(totalBedload_C_pos.Fine_pos(1:end-1), totalBedload_C_neg.Fine_neg(1:end-1), totalBedload_L_pos.Fine_pos(1:end-1), totalBedload_L_neg.Fine_neg(1:end-1), 'VariableNames',varNames, 'RowNames',instruLocs);
BLV_Fine = flipud(abs(BLV_Fine));

BLV_Coarse = table(totalBedload_C_pos.Coarse_pos(1:end-1), totalBedload_C_neg.Coarse_neg(1:end-1), totalBedload_L_pos.Coarse_pos(1:end-1), totalBedload_L_neg.Coarse_neg(1:end-1), 'VariableNames',varNames, 'RowNames',instruLocs);
BLV_Coarse = flipud(abs(BLV_Coarse));

% Define the variables in a cell array for easy iteration
dataVars = {BLV_Mg, BLV_d10, BLV_d50, BLV_d90, BLV_Fine, BLV_Coarse};

% Convert total Qb during SEDMEX to annual rate
numYears = years(totalDuration_noNaN);
beachWidth = 25;  % average intertidal beach width

newTables = cellfun(@(x) x .* beachWidth ./ numYears, dataVars, 'UniformOutput',false); % m3/y


%% Visualisation: Net BL Transport Rate (bar diagram: m3/m)

% Create a figure and define the layout
f7a = figure('Position',[1281, 957, 1178, 1336]);
tl = tiledlayout(6, 1, 'TileSpacing', 'compact');

% Loop through each variable to create the plots
for i = 1:6
    ax(i) = nexttile;
    bar(table2array(dataVars{i}), 'grouped', 'BarWidth',1, 'GroupWidth',.8); % Create the bar plot
    text(.55, 3.6, GS_headers{i}, 'FontSize',fontsize, 'HorizontalAlignment','left')

    grid off
    ax(i).YGrid = 'on';
    ax(i).GridAlpha = 1;  % Ensure grid lines are fully opaque
    ax(i).GridColor = [0, 0, 0]; % Set grid color (optional, black in this case)
    ax(i).GridLineWidth = 1.5; % Increase grid line width without changing axes
    % ax(i).YColor = 'none'; % Make y-axis invisible
end

% Customize the axes
% set(ax, 'XDir', 'reverse');
set(ax, 'YScale', 'log') % Set y-axis to logarithmic scale
ylim(ax, [1e-4, 10])
yticks(ax, [1e-4, 1e-3, 1e-2, 1e-1, 1, 1e1])
yticklabels(ax, [0.0001, 0.001, 0.01, 0.1, 1, 10])
xticklabels(ax(1:end-1), [])
xticklabels(ax(end), flip(instruLocs))
ylabel(tl, 'total Q_b (m^3 m^{-1})', 'FontSize',fontsize)
legend(ax(1), varNames, 'Location','northoutside', 'NumColumns',4)
box(ax, 'off')


%% Visualisation: Net BL Transport Rate (bar diagram: m3/y)

% Create a figure and define the layout
f7b = figure('Position',[1281, 957, 1178, 1336]);
tl = tiledlayout(6, 1, 'TileSpacing', 'compact');

% Loop through each variable to create the plots
for i = 1:6
    ax(i) = nexttile;
    bar(table2array(newTables{i}), 'grouped', 'BarWidth',1, 'GroupWidth',.8); % Create the bar plot
    text(6.55, 4e3, GS_headers{i}, 'FontSize',fontsize, 'HorizontalAlignment','right')

    grid off
    ax(i).YGrid = 'on';
    ax(i).GridAlpha = 1;  % Ensure grid lines are fully opaque
    ax(i).GridColor = [0, 0, 0]; % Set grid color (optional, black in this case)
    ax(i).GridLineWidth = 1.5; % Increase grid line width without changing axes
    % ax(i).YColor = 'none'; % Make y-axis invisible
end

% Customize the axes
% set(ax, 'XDir', 'reverse');
set(ax, 'YScale', 'log') % Set y-axis to logarithmic scale
ylim(ax, [1, 1e4])
yticks(ax, [1e0, 1e1, 1e2, 1e3, 1e4])
yticklabels(ax, [1, 10, 100, 1000, 10000])
xticklabels(ax(1:end-1), [])
xticklabels(ax(end), flip(instruLocs))
ylabel(tl, 'annual Q_b (m^3 y^{-1})', 'FontSize',fontsize)
legend(ax(1), varNames, 'Location','northoutside', 'NumColumns',4)
box(ax, 'off')

colororder(brewermap(4,'Set3'))


%% CHECK

f = figure('Position',[740, 957, 1719, 1336]); hold on
plot(TT_q_bMg_zeros.Time, TT_q_bMg_zeros.L1C1VEC, 'k', 'LineWidth',3)
plot(TT_q_bMg_zeros.Time, TT_qL_bMg_net_zeros.L1C1VEC, 'r', 'LineWidth',3)
plot(TT_q_bMg_zeros.Time, TT_qC_bMg_net_zeros.L1C1VEC, 'b', 'LineWidth',3)

%%
f = figure('Position',[740, 957, 1719, 1336]); hold on
plot(t, uLong_z, 'k', 'LineWidth',3)
plot(t, uCross_z, 'r', 'LineWidth',3)

%%
U = sqrt( uLong_z.^2 + uCross_z.^2 );