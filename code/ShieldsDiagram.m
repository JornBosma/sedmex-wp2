%% Initialization
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

dataPath = fullfile(filesep, 'Volumes', 'T7 Shield', 'DataDescriptor', 'hydrodynamics', 'ADV', filesep);

instruments = {'L1C1VEC', 'L2C5SONTEK1', 'L3C1VEC', 'L4C1VEC', 'L5C1VEC', 'L6C1VEC'};
sampleLocs = {'L0', 'L2', 'L3.5', 'L4', 'Tmb', 'L6'};
instruLocs = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'};
GS_fractions = {'Mg', 'd10', 'd50', 'd90', 'Fine', 'Coarse'};

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
% Coarse = repmat(30*1e-3, length(Mg), 1); % when waves become indispensible

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
load('hab_measured.mat');

% Iterate over each variable in the timetable
for varName = hab.Properties.VariableNames(1:end)
    % Get the variable data
    data = hab.(varName{1})./100; % [cm] to [m]
    
    % Apply the condition: if value < 0.1, set it to 0.1 (to prevent
    % extremely high depth-average velocities)
    data(data < 0.12) = 0.12;
    
    % Update the variable in the timetable
    hab.(varName{1}) = data;
end

clearvars varName data


%% Computations

% Defined starting time
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');

% Create a datetime array for the time axis
dateStart = datetime(2021, 9, 10);
dateEnd = datetime(2021, 10, 19);
dt = minutes(10);
timeAxis = dateStart:dt:dateEnd;
timeAxis = timeAxis';

% Preallocate the data matrix with NaNs
n = length(timeAxis);
emptyData = nan(n, length(instruments));

% List of variable suffixes
suffixes = {'tau_c', 'tau_w', 'tau_cw', ...
    'theta_c10', 'theta_w10', 'theta_cw10', 'phi_b10', 'q_b10', ...
    'theta_c50', 'theta_w50', 'theta_cw50', 'phi_b50', 'q_b50', ...
    'theta_c90', 'theta_w90', 'theta_cw90', 'phi_b90', 'q_b90', ...
    'theta_cMg', 'theta_wMg', 'theta_cwMg', 'phi_bMg', 'q_bMg', ...
    'theta_cFine', 'theta_wFine', 'theta_cwFine', 'phi_bFine', 'q_bFine', ...
    'theta_cCoarse', 'theta_wCoarse', 'theta_cwCoarse', 'phi_bCoarse', 'q_bCoarse', ...
    'long_sign', 'long_frac','cross_sign', 'cross_frac', 'U'};

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

    uLong_z = ncread(filename, 'ulm');  % longshore flow velocity [m/s] at depth z
    uLong_sign = sign(uLong_z);         % direction sign of longshore current
    % Positive onshore direction at this cross-section is 135◦
    % counterclockwise from east, and positive alongshore direction is
    % 225◦ counterclockwise from east.
    uLong_sign = -uLong_sign;           % reverse the sign so that NE (flood) is positive
    uLong_sign(isnan(uLong_sign)) = 0; % NaN -> 0

    uCross_z = ncread(filename, 'ucm');  % cross-shore flow velocity [m/s] at depth z
    uCross_sign = sign(uCross_z);        % direction sign of cross-shore current
    uCross_sign(isnan(uCross_sign)) = 0; % NaN -> 0

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

    % Correct SONTEK wave height based on values L2C4ADV and L2C6OSSI
    if strcmp(instruments{i}, 'L2C5SONTEK1')
        [~, H] = correct_sontek(t, H, 'Hm0');
    end

    % Interpolate the measurements to the new time vector
    time = t0 + seconds(t);
    habCVnew = retime(hab, time, 'pchip');
    z = habCVnew.(instruments{i});

    % Estimate the depth-averaged current velocity
    k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)
    [u_c, ~] = compute_DAV(abs(u_z), z, k_sc, h);

    % Calculate the longshore and cross-shore fractions of the total velocity vector
    uLong_frac = abs(uLong_z).^2 ./ (abs(uLong_z).^2 + abs(uCross_z).^2);
    % uLong_frac(isnan(uLong_frac)) = 0; % NaN -> 0

    uCross_frac = abs(uCross_z).^2 ./ (abs(uLong_z).^2 + abs(uCross_z).^2);
    % uCross_frac(isnan(uCross_frac)) = 0; % NaN -> 0

    % Compute the shear stress components
    [tau_c, tau_w, tau_cw] = compute_BSS_orbital(u_c, h, Urms, T, rho_w, phi_c, phi_w, GS_stats.d50(i), g);
    % [tau_c, tau_w, tau_cw] = compute_BSS_linearTheory(u_c, H, k, T, h, rho_w, phi_c, phi_w, GS_stats.d50(i), g);

    % Append to timetable
    firstRow = find(TT_tau_c.Time == tStart);
    lastRow = find(TT_tau_c.Time == tEnd);

    TT_U.(i)(firstRow:lastRow) = u_c;

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

end

clearvars dataPath dt emptyData filename firstRow h H i info k lastRow n phi_c phi_w t T tau_c tau_w tau_cw theta_c10 theta_c50 theta_c90 theta_cFine theta_cCoarse theta_w10 theta_w50 theta_w90 theta_wFine theta_wCoarse theta_cw10 theta_cw50 theta_cw90 theta_cwFine theta_cwCoarse timeAxis rho_w u_z u_c h Urms emptyData firstRow lastRow z z_new z_measured nonNaNValues nonNaNIndices interpolatedValues phi_b10 phi_b50 phi_b90 phi_bFine phi_bCoarse q_b10 q_b50 q_b90 q_bFine q_bCoarse uLong_z uLong_sign uLong_frac uCross_z uCross_sign uCross_frac q_b10 q_b50 q_b90 q_bFine q_bCoarse suffix suffixes timetableName tStart tEnd time t0 habCVnew theta_cMg theta_wMg theta_cwMg phi_bMg q_bMg q_bMg locations_c locations_w qL_bMg qL_b10 qL_b50 qL_b90 qL_bFine qL_bCoarse qC_bMg qC_b10 qC_b50 qC_b90 qC_bFine qC_bCoarse


%% Visualisation: BSS timeseries overview
f1 = figure('Position',[988, 1665, 1301, 628]);
tiledlayout(3,1, 'TileSpacing','tight')

nexttile
hold on
for i = 1:width(TT_tau_c)
    plot(TT_tau_c.Time, TT_tau_c.(i), 'LineWidth',2)
end
hold off

yline(mean(GS_stats.tau_cr10_McCarron), '-.', 'LineWidth',2)
yline(mean(GS_stats.tau_cr50_McCarron), ':', 'LineWidth',2)
yline(mean(GS_stats.tau_cr90_McCarron), '--', 'LineWidth',2)
 
legend([instruLocs, '\tau_{cr,10}', '\tau_{cr,50}', '\tau_{cr,90}'], 'Box','on', 'Location','northoutside', 'NumColumns',9)
xticks(dateStart+1:days(4):dateEnd)
xticklabels({})
ylabel('\tau_{c} (Pa)')
xlim([dateStart, dateEnd])
ylim([0, 3])
grid on

nexttile
hold on
for i = 1:width(TT_tau_w)
    plot(TT_tau_w.Time, TT_tau_w.(i), 'LineWidth',2)
end
hold off

yline(mean(GS_stats.tau_cr10_McCarron), '-.', 'LineWidth',2)
yline(mean(GS_stats.tau_cr50_McCarron), ':', 'LineWidth',2)
yline(mean(GS_stats.tau_cr90_McCarron), '--', 'LineWidth',2)

% text(tEnd+hours(5), GS_stats.tau_cr10_McCarron(i), '\tau_{cr,10}', 'FontSize',fontsize*.8)
% text(tEnd+hours(5), GS_stats.tau_cr50_McCarron(i), '\tau_{cr,50}', 'FontSize',fontsize*.8)
% text(tEnd+hours(5), GS_stats.tau_cr90_McCarron(i), '\tau_{cr,90}', 'FontSize',fontsize*.8)

% lgnd = legend('', '', '', '', '', '', '\tau_{cr,10}', '\tau_{cr,50}', '\tau_{cr,90}', 'Location','northeast', 'NumColumns',3);

xticks(dateStart+1:days(4):dateEnd)
xticklabels({})
ylabel('\tau_{w} (Pa)')
xlim([dateStart, dateEnd])
ylim([0, 3])
grid on

nexttile
hold on
for i = 1:length(instruLocs)
    plot(TT_tau_cw.Time, TT_tau_cw.(i), 'LineWidth',2)
end
hold off

yline(mean(GS_stats.tau_cr10_McCarron), '-.', 'LineWidth',2)
yline(mean(GS_stats.tau_cr50_McCarron), ':', 'LineWidth',2)
yline(mean(GS_stats.tau_cr90_McCarron), '--', 'LineWidth',2)

% text(tEnd+hours(5), GS_stats.tau_cr10_McCarron(i), '\tau_{cr,10}', 'FontSize',fontsize*.8)
% text(tEnd+hours(5), GS_stats.tau_cr50_McCarron(i), '\tau_{cr,50}', 'FontSize',fontsize*.8)
% text(tEnd+hours(5), GS_stats.tau_cr90_McCarron(i), '\tau_{cr,90}', 'FontSize',fontsize*.8)

colororder(f1, newcolors)
xticks(dateStart+1:days(4):dateEnd)
ylabel('\tau_{cw} (Pa)')
xlim([dateStart, dateEnd])
ylim([0, 3])
zoom xon
linkaxes
grid on

clearvars i


%% Calculations: relative importance current and waves (1/2)
threshold = GS_stats.theta_crMg_McCarron;
mixedBoundary = threshold./3;
thresholdMean = mean(threshold);
mixedBoundaryMean = thresholdMean/3;

data_c = TT_theta_cMg{:,:};
data_w = TT_theta_wMg{:,:};
data_cw = TT_theta_cwMg{:,:};

% Find indices within the date range
idx = (TT_theta_cMg.Time >= stormy_period(1)) & (TT_theta_cMg.Time <= stormy_period(2));
% idx = (TT_theta_cMg.Time >= storm_period(1)) & (TT_theta_cMg.Time <= storm_period(2));

% Get the first and last index
firstIdx = find(idx, 1, 'first');
lastIdx = find(idx, 1, 'last');

nonStormIdx = [1:firstIdx-1, lastIdx+1:height(TT_theta_cMg)];
stormIdx = firstIdx:lastIdx;

% % Count the number of non-NaN elements in each column
% elementCounts_c = sum(isfinite(data_c), 1);
% elementCounts_w = sum(isfinite(data_w), 1);
% 
% fractionExceed_c = sum(data_c>meanThreshold) ./ elementCounts_c * 100;
% fractionExceed_w = sum(data_w>meanThreshold) ./ elementCounts_w * 100;

% Create axis for combined current and wave-related Shields
x_cw = linspace(1e-6, 1, 1e4);

clearvars elementCounts_c elementCounts_w fractionExceed_c fractionExceed_w idx firstIdx lastIdx


%% Calculations: relative importance current and waves (BULK)

% Flatten the arrays
data_c_flat = data_c(:);
data_w_flat = data_w(:);

% Create a mask for non-NaN elements and apply to arrays
nonNaN_mask = ~isnan(data_c_flat) & ~isnan(data_w_flat);
num_nonNaN = sum(nonNaN_mask);
data_c_flat = data_c_flat(nonNaN_mask);
data_w_flat = data_w_flat(nonNaN_mask);

% Calculate the sum of data_c and data_w
sum_data = data_c_flat + data_w_flat;

% Condition 1: waves only
cond1 = sum_data > thresholdMean & data_c_flat <= mixedBoundaryMean;
percent_cond1 = (sum(cond1) / num_nonNaN) * 100;

% Condition 2: wave-dominated
cond2 = sum_data > thresholdMean & data_w_flat > data_c_flat & data_c_flat > mixedBoundaryMean;
percent_cond2 = (sum(cond2) / num_nonNaN) * 100;

% Condition 3: current-dominated
cond3 = sum_data > thresholdMean & data_c_flat > data_w_flat & data_w_flat > mixedBoundaryMean;
percent_cond3 = (sum(cond3) / num_nonNaN) * 100;

% Condition 4: current only
cond4 = sum_data > thresholdMean & data_w_flat <= mixedBoundaryMean;
percent_cond4 = (sum(cond4) / num_nonNaN) * 100;

% Condition 5: no motion
cond5 = sum_data <= thresholdMean;
percent_cond5 = (sum(cond5) / num_nonNaN) * 100;

percent_tot = percent_cond1+percent_cond2+percent_cond3+percent_cond4+percent_cond5;

% Display the results
fprintf('Percentage wave only: %.2f%%\n', percent_cond1);
fprintf('Percentage wave-dominated: %.2f%%\n', percent_cond2);
fprintf('Percentage current-dominated: %.2f%%\n', percent_cond3);
fprintf('Percentage current only: %.2f%%\n', percent_cond4);
fprintf('Percentage no motion: %.2f%%\n', percent_cond5);
fprintf('Total percentage: %.2f%%\n', percent_tot);

clearvars data_c_flat data_w_flat sum_data non_nan_mask cond1 cond2 cond3 cond4 cond5 percent_tot


%% Visualisation: relative importance current and waves (BULK)
f2a = figure('Position',[740, 957, 1719, 628]);

% Plot each column with different color and symbol
hold on
for i = 1:size(data_c, 2)
    scatter(data_c(:, i), data_w(:, i), 20, 'Marker',newsymbols{i},...
        'MarkerEdgeColor',newcolors(i,:), 'LineWidth',1);
end

% Plot threshold lines
plot(x_cw, thresholdMean - x_cw, '-k', 'LineWidth',2); % incipient motion
plot(x_cw, 1 - x_cw, '-k', 'LineWidth',2); % sheet flow (Kleinhans, 2002)
a = 1;
b = mixedBoundaryMean;
c = thresholdMean-b;
d = thresholdMean/2;
line([a-.5, d], [a-.5, d], 'Color','k', 'LineStyle','--', 'LineWidth',3)
line([a, c], [b, b], 'Color','k', 'LineStyle','--', 'LineWidth',3)
line([b, b], [a, c], 'Color','k', 'LineStyle','--', 'LineWidth',3)

% Add region descriptions
text(4e-4, .6, 'waves only', 'FontSize',fontsize*.6)
text(.014, .6, 'wave-dominated', 'FontSize',fontsize*.6)
text(.6, .25, 'current-dominated', 'FontSize',fontsize*.6, 'Rotation',-90)
text(.6, 2e-3, 'current only', 'FontSize',fontsize*.6, 'Rotation',-90)
text(4e-4, 2e-4, 'no motion', 'FontSize',fontsize*.6)

% Add region fractions
text(8e-4, .2, [mat2str(round(percent_cond1), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(.03, .2, [mat2str(round(percent_cond2), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(.1, .03, [mat2str(round(percent_cond3), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(.1, 1e-3, [mat2str(round(percent_cond4), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(8e-4, 1e-3, [mat2str(round(percent_cond5), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')

% Additional plotting commands
% plot(x_cw, mean_threshold - x_cw, '--w', 'LineWidth', 2)
% text(mean_threshold - .027, mean_threshold - .019, '\theta_{cr}', 'FontSize', fontsize, 'Color', 'w')
% xline(mean_threshold, '-', '\theta_{cr}', 'LineWidth', 2, 'FontSize', fontsize)
% yline(mean_threshold, '-', '\theta_{cr}', 'LineWidth', 2, 'FontSize', fontsize)
hold off

% Set axes properties
set(gca, 'XScale','log');
set(gca, 'YScale','log');

% Label axes
xlabel('\theta_{c}')
ylabel('\theta_{w}')

% Add legend
legend(instruLocs, 'Location','eastoutside')

% Set axis limits
xlim([1e-4, 1])
ylim([1e-4, 1])

% Enable grid and set axis properties
zoom xon
linkaxes
grid on
axis square

clearvars a b c d percent_cond1 percent_cond2 percent_cond3 percent_cond4 percent_cond5


%% Calculations: relative importance current and waves (2/2)

% Initialize matrices to store percentages for each condition
nCols = width(data_c);
percent_cond1 = zeros(1, nCols);
percent_cond2 = zeros(1, nCols);
percent_cond3 = zeros(1, nCols);
percent_cond4 = zeros(1, nCols);
percent_cond5 = zeros(1, nCols);

for col = 1:nCols
    % Extract the current column
    data_c_col = data_c(:, col);
    data_w_col = data_w(:, col);
    
    % Create a mask for non-NaN elements and apply to arrays
    nonNaN_mask = ~isnan(data_c_col) & ~isnan(data_w_col);
    num_nonNaN = sum(nonNaN_mask);
    data_c_col = data_c_col(nonNaN_mask);
    data_w_col = data_w_col(nonNaN_mask);
    
    % Calculate the sum of data_c and data_w
    sum_data = data_c_col + data_w_col;
    
    % Condition 1: waves only
    cond1 = sum_data > threshold(col) & data_c_col <= mixedBoundary(col);
    percent_cond1(col) = (sum(cond1) / num_nonNaN) * 100;
    
    % Condition 2: wave-dominated
    cond2 = sum_data > threshold(col) & data_w_col > data_c_col & data_c_col > mixedBoundary(col);
    percent_cond2(col) = (sum(cond2) / num_nonNaN) * 100;
    
    % Condition 3: current-dominated
    cond3 = sum_data > threshold(col) & data_c_col > data_w_col & data_w_col > mixedBoundary(col);
    percent_cond3(col) = (sum(cond3) / num_nonNaN) * 100;
    
    % Condition 4: current only
    cond4 = sum_data > threshold(col) & data_w_col <= mixedBoundary(col);
    percent_cond4(col) = (sum(cond4) / num_nonNaN) * 100;
    
    % Condition 5: no motion
    cond5 = sum_data <= threshold(col);
    percent_cond5(col) = (sum(cond5) / num_nonNaN) * 100;

end

percent_tot = percent_cond1 + percent_cond2 + percent_cond3 + percent_cond4 + percent_cond5;

% Create table
colNames = {'noMotion', 'currentOnly', 'currentDom', 'waveDom', 'waveOnly', 'total'};
percentTable = table(percent_cond5', percent_cond4', percent_cond3', percent_cond2', percent_cond1', percent_tot', 'RowNames', flipud(instruLocs), 'VariableNames',colNames);

% Calculate mean row and column
meanRow = mean(percentTable.Variables, 1);

% Add mean row and column to the table
percentTable{'Mean', :} = meanRow;

clearvars data_c_col data_w_col sum_data nonNaN_mask num_nonNaN cond1 cond2 cond3 cond4 cond5 percent_tot nCols col meanRow


%% Visualisation: relative importance current and waves
f2b = figure('Position',[740, 1665, 1719, 628]);
tiledlayout(1,3, 'TileSpacing','compact')

j = 1;
ax = gobjects(1,3);
for i = [1, 4, 6]
    ax(j) = nexttile; hold on
    scatter(data_c(:, i), data_w(:, i), 50, 'Marker','o',...
        'MarkerEdgeColor',newcolors(i,:), 'LineWidth',1);
    % scatter(data_c(nonStormIdx, i), data_w(nonStormIdx, i), 50, 'Marker','o',...
    %     'MarkerEdgeColor',newcolors(i,:), 'LineWidth',1);
    % scatter(data_c(stormIdx, i), data_w(stormIdx, i), 50, 'Marker','+',...
    %     'MarkerFaceColor','none', 'MarkerEdgeColor','r');
    text(.47, .7, instruLocs{i}, 'FontSize',fontsize, 'FontWeight','bold', 'Color','k')

    % Plot threshold lines
    plot(x_cw, threshold(i) - x_cw, '-k', 'LineWidth',2); % incipient motion
    plot(x_cw, 1 - x_cw, '-k', 'LineWidth',2); % sheet flow (Kleinhans, 2002)
    a = 1;
    b = mixedBoundary(i);
    c = threshold(i)-b;
    d = threshold(i)/2;
    line([a-.5, d], [a-.5, d], 'Color','k', 'LineStyle','--', 'LineWidth',3)
    line([a, c], [b, b], 'Color','k', 'LineStyle','--', 'LineWidth',3)
    line([b, b], [a, c], 'Color','k', 'LineStyle','--', 'LineWidth',3)
    hold off

    % Add region descriptions
    text(4e-4, .6, 'waves only', 'FontSize',fontsize*.6)
    text(.014, .6, 'wave-dominated', 'FontSize',fontsize*.6)
    text(.6, .25, 'current-dominated', 'FontSize',fontsize*.6, 'Rotation',-90)
    text(.6, 2e-3, 'current only', 'FontSize',fontsize*.6, 'Rotation',-90)
    text(4e-4, 2e-4, 'no motion', 'FontSize',fontsize*.6)
    
    % Add region fractions
    text(8e-4, .2, [mat2str(round(percent_cond1(i)), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(.03, .2, [mat2str(round(percent_cond2(i)), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(.1, .03, [mat2str(round(percent_cond3(i)), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(.1, 1e-3, [mat2str(round(percent_cond4(i)), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(8e-4, 1e-3, [mat2str(round(percent_cond5(i)), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
        
    % Set axes properties
    set(gca, 'XScale','log');
    set(gca, 'YScale','log');
    
    % Set axis limits
    xlim([1e-4, 1])
    ylim([1e-4, 1])

    % Enable grid and set axis properties
    zoom xon
    linkaxes
    grid on
    axis square

    j = j+1;
end

% Set axis properties
xlabel(ax(2), '\theta_{c}')
ylabel(ax(1), '\theta_{w}')
yticklabels(ax(2:3), {})

clearvars data_c data_w x_cw a b c d ax j percent_cond1 percent_cond2 percent_cond3 percent_cond4 percent_cond5