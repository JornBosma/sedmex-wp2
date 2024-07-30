%% Initialization
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

dataPath = fullfile(filesep, 'Volumes', 'T7 Shield', 'DataDescriptor', 'hydrodynamics', 'ADV', filesep);

instruments = ["L1C1VEC", "L2C4VEC", "L3C1VEC", "L4C1VEC", "L5C1VEC", "L6C1VEC"];
sampleLocs = {'L1', 'L2', 'L3.5', 'L4', 'Tmb', 'L6'};
instruLocs = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'};
GS_fractions = {'Mg', 'd10', 'd50', 'd90'};

% newcolors = [cbf.vermilion; cbf.blue; cbf.bluegreen; cbf.yellow; cbf.redpurp; cbf.skyblue; cbf.orange];
newcolors = cbf.six([2, 1, 3:5, 7:end], :);
newsymbols = {'o', 's', '^', 'd', 'v', 'p'};

g = 9.81;     
rho_s = 2650;  
rho_w = 1025; 


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

GS_stats = table(Mg, d10, d50, d90, fg, 'RowNames', flipud(sampleLocs));

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
percentiles = {'Mg', '10', '50', '90'};

% Loop over prefixes, methods, and percentiles to assign values
for p = 1:length(prefixes)
    for j = 1:num_methods
        method = crit_methods{j};
        for perc = 1:length(percentiles)
            fieldName = sprintf('%s%s_%s', prefixes{p}, percentiles{perc}, method);
            GS_stats.(fieldName) = eval(sprintf('%s(:, %d, %d)', prefixes{p}, perc, j));
        end
    end
end

clearvars d10 d50 d90 fg theta_cr tau_cr crit_methods prefixes percentiles i j p perc num_methods num_percentiles method fieldName Mg Mg_a Mg_b d10_a d10_b d50_a d50_b d90_a d90_b fg_a fg_b


%% Correct measured height of control volume above bed
hab_measured = load('hab_C4.mat');
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


%% Computations: bed-shear stress & Shields
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
suffixes = {'tau_c', 'tau_w', 'tau_cw', 'theta_c10', 'theta_w10', 'theta_cw10', ...
            'phi_b10', 'q_b10', 'theta_c50', 'theta_w50', 'theta_cw50', ...
            'phi_b50', 'q_b50', 'theta_c90', 'theta_w90', 'theta_cw90', ...
            'phi_b90', 'q_b90', 'long_sign', 'q_b10_net', 'q_b50_net', 'q_b90_net', ...
            'theta_cMg', 'theta_wMg', 'theta_cwMg', 'phi_bMg', 'q_bMg', 'q_bMg_net', ...
            'long_frac'};

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

    % Calculate the longshore fraction of the total velocity vector
    uLong_frac = abs(uLong_z) ./ u_z;

    % Estimate the depth-averaged current velocity
    k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)
    [u_c, ~] = compute_DAV(abs(u_z), z_new, k_sc, h);

    % Compute the shear stress components
    [tau_c, tau_w, tau_cw] = compute_BSS_orbital(u_c, h, Urms, T, rho_w, phi_c, phi_w, GS_stats.d90(i), g);
    % [tau_c, tau_w, tau_cw] = compute_BSS_linearTheory(u_c, H, k, T, h, rho_w, phi_c, phi_w, GS_stats.d90(i), g);

    % Append to timetable
    firstRow = find(TT_tau_c.Time == tStart);
    lastRow = find(TT_tau_c.Time == tEnd);

    TT_long_sign.(i)(firstRow:lastRow) = uLong_sign;
    TT_long_frac.(i)(firstRow:lastRow) = uLong_frac;

    TT_tau_c.(i)(firstRow:lastRow) = tau_c;
    TT_tau_w.(i)(firstRow:lastRow) = tau_w;
    TT_tau_cw.(i)(firstRow:lastRow) = tau_cw;
    
    % Fraction-specific Shields numbers
    [theta_cMg, theta_wMg, theta_cwMg] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.Mg(i), g);
    [theta_c10, theta_w10, theta_cw10] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d10(i), g);
    [theta_c50, theta_w50, theta_cw50] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d50(i), g);
    [theta_c90, theta_w90, theta_cw90] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d90(i), g);

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

    % Nondimensional bedload predictors (gross)
    alpha = 11;  % calibration coefficient of Ribberink (1998)
    beta = 1.65; % calibration exponent

    phi_bMg = compute_Einstein_parameter(TT_theta_cwMg.(i), GS_stats.theta_crMg_McCarron(i), alpha, beta);
    phi_b10 = compute_Einstein_parameter(TT_theta_cw10.(i), GS_stats.theta_cr10_McCarron(i), alpha, beta);
    phi_b50 = compute_Einstein_parameter(TT_theta_cw50.(i), GS_stats.theta_cr50_McCarron(i), alpha, beta);
    phi_b90 = compute_Einstein_parameter(TT_theta_cw90.(i), GS_stats.theta_cr90_McCarron(i), alpha, beta);

    TT_phi_bMg.(i) = phi_bMg;
    TT_phi_b10.(i) = phi_b10;
    TT_phi_b50.(i) = phi_b50;
    TT_phi_b90.(i) = phi_b90;

    % Dimensional bedload transport rate (gross)
    q_bMg = compute_transport_rate(phi_bMg, rho_w, rho_s, GS_stats.Mg(i), g);
    q_b10 = compute_transport_rate(phi_b10, rho_w, rho_s, GS_stats.d10(i), g);
    q_b50 = compute_transport_rate(phi_b50, rho_w, rho_s, GS_stats.d50(i), g);
    q_b90 = compute_transport_rate(phi_b90, rho_w, rho_s, GS_stats.d90(i), g);

    TT_q_bMg.(i) = q_bMg;
    TT_q_b10.(i) = q_b10;
    TT_q_b50.(i) = q_b50;
    TT_q_b90.(i) = q_b90;

    % Dimensional bedload transport rate (net: longshore direction)
    q_bMg_net = TT_long_sign.(i) .* q_bMg .* TT_long_frac.(i);
    q_b10_net = TT_long_sign.(i) .* q_b10 .* TT_long_frac.(i);
    q_b50_net = TT_long_sign.(i) .* q_b50 .* TT_long_frac.(i);
    q_b90_net = TT_long_sign.(i) .* q_b90 .* TT_long_frac.(i);

    TT_q_bMg_net.(i) = q_bMg_net;
    TT_q_b10_net.(i) = q_b10_net;
    TT_q_b50_net.(i) = q_b50_net;
    TT_q_b90_net.(i) = q_b90_net;

end

clearvars dataPath dt emptyData filename firstRow h H i info k lastRow n phi_c phi_w t T tau_c tau_w tau_cw theta_c10 theta_c50 theta_c90 theta_w10 theta_w50 theta_w90 theta_cw10 theta_cw50 theta_cw90 timeAxis rho_w u_z u_c h Urms emptyData firstRow lastRow z z_new z_measured nonNaNValues nonNaNIndices interpolatedValues phi_b10 phi_b50 phi_b90 q_b10 q_b50 q_b90 uLong_z uLong_sign uLong_frac q_b10_net q_b50_net q_b90_net suffix suffixes timetableName tStart tEnd time t0 habCVnew theta_cMg theta_wMg theta_cwMg phi_bMg q_bMg q_bMg_net


%% Quality check
variableNames = TT_tau_w.Properties.VariableNames;

% Iterate over each ADV
for i = 1:numel(variableNames)
    % Find locations of BSS spikes
    locations_c = find(TT_tau_c.(variableNames{i}) > 3.0);
    locations_w = find(TT_tau_w.(variableNames{i}) > 3.5);

    % Change the values to NaN
    TT_tau_c.(variableNames{i})(locations_c) = NaN;
    TT_tau_w.(variableNames{i})(locations_w) = NaN;
    TT_tau_cw.(variableNames{i})(locations_c) = NaN;
    TT_tau_cw.(variableNames{i})(locations_w) = NaN;
end

clearvars variableNames i locations_c locations_w


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


%% Calculations: relative importance current and waves (1/3)

% threshold_50_Soulsby = mean(GS_stats.theta_cr50_Soulsby);
% threshold_50_Egiazaroff = mean(GS_stats.theta_cr50_Egiazaroff);
threshold_50_McCarron = mean(GS_stats.theta_cr50_McCarron);

threshold = threshold_50_McCarron;

data_c = TT_theta_c50{:,:};
data_w = TT_theta_w50{:,:};

% Count the number of non-NaN elements in each column
elementCounts_c = sum(isfinite(data_c), 1);
elementCounts_w = sum(isfinite(data_w), 1);

fractionExceed_c = sum(data_c>threshold) ./ elementCounts_c * 100;
fractionExceed_w = sum(data_w>threshold) ./ elementCounts_w * 100;

% Create axis for combined current and wave-related Shields
x_cw = linspace(1e-6, 1, 1e4);

clearvars threshold_50_McCarron elementCounts_c elementCounts_w fractionExceed_c fractionExceed_w


%% Calculations: relative importance current and waves (2/3)
% Flatten the arrays
data_c_flat = data_c(:);
data_w_flat = data_w(:);

% Calculate the sum of data_c and data_w
sum_data = data_c_flat + data_w_flat;

% Create a mask for non-NaN elements
non_nan_mask = ~isnan(data_c_flat) & ~isnan(data_w_flat);

% Condition 1: waves only
cond1 = sum_data > threshold & data_c_flat < 0.01 & non_nan_mask;
percentage_cond1 = (sum(cond1) / sum(non_nan_mask)) * 100;

% Condition 2: wave-dominated
cond2 = sum_data > threshold & data_w_flat > data_c_flat & data_c_flat > 0.01 & data_w_flat > 0.01 & non_nan_mask;
percentage_cond2 = (sum(cond2) / sum(non_nan_mask)) * 100;

% Condition 3: current-dominated
cond3 = sum_data > threshold & data_c_flat > data_w_flat & data_c_flat > 0.01 & data_w_flat > 0.01 & non_nan_mask;
percentage_cond3 = (sum(cond3) / sum(non_nan_mask)) * 100;

% Condition 4: current only
cond4 = sum_data > threshold & data_w_flat < 0.01 & non_nan_mask;
percentage_cond4 = (sum(cond4) / sum(non_nan_mask)) * 100;

% Condition 5: no motion
cond5 = sum_data <= threshold & non_nan_mask;
percentage_cond5 = (sum(cond5) / sum(non_nan_mask)) * 100;

percentage_tot = percentage_cond1+percentage_cond2+percentage_cond3+percentage_cond4+percentage_cond5;

% Display the results
fprintf('Percentage wave only: %.2f%%\n', percentage_cond1);
fprintf('Percentage wave-dominated: %.2f%%\n', percentage_cond2);
fprintf('Percentage current-dominated: %.2f%%\n', percentage_cond3);
fprintf('Percentage current only: %.2f%%\n', percentage_cond4);
fprintf('Percentage no motion: %.2f%%\n', percentage_cond5);
fprintf('Total percentage: %.2f%%\n', percentage_tot);

clearvars data_c_flat data_w_flat sum_data non_nan_mask cond1 cond2 cond3 cond4 cond5 percentage_tot


%% Visualisation: relative importance current and waves (1/2)
f2a = figure('Position',[740, 957, 1719, 628]);

% Plot each column with different color and symbol
hold on
for i = 1:size(data_c, 2)
    scatter(data_c(:, i), data_w(:, i), 20, 'Marker',newsymbols{i},...
        'MarkerEdgeColor',newcolors(i,:), 'LineWidth',1);
end

% Plot threshold lines
plot(x_cw, threshold - x_cw, '-k', 'LineWidth',2); % incipient motion
plot(x_cw, 1 - x_cw, '-k', 'LineWidth',2); % sheet flow (Kleinhans, 2002)
a = 1;
b = .01;
c = threshold-b;
line([a-.5, .015325], [a-.5, .015325], 'Color','k', 'LineStyle','--', 'LineWidth',3)
line([a, c], [b, b], 'Color','k', 'LineStyle','--', 'LineWidth',3)
line([b, b], [a, c], 'Color','k', 'LineStyle','--', 'LineWidth',3)

% Add region descriptions
text(4e-4, .6, 'waves only', 'FontSize',fontsize*.6)
text(.014, .6, 'wave-dominated', 'FontSize',fontsize*.6)
text(.6, .25, 'current-dominated', 'FontSize',fontsize*.6, 'Rotation',-90)
text(.6, 2e-3, 'current only', 'FontSize',fontsize*.6, 'Rotation',-90)
text(4e-4, 2e-4, 'no motion', 'FontSize',fontsize*.6)

% Add region fractions
text(8e-4, .2, [mat2str(percentage_cond1, 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(.03, .2, [mat2str(percentage_cond2, 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(.1, .03, [mat2str(percentage_cond3, 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(.1, 1e-3, [mat2str(percentage_cond4, 1), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
text(8e-4, 1e-3, [mat2str(percentage_cond5, 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')

% Additional plotting commands
% plot(x_cw, threshold - x_cw, '--w', 'LineWidth', 2)
% text(threshold - .027, threshold - .019, '\theta_{cr}', 'FontSize', fontsize, 'Color', 'w')
% xline(threshold, '-', '\theta_{cr}', 'LineWidth', 2, 'FontSize', fontsize)
% yline(threshold, '-', '\theta_{cr}', 'LineWidth', 2, 'FontSize', fontsize)
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

clearvars a b c threshold


%% Calculations: relative importance current and waves (3/3)
threshold = GS_stats.theta_cr50_McCarron;

% Initialize matrices to store percentages for each condition
nCols = width(data_c);
percentage_cond1 = zeros(1, nCols);
percentage_cond2 = zeros(1, nCols);
percentage_cond3 = zeros(1, nCols);
percentage_cond4 = zeros(1, nCols);
percentage_cond5 = zeros(1, nCols);

for col = 1:nCols
    % Extract the current column
    data_c_col = data_c(:, col);
    data_w_col = data_w(:, col);
    
    % Calculate the sum of data_c and data_w for the current column
    sum_data = data_c_col + data_w_col;
    
    % Create a mask for non-NaN elements for the current column
    non_nan_mask = ~isnan(data_c_col) & ~isnan(data_w_col);
    
    % Condition 1: waves only
    cond1 = sum_data > threshold(i) & data_c_col < 0.01 & non_nan_mask;
    percentage_cond1(col) = (sum(cond1) / sum(non_nan_mask)) * 100;
    
    % Condition 2: wave-dominated
    cond2 = sum_data > threshold(i) & data_w_col > data_c_col & data_c_col > 0.01 & data_w_col > 0.01 & non_nan_mask;
    percentage_cond2(col) = (sum(cond2) / sum(non_nan_mask)) * 100;
    
    % Condition 3: current-dominated
    cond3 = sum_data > threshold(i) & data_c_col > data_w_col & data_c_col > 0.01 & data_w_col > 0.01 & non_nan_mask;
    percentage_cond3(col) = (sum(cond3) / sum(non_nan_mask)) * 100;
    
    % Condition 4: current only
    cond4 = sum_data > threshold(i) & data_w_col < 0.01 & non_nan_mask;
    percentage_cond4(col) = (sum(cond4) / sum(non_nan_mask)) * 100;
    
    % Condition 5: no motion
    cond5 = sum_data <= threshold(i) & non_nan_mask;
    percentage_cond5(col) = (sum(cond5) / sum(non_nan_mask)) * 100;
end

clearvars data_c_col data_w_col sum_data non_nan_mask cond1 cond2 cond3 cond4 cond5 percentage_tot nCols col


%% Visualisation: relative importance current and waves (2/2)
f2b = figure('Position',[740, 1665, 1719, 628]);
tiledlayout(1,3, 'TileSpacing','tight')

j = 1;
ax = gobjects(1,3);
for i = [1, 4, 6]
    ax(j) = nexttile; hold on
    scatter(data_c(:, i), data_w(:, i), 50, 'Marker','o',...
        'MarkerEdgeColor',newcolors(i,:), 'LineWidth',1);
    text(.47, .7, instruLocs{i}, 'FontSize',fontsize, 'FontWeight','bold', 'Color','r')

    % Plot threshold lines
    plot(x_cw, threshold(i) - x_cw, '-k', 'LineWidth',2); % incipient motion
    plot(x_cw, 1 - x_cw, '-k', 'LineWidth',2); % sheet flow (Kleinhans, 2002)
    a = 1;
    b = .01;
    c = threshold(i)-b;
    line([a-.5, threshold(i)/2], [a-.5, threshold(i)/2], 'Color','k', 'LineStyle','--', 'LineWidth',3)
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
    text(8e-4, .2, [mat2str(percentage_cond1(i), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(.03, .2, [mat2str(percentage_cond2(i), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(.1, .03, [mat2str(percentage_cond3(i), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(.1, 1e-3, [mat2str(percentage_cond4(i), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
    text(8e-4, 1e-3, [mat2str(percentage_cond5(i), 2), '%'], 'FontSize',fontsize*.8, 'FontWeight','bold')
        
    % Set axes properties
    set(gca, 'XScale','log');
    set(gca, 'YScale','log');
    
    % Set axis limits
    xlim([1e-4, 1])
    ylim([1e-4, 1])

    % Label axes
    xlabel('\theta_{c}')

    % Enable grid and set axis properties
    zoom xon
    linkaxes
    grid on
    axis square

    j = j+1;
end

% Set axis properties
ylabel(ax(1), '\theta_{w}')
yticklabels(ax(2:3), {})

clearvars data_c data_w x_cw a b c ax threshold j


%% Only consider time steps without NaNs
rowsWithoutNaNs = any(ismissing(TT_tau_cw), 2);

TT_theta_cwMg_NaN = TT_theta_cwMg;
TT_theta_cw10_NaN = TT_theta_cw10;
TT_theta_cw50_NaN = TT_theta_cw50;
TT_theta_cw90_NaN = TT_theta_cw90;

% Get the size of the timetable
[numRows, numCols] = size(TT_theta_cw50_NaN);

% Loop through each row to check for NaNs
for i = 1:numRows
    % Check if there are any NaNs in the current row
    if ~rowsWithoutNaNs
        % Assign NaNs to the entire row
        TT_theta_cwMg_NaN{i, :} = NaN(1, numCols);
        TT_theta_cw10_NaN{i, :} = NaN(1, numCols);
        TT_theta_cw50_NaN{i, :} = NaN(1, numCols);
        TT_theta_cw90_NaN{i, :} = NaN(1, numCols);
    end
end

clearvars numRows numCols


%% Computations: mobilisation duration (Shields)

% Create an empty table
percentExceed_Soulsby = array2table(NaN(length(instruLocs), length(GS_fractions(2:end))),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions(2:end));
percentExceed_Egiazaroff = array2table(NaN(length(instruLocs), length(GS_fractions(2:end))),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions(2:end));
percentExceed_McCarron = array2table(NaN(length(instruLocs), length(GS_fractions(2:end))),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions(2:end));

% Get the frequency of the time axis
dt = TT_theta_cw50_NaN.Properties.TimeStep;
totalDuration = dateEnd-dateStart;

% Find rows without any NaNs
numRowsWithoutNaNs = sum(rowsWithoutNaNs);
totalDuration_noNaN = numRowsWithoutNaNs * dt;
totalDurationFraction = totalDuration_noNaN / totalDuration;

% Calculate the time percentages of exceedance
disp('McCarron:')
% for i = 1:width(TT_theta_cwMg_NaN)
% 
%     % Find the moments when the threshold is exceeded
%     exceededTimes_Soulsby = sum(TT_theta_cwMg_NaN.(i) > GS_stats.theta_crMg_Soulsby(i));
%     exceededTimes_Egiazaroff = sum(TT_theta_cwMg_NaN.(i) > GS_stats.theta_crMg_Soulsby(i));
%     exceededTimes_McCarron = sum(TT_theta_cwMg_NaN.(i) > GS_stats.theta_crMg_McCarron(i));
% 
%     % Find the duration when the threshold is exceeded
%     exceededDuration_Soulsby = exceededTimes_Soulsby*dt;
%     percent_Mg_Soulsby = exceededDuration_Soulsby/totalDuration_noNaN*100;
%     exceededDuration_Egiazaroff = exceededTimes_Egiazaroff*dt;
%     percent_Mg_Egiazaroff = exceededDuration_Egiazaroff/totalDuration_noNaN*100;
%     exceededDuration_McCarron = exceededTimes_McCarron*dt;
%     percent_Mg_McCarron = exceededDuration_McCarron/totalDuration_noNaN*100;
% 
%     % Append the percentage value to the table
%     percentExceed_Soulsby.Mg(i) = percent_Mg_Soulsby;
%     percentExceed_Egiazaroff.Mg(i) = percent_Mg_Egiazaroff;
%     percentExceed_McCarron.Mg(i) = percent_Mg_McCarron;
% 
%     % Display the total duration
%     disp(['Threshold Mg at ', char(instruLocs(i)), ' exceeded: ', char(exceededDuration_McCarron),...
%         ' h (', num2str(percent_Mg_McCarron, 2),'%) of ', char(totalDuration_noNaN), ' h']);
% 
% end

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

% Calculate mean row and column
meanRow_Soulsby = mean(percentExceed_Soulsby.Variables, 1);
meanCol_Soulsby = mean(percentExceed_Soulsby.Variables, 2);

meanRow_Egiazaroff = mean(percentExceed_Egiazaroff.Variables, 1);
meanCol_Egiazaroff = mean(percentExceed_Egiazaroff.Variables, 2);

meanRow_McCarron = mean(percentExceed_McCarron.Variables, 1);
meanCol_McCarron = mean(percentExceed_McCarron.Variables, 2);

% Add mean row and column to the table
percentExceed_Soulsby{'Mean', :} = meanRow_Soulsby;
percentExceed_Soulsby.Mean = [meanCol_Soulsby; mean(meanRow_Soulsby)];

percentExceed_Egiazaroff{'Mean', :} = meanRow_Egiazaroff;
percentExceed_Egiazaroff.Mean = [meanCol_Egiazaroff; mean(meanRow_Egiazaroff)];

percentExceed_McCarron{'Mean', :} = meanRow_McCarron;
percentExceed_McCarron.Mean = [meanCol_McCarron; mean(meanRow_McCarron)];

clearvars dt exceededDuration_Soulsby exceededDuration_Egiazaroff exceededDuration_McCarron exceededTimes_Soulsby exceededTimes_Egiazaroff exceededTimes_McCarron i percent_d10_Soulsby percent_d10_Egiazaroff percent_d10_McCarron percent_d50_Soulsby percent_d50_Egiazaroff percent_d50_McCarron percent_d90_Soulsby percent_d90_Egiazaroff percent_d90_McCarron meanCol_Soulsby meanCol_Egiazaroff meanCol_McCarron meanRow_Soulsby meanRow_Egiazaroff meanRow_McCarron totalDuration totalDuration_noNaN totalDurationFraction numRowsWithoutNaNs rowsWithoutNaNs


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
f3b = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_tau_cw.(i);

     % Remove NaN values
    data = data(~isnan(data));
    
    % Create the KDE plot
    [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
xlim([0, 4])

% Add title, labels and legend
title('Kernel Density Estimates of Bed Shear Stress')
xlabel('\tau_{cw} (Pa)')
ylabel('Density')
legend('show')

clearvars data xi f i


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
f4b = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_phi_b50.(i);

     % Remove NaN values
    data = data(~isnan(data));
    
    % Create the KDE plot
    [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
xlim([0, 1])

% Add title, labels and legend
title('Kernel Density Estimates of Nondimensional Bedload Transport Rate')
xlabel('\phi_{b}')
ylabel('Density')
legend('show')

clearvars data xi f i


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
f5b = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_b50.(i);

     % Remove NaN values
    data = data(~isnan(data));
    
    % Create the KDE plot
    [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
end
hold off
xlim([0, 1e-4])

% Add title, labels and legend
title('Kernel Density Estimates of Bedload Transport Rate')
xlabel('q_{b} (m^2 s^{-1})')
ylabel('Density')
legend('show')

clearvars data xi f i


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

TT_theta_wMg_NaN = TT_theta_wMg;
TT_theta_w10_NaN = TT_theta_w10;
TT_theta_w50_NaN = TT_theta_w50;
TT_theta_w90_NaN = TT_theta_w90;

TT_theta_cwMg_NaN = TT_theta_cwMg;
TT_theta_cw10_NaN = TT_theta_cw10;
TT_theta_cw50_NaN = TT_theta_cw50;
TT_theta_cw90_NaN = TT_theta_cw90;

TT_phi_bMg_NaN = TT_phi_bMg;
TT_phi_b10_NaN = TT_phi_b10;
TT_phi_b50_NaN = TT_phi_b50;
TT_phi_b90_NaN = TT_phi_b90;

TT_q_bMg_NaN = TT_q_bMg;
TT_q_b10_NaN = TT_q_b10;
TT_q_b50_NaN = TT_q_b50;
TT_q_b90_NaN = TT_q_b90;

TT_q_bMg_zeros = TT_q_bMg;
TT_q_b10_zeros = TT_q_b10;
TT_q_b50_zeros = TT_q_b50;
TT_q_b90_zeros = TT_q_b90;

TT_q_bMg_net_zeros = TT_q_bMg_net;
TT_q_b10_net_zeros = TT_q_b10_net;
TT_q_b50_net_zeros = TT_q_b50_net;
TT_q_b90_net_zeros = TT_q_b90_net;

% Assign NaNs to the entire row
TT_tau_c_NaN{nanRowIndices, :} = NaN;
TT_tau_w_NaN{nanRowIndices, :} = NaN;
TT_tau_cw_NaN{nanRowIndices, :} = NaN;

TT_theta_cMg_NaN{nanRowIndices, :} = NaN;
TT_theta_c10_NaN{nanRowIndices, :} = NaN;
TT_theta_c50_NaN{nanRowIndices, :} = NaN;
TT_theta_c90_NaN{nanRowIndices, :} = NaN;

TT_theta_wMg_NaN{nanRowIndices, :} = NaN;
TT_theta_w10_NaN{nanRowIndices, :} = NaN;
TT_theta_w50_NaN{nanRowIndices, :} = NaN;
TT_theta_w90_NaN{nanRowIndices, :} = NaN;

TT_theta_cwMg_NaN{nanRowIndices, :} = NaN;
TT_theta_cw10_NaN{nanRowIndices, :} = NaN;
TT_theta_cw50_NaN{nanRowIndices, :} = NaN;
TT_theta_cw90_NaN{nanRowIndices, :} = NaN;

TT_phi_bMg_NaN{nanRowIndices, :} = NaN;
TT_phi_b10_NaN{nanRowIndices, :} = NaN;
TT_phi_b50_NaN{nanRowIndices, :} = NaN;
TT_phi_b90_NaN{nanRowIndices, :} = NaN;

TT_q_bMg_NaN{nanRowIndices, :} = NaN;
TT_q_b10_NaN{nanRowIndices, :} = NaN;
TT_q_b50_NaN{nanRowIndices, :} = NaN;
TT_q_b90_NaN{nanRowIndices, :} = NaN;

% Assign zeros to the entire row
TT_q_bMg_zeros{nanRowIndices, :} = 0;
TT_q_b10_zeros{nanRowIndices, :} = 0;
TT_q_b50_zeros{nanRowIndices, :} = 0;
TT_q_b90_zeros{nanRowIndices, :} = 0;

TT_q_bMg_net_zeros{nanRowIndices, :} = 0;
TT_q_b10_net_zeros{nanRowIndices, :} = 0;
TT_q_b50_net_zeros{nanRowIndices, :} = 0;
TT_q_b90_net_zeros{nanRowIndices, :} = 0;

clearvars nanRows nanRowIndices


%% Calculation: Total Bedload Transport Volume (Gross)
Mg = NaN(length(instruLocs), 1);
d10 = NaN(length(instruLocs), 1);
d50 = NaN(length(instruLocs), 1);
d90 = NaN(length(instruLocs), 1);

% Convert datetime to duration and then to seconds
time_seconds = seconds(TT_q_b10_zeros.Time - TT_q_b10_zeros.Time(1));

for i = 1:length(instruLocs)
    
    % Integrate using the trapezoidal rule
    Mg(i) = trapz(time_seconds, TT_q_bMg_zeros.(i));
    d10(i) = trapz(time_seconds, TT_q_b10_zeros.(i));
    d50(i) = trapz(time_seconds, TT_q_b50_zeros.(i));
    d90(i) = trapz(time_seconds, TT_q_b90_zeros.(i));
    
    % Display the result
    % fprintf('Total amount of bedload transport: %.2f m^3/m width\n', totalBedload);
end

totalBedload_gross = table(d10, d50, d90, 'RowNames', flipud(instruLocs));

% Calculate mean row and column
meanRow = mean(totalBedload_gross.Variables, 1);
meanCol = mean(totalBedload_gross.Variables, 2);

% Add mean row and column to the table
totalBedload_gross{'Mean', :} = meanRow;
totalBedload_gross.Mean = [meanCol; mean(meanRow)];

clearvars Mg d10 d50 d90 i meanCol meanRow time_seconds


%% Calculation: Total Bedload Transport Volume (Net)
Mg = NaN(length(instruLocs), 1);
d10 = NaN(length(instruLocs), 1);
d50 = NaN(length(instruLocs), 1);
d90 = NaN(length(instruLocs), 1);

% Convert datetime to duration and then to seconds
time_seconds = seconds(TT_q_b10_net_zeros.Time - TT_q_b10_net_zeros.Time(1));

for i = 1:length(instruLocs)
    
    % Integrate using the trapezoidal rule
    Mg(i) = trapz(time_seconds, TT_q_bMg_net_zeros.(i));
    d10(i) = trapz(time_seconds, TT_q_b10_net_zeros.(i));
    d50(i) = trapz(time_seconds, TT_q_b50_net_zeros.(i));
    d90(i) = trapz(time_seconds, TT_q_b90_net_zeros.(i));
    
    % Display the result
    % fprintf('Total amount of bedload transport: %.2f m^3/m width\n', totalBedload);
end

totalBedload_net = table(d10, d50, d90, 'RowNames', flipud(instruLocs));

% Calculate mean row and column
meanRow = mean(totalBedload_net.Variables, 1);
meanCol = mean(totalBedload_net.Variables, 2);

% Add mean row and column to the table
totalBedload_net{'Mean', :} = meanRow;
totalBedload_net.Mean = [meanCol; mean(meanRow)];

clearvars Mg d10 d50 d90 i meanCol meanRow time_seconds


%% Visualisation: Bed Shear Stress (timeseries)
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


%% Visualisation: Bed Shear Stress (pdf)
f3d = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_tau_cw_NaN.(i);

     % Remove NaN values
    data = data(~isnan(data));
    
    % Create the KDE plot
    [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
end
hold off
xlim([0, 3])

% Add title, labels and legend
title('Kernel Density Estimates of Bed Shear Stress')
xlabel('\tau_{cw} (Pa)')
ylabel('Density')
legend('show')

clearvars data xi f i


%% Visualisation: Nondimensional Bedload Transport Rate (timeseries)
ax = gobjects(1, length(instruLocs));

f4c = figure('Position',[740, 957, 1719, 1336]);
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

clearvars ax tl i


%% Visualisation: Nondimensional Bedload Transport Rate (pdf)
f4d = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)
    
    % Extract the variable data for the current location
    data = TT_phi_b50_NaN.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Create the KDE plot
    [f, xi] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
xlim([0, .2])

% Add title, labels and legend
title('Kernel Density Estimates of Nondimensional Bedload Transport Rate')
xlabel('\phi_{b}')
ylabel('Density')
legend('show')

clearvars data xi f i


%% Visualisation: Bedload Transport Rate (timeseries)
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


%% Visualisation: Bedload Transport Rate (pdf)
f5d = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_b50_zeros.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Create the KDE plot
    [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
xlim([0, 6e-5])

% Add title, labels and legend
title('Kernel Density Estimates of Bedload Transport')
xlabel('q_{b} (m^2 s^{-1})')
ylabel('Density')
legend('show')

clearvars data xi f i


%% Visualisation: Bedload Transport Rate (pdf)
p = gobjects(size(instruLocs));

f5e = figure('Position',[740, 957, 1719, 1336]);

hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_b90_net_zeros.(i);

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

    % Find the location of the peak frequency (mode) and 
    [max_f, max_index] = max(smoothed_f);
    peak_location = xi(max_index);
    meanValue = mean(data);

    % Create the KDE plot
    plot(xi, smoothed_f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
    xline(peak_location, 'LineStyle','-', 'LineWidth',2, 'Color',newcolors(i,:), 'HandleVisibility','off')
    xline(meanValue, 'LineStyle','--', 'LineWidth',2, 'Color',newcolors(i,:), 'HandleVisibility','off')

end
hold off
xlim([-.6e-4, .6e-4])

% Add title, labels and legend
title('Kernel Density Estimates of Bedload Transport')
xlabel('< NE        q_{b,50} (m^2 s^{-1})        SW >')
ylabel('Density')
legend('show')

clearvars data xi f i max_f max_index peak_location smoothed_f gaussianWindowSize numPoints meanValue p


%% Visualisation: Bed Shear Stress / Gross/Net BL Transport Rate (pdf)
f6a = figure('Position',[740, 1665, 1719, 628]);
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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

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
    data = TT_q_b10_zeros.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Create the KDE plot
    [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
axis tight
xlim([0, 5e-5])

% Add title, labels and legend
xlabel('q_{b,10,gross} (m^2 s^{-1})')
legend('show', 'Location','east')
yticklabels({})

nexttile; hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_b10_net_zeros.(i);

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
    plot(xi, smoothed_f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
    xline(peak_location, 'LineStyle','-', 'LineWidth',3, 'Color',newcolors(i,:), 'HandleVisibility','off')
    xline(meanValue, 'LineStyle','--', 'LineWidth',3, 'Color',newcolors(i,:), 'HandleVisibility','off')

end
hold off
axis tight
xlim([-.4e-4, .4e-4])

% Add title, labels and legend
xlabel('q_{b,10,net} (m^2 s^{-1})')
yticklabels({})


% Add figure annotations
annotation('textbox', [0.31, 0.8, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.605, 0.8, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.895, 0.8, 0.1, 0.1], 'String','(c)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.66, 0.6, 0.1, 0.1], 'String','SW direction', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.83, 0.6, 0.1, 0.1], 'String','NE direction', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');

clearvars data xi f i numPoints bw gaussianWindowSize smoothed_f max_f max_index peak_location meanValue xLim yLim


%% Visualisation: Bed Shear Stress / Gross/Net BL Transport Rate (pdf)
f6b = figure('Position',[740, 1665, 1719, 628]);
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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

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
    data = TT_q_b50_zeros.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Create the KDE plot
    [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
axis tight
xlim([0, 5e-5])

% Add title, labels and legend
xlabel('q_{b,50,gross} (m^2 s^{-1})')
legend('show', 'Location','east')
yticklabels({})

nexttile; hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_b50_net_zeros.(i);

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
    plot(xi, smoothed_f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
    xline(peak_location, 'LineStyle','-', 'LineWidth',3, 'Color',newcolors(i,:), 'HandleVisibility','off')
    xline(meanValue, 'LineStyle','--', 'LineWidth',3, 'Color',newcolors(i,:), 'HandleVisibility','off')

end
hold off
axis tight
xlim([-.4e-4, .4e-4])

% Add title, labels and legend
xlabel('q_{b,50,net} (m^2 s^{-1})')
yticklabels({})

% Add figure annotations
annotation('textbox', [0.31, 0.8, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.605, 0.8, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.895, 0.8, 0.1, 0.1], 'String','(c)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.66, 0.6, 0.1, 0.1], 'String','SW direction', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.83, 0.6, 0.1, 0.1], 'String','NE direction', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');

clearvars data xi f i numPoints bw gaussianWindowSize smoothed_f max_f max_index peak_location meanValue


%% Visualisation: Bed Shear Stress / Gross/Net BL Transport Rate (pdf)
f6c = figure('Position',[740, 1665, 1719, 628]);
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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

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
    data = TT_q_b90_zeros.(i);

     % Remove NaN values
    data = data(~isnan(data));

    % Create the KDE plot
    [f, xi, bw] = ksdensity(data, 'BoundaryCorrection','reflection',...
        'Support','nonnegative');
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))

end
hold off
axis tight
xlim([0, 5e-5])

% Add title, labels and legend
xlabel('q_{b,90,gross} (m^2 s^{-1})')
legend('show', 'Location','east')
yticklabels({})

nexttile; hold on
% Loop through each location to create KDE plots
for i = 1:length(instruLocs)

    % Extract the variable data for the current location
    data = TT_q_b90_net_zeros.(i);

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
    plot(xi, smoothed_f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',newcolors(i,:))
    xline(peak_location, 'LineStyle','-', 'LineWidth',3, 'Color',newcolors(i,:), 'HandleVisibility','off')
    xline(meanValue, 'LineStyle','--', 'LineWidth',3, 'Color',newcolors(i,:), 'HandleVisibility','off')

end
hold off
axis tight
xlim([-.4e-4, .4e-4])

% Add title, labels and legend
xlabel('q_{b,90,net} (m^2 s^{-1})')
yticklabels({})

% Add figure annotations
annotation('textbox', [0.31, 0.8, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.605, 0.8, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.895, 0.8, 0.1, 0.1], 'String','(c)', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.66, 0.6, 0.1, 0.1], 'String','SW direction', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');
annotation('textbox', [0.83, 0.6, 0.1, 0.1], 'String','NE direction', 'EdgeColor','none', 'FontSize',fontsize, 'FontWeight','bold');

clearvars data xi f i numPoints bw gaussianWindowSize smoothed_f max_f max_index peak_location meanValue


%% Visualisation: Bed Shear Stress / Gross/Net BL Transport Rate (pdf)
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
    data = TT_q_bMg_net_zeros.(i);

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
    xline(peak_location, 'LineStyle',':', 'LineWidth',2, 'Color','k', 'HandleVisibility','off')

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
