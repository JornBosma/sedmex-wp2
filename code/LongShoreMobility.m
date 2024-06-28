%% Initialization
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;

dataPath = fullfile(filesep, 'Volumes', 'T7 Shield', 'DataDescriptor', 'hydrodynamics', 'ADV', filesep);

instruments = ["L1C1VEC", "L2C4VEC", "L3C1VEC", "L4C1VEC", "L5C1VEC", "L6C1VEC"];
sampleLocs = {'L1', 'L2', 'L3.5', 'L4', 'Tmb', 'L6'};
instruLocs = {'L1', 'L2', 'L3', 'L4', 'L5', 'L6'};
GS_fractions = {'d10', 'd50', 'd90'};

newcolors = [cbf.vermilion; cbf.blue; cbf.bluegreen; cbf.yellow; cbf.redpurp; cbf.skyblue; cbf.orange];
newsymbols = {'o', 's', '^', 'd', 'v', 'p'};

g = 9.81;     
rho_s = 2650;  
rho_w = 1025; 


%% Computations: Critical Bed-Shear Stress & Critical Shields

% Sample data L1, L2...L6
fg = [125e-6; 175e-6; 123e-6; 107e-6; 66e-6; 157e-6];
d10 = [277e-6; 292e-6; 315e-6; 277e-6; 296e-6; 244e-6];
d50 = [672e-6; 747e-6; 776e-6; 522e-6; 576e-6; 629e-6];
d90 = [2427e-6; 2844e-6; 2371e-6; 2138e-6; 1746e-6; 3102e-6];

GS_stats = table(d10, d50, d90, fg, 'RowNames', flipud(sampleLocs));

% Initialize arrays for critical Shields and bed-shear stress computations
crit_methods = {'Soulsby', 'Egiazaroff', 'McCarron'};
num_methods = length(crit_methods);
num_percentiles = 3;
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
percentiles = {'10', '50', '90'};

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

clearvars d10 d50 d90 fg theta_cr tau_cr crit_methods prefixes percentiles i j p perc num_methods num_percentiles method


%% Computations: bed-shear stress & Shields
% GS_stats(2,:) = [];  % Exclude L2

% Defined starting time
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');

% Create a datetime array for the time axis
dateStart = datetime(2021, 9, 10);
dateEnd = datetime(2021, 10, 19);
dt = minutes(10);
timeAxis = dateStart:dt:dateEnd;

% Preallocate the data matrix with NaNs
n = length(timeAxis);
emptyData = nan(n, length(instruments));

% List of variable suffixes
suffixes = {'tau_c', 'tau_w', 'tau_cw', 'theta_c10', 'theta_w10', 'theta_cw10', ...
            'phi_b10', 'q_b10', 'theta_c50', 'theta_w50', 'theta_cw50', ...
            'phi_b50', 'q_b50', 'theta_c90', 'theta_w90', 'theta_cw90', ...
            'phi_b90', 'q_b90'};

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

    if i == 1 % L1C1VEC
        % Refine measuring heights L1C1VEC
        % index_measurement = [1, 258, 478, 774, 927, 1079, 1233, 1382, 1968, 2411, 2932, 3236, 3835, 3992, 4575, 4867, 5099, 5543]';
        % z_measured = [.15, .15, .14, .15, .16, .17, .17, .20, .19, .19, .20, .11, .06, .06, .13, .17, .14, .16]';
        index_measurement = [1, 258, 478, 774, 927, 1079, 1233, 1382, 1968, 2411, 2932, 3236, 4575, 4867, 5099, 5543]';
        z_measured = [.15, .15, .14, .15, .16, .17, .17, .20, .19, .19, .20, .11, .13, .17, .14, .16]';
        z_new = nan(max(index_measurement), 1);
        z_new(index_measurement) = z_measured;
        nonNaNIndices = find(~isnan(z_new));
        nonNaNValues = z_new(nonNaNIndices);
        interpolatedValues = interp1(nonNaNIndices, nonNaNValues, 1:length(z_new), 'pchip');
        z_new(isnan(z_new)) = interpolatedValues(isnan(z_new));
        z_old = z;
        z = z_new(1:length(z));
    end

    % Estimate the depth-averaged current velocity
    k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)
    [u_c, ~] = compute_DAV(u_z, z, k_sc, h);

    % % Check new measuring heights
    % figure
    % yyaxis left
    % plot(z_new,'b'); hold on
    % plot(index_measurement, z_measured, 'ob')
    % plot(z_old,'k')
    % yyaxis right
    % plot(u_z, 'r'); hold on
    % plot(u_c, 'y')
    % hold off

    % Compute the shear stress components
    [tau_c, tau_w, tau_cw] = compute_BSS_orbital(u_c, h, Urms, T, rho_w, phi_c, phi_w, GS_stats.d90(i), g);
    % [tau_c, tau_w, tau_cw] = compute_BSS_linearTheory(u_c, H, k, T, h, rho_w, phi_c, phi_w, GS_stats.d90(i), g);

    % % Define two arrays with NaN values
    % array1 = tau_w3;
    % array2 = tau_w4;
    % 
    % % Find the indices where the arrays are not equal and the differences
    % [unequal_indices, differences] = find_unequal_indices_and_differences(array1, array2);
    % 
    % % Display the results
    % disp('Indices where the arrays are not equal:');
    % disp(unequal_indices);
    % 
    % disp('Differences at the unequal indices:');
    % disp(differences);
    % 
    % % Visualise the results
    % scatter(array1, array2)

    % Append to timetable
    firstRow = find(TT_tau_c.Time == tStart);
    lastRow = find(TT_tau_c.Time == tEnd);

    TT_tau_c.(i)(firstRow:lastRow) = tau_c;
    TT_tau_w.(i)(firstRow:lastRow) = tau_w;
    TT_tau_cw.(i)(firstRow:lastRow) = tau_cw;
    
    % Fraction-specific Shields numbers
    [theta_c10, theta_w10, theta_cw10] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d10(i), g);
    [theta_c50, theta_w50, theta_cw50] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d50(i), g);
    [theta_c90, theta_w90, theta_cw90] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, GS_stats.d90(i), g);

    TT_theta_c10.(i)(firstRow:lastRow) = theta_c10;
    TT_theta_w10.(i)(firstRow:lastRow) = theta_w10;
    TT_theta_cw10.(i)(firstRow:lastRow) = theta_cw10;

    TT_theta_c50.(i)(firstRow:lastRow) = theta_c50;
    TT_theta_w50.(i)(firstRow:lastRow) = theta_w50;
    TT_theta_cw50.(i)(firstRow:lastRow) = theta_cw50;
    
    TT_theta_c90.(i)(firstRow:lastRow) = theta_c90;
    TT_theta_w90.(i)(firstRow:lastRow) = theta_w90;
    TT_theta_cw90.(i)(firstRow:lastRow) = theta_cw90;

    % Nondimensional bedload predictors
    alpha = 11;  % calibration coefficient of Ribberink (1998)
    beta = 1.65; % calibration exponent

    phi_b10 = compute_Einstein_parameter(TT_theta_cw10.(i), GS_stats.theta_cr10_McCarron(i), alpha, beta);
    phi_b50 = compute_Einstein_parameter(TT_theta_cw50.(i), GS_stats.theta_cr50_McCarron(i), alpha, beta);
    phi_b90 = compute_Einstein_parameter(TT_theta_cw90.(i), GS_stats.theta_cr90_McCarron(i), alpha, beta);

    TT_phi_b10.(i) = phi_b10;
    TT_phi_b50.(i) = phi_b50;
    TT_phi_b90.(i) = phi_b90;

    % Dimensional bedload transport rate
    q_b10 = compute_transport_rate(phi_b10, rho_w, rho_s, GS_stats.d10(i), g);
    q_b50 = compute_transport_rate(phi_b50, rho_w, rho_s, GS_stats.d50(i), g);
    q_b90 = compute_transport_rate(phi_b90, rho_w, rho_s, GS_stats.d90(i), g);

    TT_q_b10.(i) = q_b10;
    TT_q_b50.(i) = q_b50;
    TT_q_b90.(i) = q_b90;

end

clearvars dataPath dt emptyData filename firstRow h H i info k lastRow n phi_c phi_w t T tau_c tau_w tau_cw theta_c10 theta_c50 theta_c90 theta_w10 theta_w50 theta_w90 theta_cw10 theta_cw50 theta_cw90 timeAxis rho_w u_z u_c h Urms emptyData firstRow lastRow z z_new z_measured nonNaNValues nonNaNIndices interpolatedValues phi_b10 phi_b50 phi_b90 q_b10 q_b50 q_b90 z_old


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
xlim([tStart, tEnd])
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
xlim([tStart, tEnd])
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
xlim([tStart, tEnd])
ylim([0, 3])
zoom xon
linkaxes
grid on

clearvars i


%% Calculations: relative importance current and waves (1/2)

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


%% Calculations: relative importance current and waves (2/2)
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


%% Visualisation: relative importance current and waves
f2 = figure('Position',[740, 957, 1719, 628]);

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

clearvars data_c data_w x_cw threshold_50 a b c


%% Only consider time steps without NaNs
rowsWithoutNaNs = any(ismissing(TT_tau_cw), 2);

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
        TT_theta_cw10_NaN{i, :} = NaN(1, numCols);
        TT_theta_cw50_NaN{i, :} = NaN(1, numCols);
        TT_theta_cw90_NaN{i, :} = NaN(1, numCols);
    end
end

clearvars numRows numCols


%% Computations: mobilisation duration (Shields)

% Create an empty table
percentExceed_Soulsby = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);
percentExceed_Egiazaroff = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);
percentExceed_McCarron = array2table(NaN(length(instruLocs), length(GS_fractions)),...
    'RowNames',instruLocs, 'VariableNames',GS_fractions);

% Get the frequency of the time axis
dt = TT_theta_cw50_NaN.Properties.TimeStep;
totalDuration = dateEnd-dateStart;

% Find rows without any NaNs
numRowsWithoutNaNs = sum(rowsWithoutNaNs);
totalDuration_noNaN = numRowsWithoutNaNs * dt;
totalDurationFraction = totalDuration_noNaN / totalDuration;

% Calculate the time percentages of exceedance
disp('McCarron:')
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

clearvars dt exceededDuration_Soulsby exceededDuration_Egiazaroff exceededDuration_McCarron exceededTimes_Soulsby exceededTimes_Egiazaroff exceededTimes_McCarron i percent_d10_Soulsby percent_d10_Egiazaroff percent_d10_McCarron percent_d50_Soulsby percent_d50_Egiazaroff percent_d50_McCarron percent_d90_Soulsby percent_d90_Egiazaroff percent_d90_McCarron meanCol_Soulsby meanCol_Egiazaroff meanCol_McCarron meanRow_Soulsby meanRow_Egiazaroff meanRow_McCarron totalDuration totalDuration_noNaN totalDurationFraction


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

    text(tEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
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

    text(tEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
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

    text(tEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
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

TT_theta_c10_NaN = TT_theta_c10;
TT_theta_c50_NaN = TT_theta_c50;
TT_theta_c90_NaN = TT_theta_c90;

TT_theta_w10_NaN = TT_theta_w10;
TT_theta_w50_NaN = TT_theta_w50;
TT_theta_w90_NaN = TT_theta_w90;

TT_theta_cw10_NaN = TT_theta_cw10;
TT_theta_cw50_NaN = TT_theta_cw50;
TT_theta_cw90_NaN = TT_theta_cw90;

TT_phi_b10_NaN = TT_phi_b10;
TT_phi_b50_NaN = TT_phi_b50;
TT_phi_b90_NaN = TT_phi_b90;

TT_q_b10_NaN = TT_q_b10;
TT_q_b50_NaN = TT_q_b50;
TT_q_b90_NaN = TT_q_b90;

TT_q_b10_zeros = TT_q_b10;
TT_q_b50_zeros = TT_q_b50;
TT_q_b90_zeros = TT_q_b90;

% Assign NaNs to the entire row
TT_tau_c_NaN{nanRowIndices, :} = NaN;
TT_tau_w_NaN{nanRowIndices, :} = NaN;
TT_tau_cw_NaN{nanRowIndices, :} = NaN;

TT_theta_c10_NaN{nanRowIndices, :} = NaN;
TT_theta_c50_NaN{nanRowIndices, :} = NaN;
TT_theta_c90_NaN{nanRowIndices, :} = NaN;

TT_theta_w10_NaN{nanRowIndices, :} = NaN;
TT_theta_w50_NaN{nanRowIndices, :} = NaN;
TT_theta_w90_NaN{nanRowIndices, :} = NaN;

TT_theta_cw10_NaN{nanRowIndices, :} = NaN;
TT_theta_cw50_NaN{nanRowIndices, :} = NaN;
TT_theta_cw90_NaN{nanRowIndices, :} = NaN;

TT_phi_b10_NaN{nanRowIndices, :} = NaN;
TT_phi_b50_NaN{nanRowIndices, :} = NaN;
TT_phi_b90_NaN{nanRowIndices, :} = NaN;

TT_q_b10_NaN{nanRowIndices, :} = NaN;
TT_q_b50_NaN{nanRowIndices, :} = NaN;
TT_q_b90_NaN{nanRowIndices, :} = NaN;

% Assign zeros to the entire row
TT_q_b10_zeros{nanRowIndices, :} = 0;
TT_q_b50_zeros{nanRowIndices, :} = 0;
TT_q_b90_zeros{nanRowIndices, :} = 0;

clearvars nanRows nanRowIndices


%% Calculation: Total Bedload Transport Volume
d10 = NaN(width(TT_q_b10), 1);
d50 = NaN(length(instruLocs), 1);
d90 = NaN(width(TT_q_b90), 1);

% Convert datetime to duration and then to seconds
time_seconds = seconds(TT_q_b10_zeros.Time - TT_q_b10_zeros.Time(1));

for i = 1:length(instruLocs)
    
    % Integrate using the trapezoidal rule
    d10(i) = trapz(time_seconds, TT_q_b10_zeros.(i));
    d50(i) = trapz(time_seconds, TT_q_b50_zeros.(i));
    d90(i) = trapz(time_seconds, TT_q_b90_zeros.(i));
    
    % Display the result
    % fprintf('Total amount of bedload transport: %.2f m^3/m width\n', totalBedload);
end

totalBedload = table(d10, d50, d90, 'RowNames', flipud(instruLocs));

% Calculate mean row and column
meanRow = mean(totalBedload.Variables, 1);
meanCol = mean(totalBedload.Variables, 2);

% Add mean row and column to the table
totalBedload{'Mean', :} = meanRow;
totalBedload.Mean = [meanCol; mean(meanRow)];

clearvars d10 d50 d90 i meanCol meanRow


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

    text(tEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
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

    text(tEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
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

    text(tEnd-days(1.4), ax(i).YLim(2)*.8, instruLocs(i), 'FontSize',fontsize)

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
end
hold off
xlim([0, 6e-5])

% Add title, labels and legend
title('Kernel Density Estimates of Bedload Transport')
xlabel('q_{b} (m^2 s^{-1})')
ylabel('Density')
legend('show')

clearvars data xi f i


%% Visualisation: Bed Shear Stress (pdf)
f6 = figure('Position',[740, 1665, 1719, 628]);
tiledlayout(1,2, "TileSpacing","compact")

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
end
hold off
xlim([0, 3])

% Add title, labels and legend
xlabel('\tau_{cw} (Pa)')
ylabel('Density')

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
    plot(xi, f, 'LineWidth',3, 'DisplayName',instruLocs{i}, 'Color',cbf.six(i,:))
end
hold off
xlim([0, 6e-5])

% Add title, labels and legend
xlabel('q_{b,50} (m^2 s^{-1})')
legend('show')

clearvars data xi f i