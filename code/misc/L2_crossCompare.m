%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];
% instru = {'L2C4VEC', 'L2C5SONTEK1', 'L2C5SONTEK2', 'L2C5SONTEK3', 'L2C7ADCP', 'L2C10VEC'};
% instru = {'L2C4VEC', 'L2C5SONTEK1', 'L2C6OSSI', 'L2C7ADCP'};
instru = {'L2C4VEC', 'L2C5SONTEK1', 'L2C6OSSI'};
instruPath = cell(1,length(instru));

for i = 1:length(instru)
    if strcmp(instru{i}, 'L2C7ADCP')
        instruPath{i} = [dataPath 'ADCP' filesep 'tailored_' instru{i} '.nc'];
    elseif strcmp(instru{i}, 'L2C6OSSI')
        instruPath{i} = [dataPath 'pressuresensors' filesep instru{i} filesep 'tailored_' instru{i} '.nc'];
    else
        instruPath{i} = [dataPath 'ADV' filesep instru{i} filesep 'tailored_' instru{i} '.nc'];
    end
end

t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');

% Preallocate cell arrays for time and data
t = cell(1, length(instru));
time = cell(1, length(instru));
par = cell(1, length(instru));

% Read data in a loop to avoid redundancy
for i = 1:length(instru)
    t{i} = ncread(instruPath{i}, 't');
    time{i} = t0 + seconds(t{i});
    par{i} = ncread(instruPath{i}, 'Hm0');
end

% Find the common time range across all instruments
common_time = time{1};
for i = 2:length(instru)
    common_time = intersect(common_time, time{i});
end

% Filter data to keep only common times
filtered_par = cell(1, length(instru));
valid_mask = true(length(common_time), 1);  % A logical mask to keep track of valid data

for i = 1:length(instru)
    [~, idx] = ismember(common_time, time{i});
    filtered_par{i} = par{i}(idx);
    
    % Update valid_mask to remove time steps where any instrument has NaN
    % valid_mask = valid_mask & ~isnan(filtered_par{i});
end

% Apply the valid_mask to filter out invalid (NaN-containing) time steps
valid_time = common_time(valid_mask);
for i = 1:length(instru)
    filtered_par{i} = filtered_par{i}(valid_mask);
end

% Apply simple global correction factor
% par{2} = par{2}*2.6;


%% Plotting (all data)
f1 = figure('Position',[743, 1669, 1708, 617]);
hold on

% Plot each dataset with different colors
for i = 1:length(instru)
    plot(time{i}, par{i}, 'LineWidth', 2)
end

xlim([datetime('2021-09-27'), datetime('2021-10-03')])

% Improve plot with labels and legend
ylabel('Parameter', 'FontSize', fontsize)
legend(instru, 'FontSize', fontsize)
grid on
hold off


%% Plotting (no NaNs)
f2 = figure('Position',[743, 964, 1708, 617]);
hold on

% Plot each dataset with different colors
for i = 1:length(instru)
    plot(valid_time, filtered_par{i}, 'o', 'LineWidth', 2)
end

xlim([datetime('2021-09-12'), datetime('2021-10-03')])

% Improve plot with labels and legend
ylabel('Parameter', 'FontSize', fontsize)
legend(instru, 'FontSize', fontsize)
grid on
hold off


%% Calculate and display the mean for each valid dataset
for i = 1:length(instru)
    fprintf('Mean of %s: %.2f\n', instru{i}, mean(filtered_par{i}, 'omitnan'))
end


%% Correct SONTEK by just interpolation
% Initialize the corrected data for the second time series
par2_corrected = filtered_par{2};  % Start by copying original ts2 data

% Iterate through each time step to apply the correction
for i = 1:length(filtered_par{2})
    if isnan(filtered_par{1}(i)) || isnan(filtered_par{3}(i))
        % If either ts1 or ts3 has a NaN, keep the original ts2 value
        par2_corrected(i) = filtered_par{2}(i);
    else
        % If both ts1 and ts3 have values, apply interpolation
        par2_corrected(i) = (filtered_par{1}(i) + filtered_par{3}(i)) / 2;
    end
end

% Replace the second time series in the cell array with the corrected one
filtered_par{2} = par2_corrected;

fprintf('Mean of corrected %s: %.2f\n', instru{2}, mean(filtered_par{2}, 'omitnan'))


%% Plotting (no NaNs, after correction)
f3 = figure('Position',[743, 1669, 1708, 617]);
hold on

% Plot each dataset with different colors
for i = 1:length(instru)
    plot(valid_time, filtered_par{i}, 'o', 'LineWidth', 2)
end

xlim([datetime('2021-09-12'), datetime('2021-10-03')])

% Improve plot with labels and legend
ylabel('Parameter', 'FontSize', fontsize)
legend(instru, 'FontSize', fontsize)
grid on
hold off


%% Correct SONTEK by just interpolation and global correction factor
% % Initialize the corrected data for the second time series
%     par2_corrected = filtered_par{2};  % Start by copying original ts2 data
% 
%     % Correction factor arrays
%     correction_factors = nan(size(filtered_par{2}));  % To store correction factor at each time step
% 
%     % Iterate through each time step to apply the correction
%     for i = 1:length(filtered_par{2})
%         if isnan(filtered_par{1}(i)) || isnan(filtered_par{3}(i))
%             % If either ts1 or ts3 has a NaN, we will apply a correction later
%             continue;
%         else
%             % If both ts1 and ts3 have values, apply interpolation
%             interp_value = (filtered_par{1}(i) + filtered_par{3}(i)) / 2;
% 
%             % Calculate the correction factor (ratio) based on the original and corrected value
%             correction_factors(i) = interp_value / filtered_par{2}(i);
% 
%             % Apply correction for this valid time step
%             par2_corrected(i) = interp_value;
%         end
%     end
% 
%     % Compute a global correction factor: mean of valid correction factors
%     valid_factors = correction_factors(~isnan(correction_factors));
%     global_correction_factor = mean(valid_factors);
% 
%     % Apply the global correction to NaN segments
%     for i = 1:length(filtered_par{2})
%         if isnan(filtered_par{1}(i)) || isnan(filtered_par{3}(i))
%             % Apply the global correction to the original value
%             par2_corrected(i) = filtered_par{2}(i) * global_correction_factor;
%         end
%     end
% 
%     % Replace the second time series in the cell array with the corrected one
%     filtered_par{2} = par2_corrected;


%% Correct SONTEK by just interpolation and local correction factor
% % Initialize the corrected data for the second time series
% par2_corrected = filtered_par{2};  % Start by copying original ts2 data
% 
% % Correction factor arrays
% correction_factors = nan(size(filtered_par{2}));  % To store correction factor at each time step
% 
% % Iterate through each time step to calculate the correction factors
% for i = 1:length(filtered_par{2})
%     if isnan(filtered_par{1}(i)) || isnan(filtered_par{3}(i))
%         % Leave correction factor as NaN for now (will be interpolated later)
%         continue;
%     else
%         % If both ts1 and ts3 have values, calculate the interpolation
%         interp_value = (filtered_par{1}(i) + filtered_par{3}(i)) / 2;
% 
%         % Calculate the correction factor (ratio) based on the original and corrected value
%         correction_factors(i) = interp_value / filtered_par{2}(i);
% 
%         % Apply correction for this valid time step
%         par2_corrected(i) = interp_value;
%     end
% end
% 
% % Interpolate NaN correction factors using the nearest valid ones
% valid_indices = find(~isnan(correction_factors));
% correction_factors_interpolated = interp1(valid_indices, correction_factors(valid_indices), 1:length(filtered_par{2}), 'linear', 'extrap');
% 
% % Apply interpolated correction factors to the NaN segments
% for i = 1:length(filtered_par{2})
%     if isnan(filtered_par{1}(i)) || isnan(filtered_par{3}(i))
%         % Apply the interpolated correction factor to the original value
%         par2_corrected(i) = filtered_par{2}(i) * correction_factors_interpolated(i);
%     end
% end
% 
% % Replace the second time series in the cell array with the corrected one
% filtered_par{2} = par2_corrected;
