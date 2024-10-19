%% Initialisation
% close all
% clear
% clc

[~, fontsize, ~, ~, ~] = sedmex_init;

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];
instru = {'L2C4VEC', 'L2C5SONTEK1', 'L2C6OSSI'};
% instru = {'L2C4VEC', 'L2C5SONTEK1', 'L2C7ADCP'};
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
for i = 1:length(instru)
    [~, idx] = ismember(common_time, time{i});
    filtered_par{i} = par{i}(idx);
end

% Correct Instrument 2 data (L2C5SONTEK1) based on interpolation
corrected_par2 = filtered_par{2};  % Initialize corrected data for Instrument 2

for j = 1:length(common_time)
    % Get values for the current time step
    val1 = filtered_par{1}(j);  % Instrument 1
    val2 = filtered_par{2}(j);  % Instrument 2
    val3 = filtered_par{3}(j);  % Instrument 3
    
    % Handle the logic based on availability of values
    if isnan(val1) && isnan(val3)
        % Both Instrument 1 and 3 have NaN, keep Instrument 2's original value
        corrected_par2(j) = val2;
        % corrected_par2(j) = NaN;
    elseif isnan(val1)
        % Only Instrument 1 is NaN, use Instrument 3's value
        corrected_par2(j) = val3;
    elseif isnan(val3)
        % Only Instrument 3 is NaN, use Instrument 1's value
        corrected_par2(j) = val1;
    else
        % Neither is NaN, interpolate between Instrument 1 and 3
        corrected_par2(j) = (val1 + val3) / 2;
    end
end

%% Plot corrected data
f1 = figure('Position',[-973, 1668, 3424, 617]);
plot(common_time, filtered_par{2}, 'k-', 'LineWidth', 1.5);  % Original Instrument 2
hold on;
plot(common_time, corrected_par2, 'r-', 'LineWidth', 2);   % Corrected Instrument 2
legend('Original Instrument 2', 'Corrected Instrument 2');
xlabel('Time');
ylabel('Hm0 (Wave Height)');
grid on;
hold off;


%% Compare instruments
f2 = figure('Position',[-973, 965, 3424, 617]);
plot(common_time, corrected_par2, 'r-', 'LineWidth', 1.5);  % Original Instrument 2
hold on;
plot(common_time, filtered_par{1}, 'b-', 'LineWidth', 2);   % Original Instrument 1
plot(common_time, filtered_par{3}, 'y-', 'LineWidth', 2);   % Original Instrument 3
legend('Corrected Instrument 2', 'Original Instrument 1', 'Original Instrument 2');
xlabel('Time');
ylabel('Hm0 (Wave Height)');
grid on;
hold off;


%% Calculate and display the mean for each valid dataset
for i = 1:length(instru)
    fprintf('Mean of %s (original): %.2f\n', instru{i}, mean(filtered_par{i}, 'omitnan'))
end
fprintf('Mean of %s (corrected): %.2f\n', instru{2}, mean(corrected_par2, 'omitnan'));