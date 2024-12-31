%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;
% fontsize = 26; % ultra-wide screen

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];

% Coastline orientations
coastAngleNdeg = 50; % coastline angle northern beach [deg North]
coastAngleSdeg = 40; % coastline angle southern beach [deg North]

coastAngleN = deg2rad(90-coastAngleNdeg); % coastline angle northern beach [rad wrt East]
coastAngleS = deg2rad(90-coastAngleSdeg); % coastline angle southern beach [rad wrt East]

% angleDeg = -(rad2deg(angleRad)-90);


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


%% Load ADV data
instru = {'L6C1VEC', 'L5C1VEC', 'L4C1VEC', 'L3C1VEC', 'L2C5SONTEK1', 'L1C1VEC'};

t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)

TT = cell(1, length(instru));
for n = 1:length(instru)

    if strcmp(instru{n}, 'L2C7ADCP')
        ADVpath = [dataPath 'ADCP' filesep 'tailored_' instru{n} '.nc'];
    else
        ADVpath = [dataPath 'ADV' filesep instru{n} filesep 'tailored_' instru{n} '.nc'];
    end

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

    % Correct SONTEK wave height based on values L2C4ADV and L2C6OSSI
    if strcmp(instru{n}, 'L2C5SONTEK1')
        [~, Hm0] = correct_sontek(t, Hm0, 'Hm0');
    end

    % Interpolate the measurements to the new time vector
    time = t0 + seconds(t);
    habCVnew = retime(hab, time, 'pchip');
    z = habCVnew.(instru{n});

    % Estimate the depth-averaged current velocity
    [Ud, ~] = compute_DAV(abs(Uz), z, k_sc, h);
    % [Udl, ~] = compute_DAV(abs(Ulong), z, k_sc, h);
    % [Udc, ~] = compute_DAV(abs(Ucros), z, k_sc, h);

    % % Convert Cartesian direction convention into Nautical direction
    % Uang = wrapTo360(90-Uang);
    % Hdir = wrapTo360(90-Hdir);

    % % Checks if measurement volume is not buried
    % buried1 = hab < .05;
    % buried2 = U < .002;
    % 
    % t(buried1 | buried2) = NaN;
    % Ud(buried1 | buried2) = NaN;
    % Uang(buried1 | buried2) = NaN;
    % Hm0(buried1 | buried2) = NaN;
    % Hdir(buried1 | buried2) = NaN;

    TT{n} = timetable(time, z, Uz, Ud, Udir, Hm0, Hdir, Ulong, Ucros);
end

clearvars t0 t info h Uz Udir Hm0 Hdir time habCVnew z Ud n hab_measured dataPath ADVpath Ulong Ucros


%% Calculations
crossLongRatio = nan(size(TT));
for i = 1:length(TT)
    crossLongRatio(i) = mean(abs(TT{i}.Ucros), 'omitmissing') / mean(abs(TT{i}.Ulong), 'omitmissing') * 100;
end

mean(crossLongRatio)


%% Isolate the variables

% Step 1: Determine the common time vector
% Here we use the union of all time points from the timetables
commonTime = TT{1}.time;
for i = 2:numel(TT)
    commonTime = union(commonTime, TT{i}.time);
end

% Step 2: Synchronize all timetables to the common time vector
for i = 1:numel(TT)
    TT{i} = retime(TT{i}, commonTime, 'fillwithmissing');
end

% Step 3: Extract variable names (assuming all timetables have the same variables)
variableNames = TT{1}.Properties.VariableNames;

% Step 4: Initialize arrays to hold 1x6 matrices for each variable
numVars = numel(variableNames);
numTimes = numel(commonTime);
collectedVariables = cell(1, numVars);

for v = 1:numVars
    % Initialize a matrix for each variable
    varData = NaN(numTimes, numel(TT)); % Each column corresponds to a timetable
    
    for i = 1:numel(TT)
        varData(:, i) = TT{i}.(variableNames{v});
    end
    
    collectedVariables{v} = varData;
end

t = TT{1}.time;
z = collectedVariables{1};
Uz = collectedVariables{2};
Ud = collectedVariables{3};
Udir = collectedVariables{4};
Hm0 = collectedVariables{5};
Hdir = collectedVariables{6};

Udir = wrapTo360(Udir-180);
Hdir = wrapTo360(Hdir-180);

clearvars i commonTime variableNames numVars numTimes collectedVariables varData v TT


%% Reduce dataset to noNaN-only

% Step 1: Identify rows with NaNs
rowsWithNaNs = any(isnan(Ud), 2);

% Step 2: Remove these rows
t_clean = t(~rowsWithNaNs, :);
z_clean = z(~rowsWithNaNs, :);
Uz_clean = Uz(~rowsWithNaNs, :);
Ud_clean = Ud(~rowsWithNaNs, :);
Udir_clean = Udir(~rowsWithNaNs, :);
Hm0_clean = Hm0(~rowsWithNaNs, :);
Hdir_clean = Hdir(~rowsWithNaNs, :);

% Create timetable
T_clean = table(t_clean, z_clean, Uz_clean, Ud_clean, Udir_clean, Hm0_clean, Hdir_clean);
TT_clean = table2timetable(T_clean);

clearvars rowsWithNaNs T t_clean z_clean Uz_clean Ud_clean Udir_clean Hm0_clean Hdir_clean t z Uz Ud Udir Hm0 Hdir T_clean


%% Check timeseries
% figureRH;
% tiledlayout(2,1)
% 
% nexttile
% scatter(TT.time(1000:end,:), TT.Ud(1000:end,:), 200)
% ylabel('U (m/s)')
% legend(instru, "Location","northeastoutside")
% 
% nexttile
% scatter(TT.time(1000:end,:), TT.Hm0(1000:end,:), 200)
% ylabel('H_{m0} (m)')
% xlabel('minutes')
% 
% zoom xon


%% Longshore flow roses

% Preallocate a 1x6 array of Axes objects
ax = gobjects(1, size(TT_clean.Ud_clean, 2));

f1 = figure(Position=[740 1665 1719 628]);
tiledlayout(1, 6, "TileSpacing","compact");
for i = 1:size(TT_clean.Ud_clean, 2)
    ax(i) = nexttile; plot([],[]);
    Options = flowRoseOptions(fontsize, i);
    WindRose(TT_clean.Udir_clean(:,i), TT_clean.Ud_clean(:,i), Options)
    % camroll(-46)  % Rotate diagrams to match contour map orientation

    % plot coastline orientations
    line([cos(coastAngleN)*.1, cos(coastAngleN)], [sin(coastAngleN)*.1, sin(coastAngleN)], 'Color','k', 'LineStyle','-.', 'LineWidth',2, 'DisplayName','N spit')
    line([-cos(coastAngleN), -cos(coastAngleN)*.1], [-sin(coastAngleN), -sin(coastAngleN)*.1], 'Color','k', 'LineStyle','-.', 'LineWidth',2, 'HandleVisibility','off')

    line([cos(coastAngleS)*.1, cos(coastAngleS)], [sin(coastAngleS)*.1, sin(coastAngleS)], 'Color','k', 'LineStyle',':', 'LineWidth',2, 'DisplayName','S beach')
    line([-cos(coastAngleS), -cos(coastAngleS)*.1], [-sin(coastAngleS), -sin(coastAngleS)*.1], 'Color','k', 'LineStyle',':', 'LineWidth',2, 'HandleVisibility','off') 
end

clearvars i


%% Longshore wave roses

% Preallocate a 1x6 array of Axes objects
ax = gobjects(1, size(TT_clean.Hm0_clean, 2));

f2 = figure(Position=[740 957 1719 628]);
tiledlayout(1, 6, "TileSpacing","compact");
for i = 1:size(TT_clean.Hm0_clean, 2)
    ax(i) = nexttile; plot([],[]);
    Options = waveRoseOptions(fontsize, i);
    WindRose(TT_clean.Hdir_clean(:,i), TT_clean.Hm0_clean(:,i), Options)
    % camroll(-46)  % Rotate diagrams to match contour map orientation

    % plot coastline orientations
    line([cos(coastAngleN)*.1, cos(coastAngleN)], [sin(coastAngleN)*.1, sin(coastAngleN)], 'Color','k', 'LineStyle','-.', 'LineWidth',2, 'DisplayName','N spit')
    line([-cos(coastAngleN), -cos(coastAngleN)*.1], [-sin(coastAngleN), -sin(coastAngleN)*.1], 'Color','k', 'LineStyle','-.', 'LineWidth',2, 'HandleVisibility','off')

    line([cos(coastAngleS)*.1, cos(coastAngleS)], [sin(coastAngleS)*.1, sin(coastAngleS)], 'Color','k', 'LineStyle',':', 'LineWidth',2, 'DisplayName','S beach')
    line([-cos(coastAngleS), -cos(coastAngleS)*.1], [-sin(coastAngleS), -sin(coastAngleS)*.1], 'Color','k', 'LineStyle',':', 'LineWidth',2, 'HandleVisibility','off')
end

clearvars i
