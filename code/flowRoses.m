%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;
% fontsize = 22; % ultra-wide screen

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];

% Coastline orientations
coastAngleNdeg = 50; % coastline angle northern beach [deg North]
coastAngleSdeg = 40; % coastline angle southern beach [deg North]

coastAngleN = deg2rad(90-coastAngleNdeg); % coastline angle northern beach [rad wrt East]
coastAngleS = deg2rad(90-coastAngleSdeg); % coastline angle southern beach [rad wrt East]

% angleDeg = -(rad2deg(angleRad)-90);


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
adv_names = {'L6C1VEC', 'L5C1VEC', 'L4C1VEC', 'L3C1VEC', 'L2C10VEC', 'L1C1VEC'};
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');
k_sc = 0.05;  % Nikuradse current-related roughness [m] (between 5 and 10 cm)

for n = 1:length(adv_names)
    % if n == 5
    %     ADVpath = [dataPath 'ADV' filesep 'L2C5SONTEK1' filesep 'tailored_L2C5SONTEK1.nc'];
    % else
        ADVpath = [dataPath 'ADV' filesep adv_names{n} filesep 'tailored_' adv_names{n} '.nc'];
    % end
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
    z = habCVnew.(adv_names{n});

    % Estimate the depth-averaged current velocity
    [Ud, ~] = compute_DAV(abs(Uz), z, k_sc, h);
    [Udl, ~] = compute_DAV(abs(Ulong), z, k_sc, h);
    [Udc, ~] = compute_DAV(abs(Ucros), z, k_sc, h);

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

    TT{n} = timetable(time, z, Uz, Ud, Udl, Udc, Udir, Hm0, Hdir);
end

clearvars t0 t info h Uz Udir Hm0 Hdir time habCVnew z Ud Udl Udc n hab_measured dataPath ADVpath Ulong Ucros


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
Udl = collectedVariables{4};
Udc = collectedVariables{5};
Udir = collectedVariables{6};
Hm0 = collectedVariables{7};
Hdir = collectedVariables{8};

clearvars i commonTime variableNames numVars numTimes collectedVariables varData v


%% Reduce dataset to noNaN-only

% Step 1: Identify rows with NaNs
rowsWithNaNs = any(isnan(Ud), 2);

% Step 2: Remove these rows
t_clean = t(~rowsWithNaNs, :);
z_clean = z(~rowsWithNaNs, :);
Uz_clean = Uz(~rowsWithNaNs, :);
Ud_clean = Ud(~rowsWithNaNs, :);
Udl_clean = Udl(~rowsWithNaNs, :);
Udc_clean = Udc(~rowsWithNaNs, :);
Udir_clean = Udir(~rowsWithNaNs, :);
Hm0_clean = Hm0(~rowsWithNaNs, :);
Hdir_clean = Hdir(~rowsWithNaNs, :);

% Create timetable
T_clean = table(t_clean, z_clean, Uz_clean, Ud_clean, Udl_clean, Udc_clean, Udir_clean, Hm0_clean, Hdir_clean);
TT_clean = table2timetable(T_clean);

clearvars rowsWithNaNs T t_clean z_clean Uz_clean Ud_clean Udir_clean Hm0_clean Hdir_clean t z Uz Ud Udc Udl Udir Hm0 Hdir Udc_clean Udl_clean


%% Check timeseries
% figureRH;
% tiledlayout(2,1)
% 
% nexttile
% scatter(TT.time(1000:end,:), TT.Ud(1000:end,:), 200)
% ylabel('U (m/s)')
% legend(adv_names, "Location","northeastoutside")
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


%% Additional calculations (not correct!)

% Wrap wave directions to [0, 360) and convert to arrival direction
waveDir = wrapTo360(-(TT_clean.Hdir_clean-90-180));
mean(waveDir, 'omitmissing')
mean(waveDir, 'all', 'omitmissing')

% Calculate shore normals
shoreNormalS = coastAngleSdeg + 90;
shoreNormalN = coastAngleNdeg + 90;

% Initialize waveDir_shoreNormal with NaN values
waveDir_shoreNormal = nan(size(waveDir));

% Calculate angle differences
waveDir_shoreNormal(:, [1, 2]) = wrapTo360(abs(shoreNormalS - waveDir(:, [1, 2])));
waveDir_shoreNormal(:, 3:end) = wrapTo360(abs(shoreNormalN - waveDir(:, 3:end)));

% Find the smallest angle by ensuring it lies within [0, 180] degrees
waveDir_shoreNormal = min(waveDir_shoreNormal, 360 - waveDir_shoreNormal);

% Calculate the mean wave angle wrt shore normal
mean(waveDir_shoreNormal, 'omitmissing')

%%

% Cross-shore v. longshore current
ratio_crossLong = mean(TT_clean.Udc_clean, 'omitmissing') ./ mean(TT_clean.Udl_clean, 'omitmissing') .* 100
ratio_crossLong_mean = mean(ratio_crossLong)

% Mean wave height
meanHm0 = mean(TT_clean.Hm0_clean, 'omitmissing')
meanHm0S = mean(meanHm0(1:2))
meanHm0N = mean(meanHm0(3:end))
ratio_Hm0_SN = meanHm0S/meanHm0N*100
100 - ratio_Hm0_SN

% Mean storm wave height
startDate = datetime(2021, 9, 27);
endDate = datetime(2021, 10, 6);
stormTT = TT_clean(TT_clean.t_clean >= startDate & TT_clean.t_clean <= endDate, :);
stormHm0 = stormTT.Hm0_clean;

meanHm0 = mean(stormHm0, 'omitmissing');
meanHm0S = mean(meanHm0(1:2))
meanHm0N = mean(meanHm0(3:end))
ratio_Hm0_SN = meanHm0S/meanHm0N*100
100 - ratio_Hm0_SN