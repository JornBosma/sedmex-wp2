%% Initialization
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;

dataPath = fullfile(filesep, 'Volumes', 'T7 Shield', 'DataDescriptor', 'hydrodynamics', 'ADV', filesep);

instruments = ["L1C1VEC", "L2C4VEC", "L3C1VEC", "L4C1VEC", "L5C1VEC", "L6C1VEC"];

newcolors = [cbf.vermilion; cbf.blue; cbf.bluegreen; cbf.yellow; cbf.redpurp; cbf.skyblue; cbf.orange];
newsymbols = {'o', 's', '^', 'd', 'v', 'p'};

k_sc = 0.05;  % Nikuradse current-related roughness [m]

% Defined starting time
t0 = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss');

% Load measured height of control volume above bed
hab_measured = load('hab.mat');
hab_measured = hab_measured.hab./100; % [cm] to [m]

% Iterate over each variable in the timetable
for varName = hab_measured.Properties.VariableNames
    % Get the variable data
    data = hab_measured.(varName{1});
    
    % Apply the condition: if value < 0.1, set it to 0.1
    data(data < 0.1) = 0.1;
    
    % Update the variable in the timetable
    hab_measured.(varName{1}) = data;
end

for i = 1:length(instruments)
filename = [dataPath char(instruments(i)) filesep 'tailored_' char(instruments(i)) '.nc'];
    
    % Load ADV data
    info = ncinfo(filename);
    t = ncread(filename, 't');              % seconds since 2021-09-01
    timeAxis = t0 + seconds(t);

    % Interpolate the measurements to the new time vector
    habCVnew = retime(hab_measured, timeAxis, 'pchip');
    habCVnew = habCVnew.(i);
    
    habCV = ncread(filename, 'h');         % height control volume above the bed [m]
    U_z = ncread(filename, 'umag');         % flow velocity [m/s] at depth z
    waterDepth = ncread(filename, 'd');    % water depth [m]
    waterLevel = ncread(filename, 'zs');   % water level [NAP+m]
    bedLevel = ncread(filename, 'zb');     % bed level [NAP+m]
    positionCV = ncread(filename, 'zi');   % position control volume [NAP+m]

    % Estimate the depth-averaged current velocity
    [U_c, ~] = compute_DAV(abs(U_z), habCVnew, k_sc, waterDepth);

    TT = timetable(timeAxis, bedLevel, waterLevel, positionCV, habCV, habCVnew, waterDepth, U_z, U_c);

    % Check new measuring heights
    figureRH;
    yyaxis left
    plot(TT.timeAxis, TT.habCVnew, '-b', 'LineWidth',3); hold on
    plot(TT.timeAxis, TT.habCV, '--k', 'LineWidth',3); hold off
    ylim([0, .4])
    ylabel('hab (m)')

    yyaxis right
    plot(TT.timeAxis, U_z, 'r', 'LineWidth',1); hold on
    plot(TT.timeAxis, U_c, 'y', 'LineWidth',1); hold off
    ylim([0, 1])
    ylabel('velocity (m/s)')

end

clearvars varName i data filename info t timeAxis habCVnew habCV U_z U_c waterDepth waterLevel bedLevel positionCV t0 ans


%% Visualisation
figure
tiledlayout('flow')

for i = 1:width(TT)
    nexttile
    plot(TT.timeAxis, TT.(i))
    title(TT.Properties.VariableNames{i})
    zoom xon
end