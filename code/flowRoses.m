%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;
% fontsize = 22; % ultra-wide screen

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];


%% Load ADV data
adv_locs = {'L6C1', 'L5C1', 'L4C1', 'L3C1', 'L2C10', 'L1C1'};

% t_start = datetime('2021-09-01 00:00:00');

for n = 1:length(adv_locs)
    % if n == 5
    %     ADVpath{n} = [dataPath 'ADV' filesep adv_locs{n} 'SONTEK1' filesep 'tailored_' adv_locs{n} 'SONTEK1.nc'];
    % else
        ADVpath{n} = [dataPath 'ADV' filesep adv_locs{n} 'VEC' filesep 'tailored_' adv_locs{n} 'VEC.nc'];
    % end
    t_adv{n} = ncread(ADVpath{n}, 't'); % seconds since 2021-09-01 00:00:00
    info_adv.(adv_locs{n}) = ncinfo(ADVpath{n});
    hab_adv{n} = ncread(ADVpath{n}, 'h');
    U_adv{n} = ncread(ADVpath{n}, 'umag');
    Uang_adv{n} = ncread(ADVpath{n}, 'uang');
    Hm0_adv{n} = ncread(ADVpath{n}, 'Hm0');
    Hdir_adv{n} = ncread(ADVpath{n}, 'puvdir');

    % % Convert Cartesian direction convention into Nautical direction
    % Uang_adv{n} = wrapTo360(90-Uang_adv{n});
    % Hdir_adv{n} = wrapTo360(90-Hdir_adv{n});

    % Checks if measurement volume is not buried
    buried1{n} = hab_adv{n} < .05;
    buried2{n} = U_adv{n} < .002;

    t_adv{n}(buried1{n} | buried2{n}) = NaN;
    U_adv{n}(buried1{n} | buried2{n}) = NaN;
    Uang_adv{n}(buried1{n} | buried2{n}) = NaN;
    Hm0_adv{n}(buried1{n} | buried2{n}) = NaN;
    Hdir_adv{n}(buried1{n} | buried2{n}) = NaN;
end

% Hm0_adv{5} = mean([Hm0_adv{4}; Hm0_adv{6}]);


%% Select useful data from equal periods

% Step 1: Find the length of the shortest array
minLength = min(cellfun(@length, U_adv));

% Step 2: Truncate each array in the cell array
t_adv = cellfun(@(x) x(1:minLength), t_adv, 'UniformOutput',false);
U_adv = cellfun(@(x) x(1:minLength), U_adv, 'UniformOutput',false);
Uang_adv = cellfun(@(x) x(1:minLength), Uang_adv, 'UniformOutput',false);
hab_adv = cellfun(@(x) x(1:minLength), hab_adv, 'UniformOutput',false);
Hm0_adv = cellfun(@(x) x(1:minLength), Hm0_adv, 'UniformOutput',false);
% Hm0_adv{:,5} = Hm0_adv{:,5}';
Hdir_adv = cellfun(@(x) x(1:minLength), Hdir_adv, 'UniformOutput',false);

t = cell2mat(t_adv);
U = cell2mat(U_adv);
Uang = cell2mat(Uang_adv);
hab = cell2mat(hab_adv);
Hm0 = cell2mat(Hm0_adv);
Hdir =  cell2mat(Hdir_adv);

% Step 1: Identify rows with NaNs
rowsWithNaNs = any(isnan(U), 2);

% Step 2: Remove these rows
t_clean = t(~rowsWithNaNs, :);
U_clean = U(~rowsWithNaNs, :);
Uang_clean = Uang(~rowsWithNaNs, :);
hab_clean = hab(~rowsWithNaNs, :);
Hm0_clean = Hm0(~rowsWithNaNs, :);
% Hm0_clean(:, 5) = mean([Hm0_clean(:, 4), Hm0_clean(:, 6)], 2);
Hdir_clean = Hdir(~rowsWithNaNs, :);
% Hdir_clean(:, 5) = rad2deg(Hdir_clean(:, 5));
% Hm0_clean(:, 5) = Hm0_clean(:, 5).*10;


%% Check 'time series'
figureRH;
tiledlayout(2,1)

nexttile
scatter(t_clean(1000:end,:), U_clean(1000:end,:), 200)
ylabel('U (m/s)')
legend(adv_locs, "Location","northeastoutside")

nexttile
scatter(t_clean(1000:end,:), Hm0_clean(1000:end,:), 200)
ylabel('H_{m0} (m)')
xlabel('minutes')

zoom xon


%% Longshore flow roses

% Preallocate a 1x6 array of Axes objects
ax = gobjects(1, size(U_clean, 2));

f1 = figure(Position=[740 1665 1719 628]);
tiledlayout(1, 6, "TileSpacing","compact");
for i = 1:size(U_clean, 2)
    ax(i) = nexttile; plot([],[]);
    Options = flowRoseOptions(fontsize, i);
    WindRose(Uang_clean(:,i), U_clean(:,i), Options)
    % camroll(-46)  % Rotate diagrams to match contour map orientation
end


%% Longshore wave roses

% Preallocate a 1x6 array of Axes objects
ax = gobjects(1, size(Hm0_clean, 2));

f2 = figure(Position=[740 957 1719 628]);
tiledlayout(1, 6, "TileSpacing","compact");
for i = 1:size(Hm0_clean, 2)
    ax(i) = nexttile; plot([],[]);
    Options = waveRoseOptions(fontsize, i);
    WindRose(Hdir_clean(:,i), Hm0_clean(:,i), Options)
    % camroll(-46)  % Rotate diagrams to match contour map orientation
end

