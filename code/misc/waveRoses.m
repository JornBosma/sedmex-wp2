%% Initialisation
close all
clear
clc

[~, fontsize, ~, ~, ~] = sedmex_init;
% fontsize = 38; % ultra-wide screen

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor' filesep 'hydrodynamics' filesep];

locs = {'L6', 'L5', 'L4', 'L3', 'L2', 'L1'};

locs_extended_ossi = locs(1:3);
locs_extended_ossi{4} = '';
locs_extended_ossi{5} = '';
locs_extended_ossi{6} = '';
locs_extended_ossi{7} = 'L2';
locs_extended_ossi{8} = 'L1';


%% Load OSSI data
ossi_locs = {'L6C2', 'L5C2', 'L4C3', 'L2C9', 'L1C2'};

% t_start = datetime('2021-09-01 00:00:00');

for n = 1:length(ossi_locs)
    OSSIpath{n} = [dataPath 'pressuresensors' filesep ossi_locs{n} 'OSSI' filesep 'tailored_' ossi_locs{n} 'OSSI.nc'];
    info_ossi.(ossi_locs{n}) = ncinfo(OSSIpath{n});
    t_ossi{n} = ncread(OSSIpath{n}, 't'); % minutes since 2021-09-10 19:00:00
    zb_ossi{n} = ncread(OSSIpath{n}, 'zb');
    Hm0_ossi{n} = ncread(OSSIpath{n}, 'Hm0');
end

Hm0 = cell2mat(Hm0_ossi);


%% Check 'time series'
figureRH;
plot(Hm0(1000:end,:), 'LineWidth',2)
legend(ossi_locs, "Location","northeastoutside")
zoom xon


%% Box plot: longshore distribution time-average wave height
Hm0_box = [Hm0(:,1:3), NaN(length(Hm0),3), Hm0(:,4:5)];

f1 = figure('Position',[1722 892 1719 445]);
boxchart(Hm0_box, 'MarkerStyle','none', 'BoxWidth',.2)  % no outliers
set(gca,'xticklabel',[])
ylabel('H_{m0} (m)')
ylim([0 .42])
get(gca, 'XTick');
xticklabels(locs_extended_ossi)
set(gca, 'TickLength', [0 0]);  % Hide x-tick lines

