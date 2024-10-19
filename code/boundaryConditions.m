%% Initialisation
close all
clear
clc

[~, fontsize, cbf, PHZ, SEDMEX] = sedmex_init;
% fontsize = 30; % ultra-wide screenSEDMEXtime = [datetime('10-Sep-2021'), datetime('19-Oct-2021')];

dataPath = [filesep 'Volumes' filesep 'T7 Shield' filesep 'DataDescriptor'...
    filesep 'hydrodynamics' filesep 'ADV' filesep 'L2C10VEC' filesep' 'tailored_L2C10VEC.nc'];

XsamplingDates = [datetime('20-Sep-2021'); datetime('28-Sep-2021');...
    datetime('01-Oct-2021'); datetime('07-Oct-2021'); datetime('15-Oct-2021');...
    datetime('21-Sep-2021'); datetime('28-Sep-2021'); datetime('03-Oct-2021');...
    datetime('08-Oct-2021')];

YsamplingDates = [datetime('21-Sep-2021'); datetime('28-Sep-2021');...
    datetime('03-Oct-2021'); datetime('08-Oct-2021')];

stormy_period = [datetime('27-Sep-2021 12:00'), datetime('05-Oct-2021 18:00')];
storm_period = [datetime('30-Sep-2021 16:00'), datetime('03-Oct-2021 03:00')];


%% Monitoring period

% % load KNMI De Kooy weather station data
% DeKooy = readtable('DeKooy2019_2022_hourly.txt', 'VariableNamingRule','preserve');
% DeKooy.("# STN") = []; % superfluous
% 
% DeKooy = renamevars(DeKooy,'YYYYMMDD','date');
% DeKooy = renamevars(DeKooy,'H','time');
% DeKooy = renamevars(DeKooy,'DD','dir'); % mean wind direction (in degrees) during the 10-minute period preceding the time of observation (360=north; 90=east; 180=south; 270=west; 0=calm 990=variable)
% DeKooy = renamevars(DeKooy,'FF','spd'); % mean wind speed (in 0.1 m/s) during the 10-minute period preceding the time of observation
% DeKooy = renamevars(DeKooy,'T','temp'); % temperature in 0.1 degrees Celsius at 1.50 m
% 
% DeKooy.date = datetime(DeKooy.date, 'ConvertFrom','yyyyMMDD');
% DeKooy.time = hours(DeKooy.time);
% DeKooy.DateTime = DeKooy.date + DeKooy.time;
% DeKooy = removevars(DeKooy, {'date','time'});
% 
% DeKooy.dir(DeKooy.dir == 0) = NaN; % code for no wind
% DeKooy.spd(DeKooy.dir == 0) = NaN;
% DeKooy.dir(DeKooy.dir == 990) = NaN; % code for variable wind
% DeKooy.spd(DeKooy.dir == 990) = NaN;
% 
% DeKooy.spd = DeKooy.spd/10; % convert 0.1 m/s to m/s
% DeKooy.temp = DeKooy.temp/10; % convert 0.1 deg C to deg C
% TTwind = table2timetable(DeKooy);

% load and prepare water-level data
Oudeschild = readtable('20230904_029.csv', 'Range','V1:Y213210', 'VariableNamingRule','preserve');
Oudeschild.('LIMIETSYMBOOL') = []; % superfluous
Oudeschild.NUMERIEKEWAARDE(Oudeschild.NUMERIEKEWAARDE>=999) = NaN; % error code

Oudeschild = renamevars(Oudeschild,'WAARNEMINGDATUM','date');
Oudeschild = renamevars(Oudeschild,'WAARNEMINGTIJD (MET/CET)','time');
Oudeschild = renamevars(Oudeschild,'NUMERIEKEWAARDE','eta');

Oudeschild.date = datetime(Oudeschild.date, 'InputFormat','dd-MM-yyyy');
Oudeschild.DateTime = Oudeschild.date + Oudeschild.time;
Oudeschild = removevars(Oudeschild, {'date','time'});

Oudeschild.eta = Oudeschild.eta/100; % convert cm to m
TTwater = table2timetable(Oudeschild);

% Isolation
TThws = TTwater; % extract data above HWS threshold
TTlws = TTwater; % extract data below LWS threshold
TTmws = TTwater; % time series with values inside the bandwidth

TThws.eta(TTwater.eta <= PHZ.HWS) = NaN;
TTlws.eta(TTwater.eta >= PHZ.LWS) = NaN;
TTmws.eta(TTwater.eta > PHZ.HWS | TTwater.eta < PHZ.LWS) = NaN;

% Isolation (SEDMEX period)
TTmhw = TTwater; % extract data above MHW threshold during SEDMEX
TTmlw = TTwater; % extract data below MLW threshold during SEDMEX
TTmw = TTwater; % time series with values inside the bandwidth

TTmhw.eta(TTwater.eta <= SEDMEX.MeanHW) = NaN;
TTmlw.eta(TTwater.eta >= SEDMEX.MeanLW) = NaN;
TTmw.eta(TTwater.eta > SEDMEX.MeanHW | TTwater.eta < SEDMEX.MeanLW) = NaN;


%% SEDMEX
SEDMEXtime = [datetime('10-Sep-2021'), datetime('19-Oct-2021')];

info = ncinfo(dataPath);
t_seconds = ncread(dataPath, 't'); % minutes since 2021-09-10 00:00:00
t = datetime('2021-09-01 00:00:00','InputFormat','yyyy-MM-dd HH:mm:ss')+seconds(t_seconds);
eta = ncread(dataPath, 'zs'); % water depth [m]
umag = ncread(dataPath, 'umag'); % velocity magnitude [m/s]
ucm = ncread(dataPath, 'ucm'); % cross-shore velocity magnitude [m/s]
ulm = ncread(dataPath, 'ulm'); % longshore velocity magnitude [m/s]
Hm0 = ncread(dataPath, 'Hm0'); % significant wave height [m]
puvdir = ncread(dataPath, 'puvdir'); % wave propagation direction [deg]

ulm = -ulm;           % reverse the sign so that NE (flood) is positive

% Convert Cartesian direction convention into Nautical direction
% puvdir = wrapTo360(90-puvdir+180);

% Convert to the new coordinate system
puvdir = mod(-(puvdir - 90), 360); % Convert to positive clockwise with N = 0°
puvdir = mod(puvdir + 180, 360); % Convert to wave arrival direction (like wind direction)
puvdir = mod(puvdir + 360, 360); % Ensure angle is positive

% Calculate the threshold for deviation
threshold = 3 * std(puvdir, 'omitnan');

% Find the indices of values that deviate more than twice the standard deviation
indices = abs(puvdir - mean(puvdir, 'omitnan')) > threshold;

% Set those values to NaN
puvdir(indices) = NaN;

T_L2C10 = table(t, eta, umag, ucm, ulm, Hm0, puvdir);
TT_L2C10 = table2timetable(T_L2C10);

% Identify waves (SEDMEX period)
TTh_L2C10 = TT_L2C10; % extract data above mean+1std threshold during SEDMEX
TTl_L2C10 = TT_L2C10; % extract data below mean+1std threshold during SEDMEX

Hm0_threshold = mean(TT_L2C10.Hm0,'omitmissing')+2*std(TT_L2C10.Hm0,'omitmissing');
TTh_L2C10.Hm0(TT_L2C10.Hm0 <= Hm0_threshold) = NaN;
TTl_L2C10.Hm0(TT_L2C10.Hm0 > Hm0_threshold) = NaN;


%% Visualisation: SEDMEX boundary conditions
f1 = figure('Position',[902, 1434, 1470, 754]);
% tiledlayout(6,1, 'TileSpacing','loose')
tiledlayout(4,1, 'TileSpacing','loose')

% ax(1) = nexttile;
% plot(TTwind.DateTime, TTwind.spd, 'Color','k', 'LineWidth',3)
% ylabel('U_{wind} (m s^{-1})', 'FontSize',fontsize*.8)
% 
% ax(2) = nexttile;
% scatter(TTwind.DateTime, TTwind.dir, 20, 'k', 'filled'); hold on
% 
% % yline(45, '-', 'Color',cbf.redpurp, 'LineWidth',1) % longshore wind (northeastward)
% yline(45+90, '-', 'Color','r', 'LineWidth',1) % onshore wind
% % yline(45+180, '-', 'Color',cbf.redpurp, 'LineWidth',1) % longshore wind (southwestward)
% cross = fill([SEDMEXtime, fliplr(SEDMEXtime)], [[45+135, 45+135], [45+45, 45+45]], cbf.redpurp, 'FaceAlpha',0.1, 'LineStyle','none');
% longS = fill([SEDMEXtime, fliplr(SEDMEXtime)], [[45+180, 45+180], [45+135, 45+135]], cbf.vermilion, 'FaceAlpha',0.1, 'LineStyle','none');
% longN = fill([SEDMEXtime, fliplr(SEDMEXtime)], [[45, 45], [45+45, 45+45]], cbf.vermilion, 'FaceAlpha',0.1, 'LineStyle','none'); hold off
% 
% ylabel('\alpha_{wind} (\circ)', 'FontSize',fontsize*.8)
% lgd2 = legend([cross,longS],{'onshore','longshore'}, 'Position',[0.4564, 0.7681, 0.1632, 0.0240], 'NumColumns',2);

ax(3) = nexttile;
hold on
plot(TT_L2C10.t, TTl_L2C10.Hm0, 'Color','k', 'LineWidth',3)
plot(TT_L2C10.t, TTh_L2C10.Hm0, 'Color','r', 'LineWidth',3)

% Get the current y-axis limits
ylim_vals = ylim();

% Add downward-pointing triangle markers at the top, slightly above the upper y-axis limit
y_marker_position1 = ylim_vals(2) + 0.05 * (ylim_vals(2) - ylim_vals(1));  % 5% above the upper limit
y_marker_position2 = y_marker_position1 + 0.05 * (ylim_vals(2) - ylim_vals(1));  % 5% above y_marker_position1
plot(XsamplingDates, repmat(y_marker_position1, size(XsamplingDates)), 'v', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none')
plot(YsamplingDates, repmat(y_marker_position2, size(YsamplingDates)), 'v', 'MarkerSize', 10, ...
    'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'none')
hold off

yline(Hm0_threshold, '--', '            \mu+2σ', 'LineWidth',3,...
    'FontSize',fontsize*.8, 'LabelHorizontalAlignment','left')
ylabel('H_{m0} (m)', 'FontSize',fontsize*.8)
yticks(0:.2:.6)

ax(4) = nexttile;
scatter(TT_L2C10.t, TT_L2C10.puvdir, 10, 'k', 'filled'); hold on

% Create the background color gradient
normal_direction = 45+90;
n_points = 1000; % Number of points for the gradient
y_range = linspace(45, 225, n_points); % Assuming wind direction range from 45 to 225
alpha_values = 1 - abs(y_range - normal_direction) / max(y_range - normal_direction);

% Create the coloured background
for i = 1:length(y_range)-1
    y1 = y_range(i);
    y2 = y_range(i+1);
    alpha1 = alpha_values(i);
    alpha2 = alpha_values(i+1);
    
    fill([SEDMEXtime(1), SEDMEXtime(end), SEDMEXtime(end), SEDMEXtime(1)], [y1, y1, y2, y2], cbf.vermilion, ...
        'FaceAlpha',(alpha1 + alpha2)/2, 'EdgeColor','none');
end

% Redraw the scatter plot and yline to ensure they are on top
scatter(TT_L2C10.t, TT_L2C10.puvdir, 10, 'k', 'filled')
text(datetime("04-Oct-2021"), normal_direction, " onshore", 'FontSize',fontsize*.6, 'FontWeight','bold')
text(datetime("04-Oct-2021"), normal_direction-90 + 0, "longshore", 'FontSize',fontsize*.6, 'FontWeight','bold')
text(datetime("04-Oct-2021"), normal_direction+90 - 0, "longshore", 'FontSize',fontsize*.6, 'FontWeight','bold')

ylabel(ax(4), '\alpha_{wave} (\circN)', 'FontSize',fontsize*.8)

ax(5) = nexttile;
plot(TT_L2C10.t, TT_L2C10.ulm, 'Color','k', 'LineWidth',3); hold on
plot(TT_L2C10.t, TT_L2C10.ucm, 'Color','r', 'LineWidth',3); hold off
ylabel('U (m s^{-1})', 'FontSize',fontsize*.8)
lgd5 = legend('longshore', 'cross-shore', 'Position',[0.1250, 0.5140, 0.1864, 0.0240],...
    'NumColumns',2, 'FontSize',fontsize*.6);
text(SEDMEXtime(2)+hours(4), .3, 'NE', 'FontSize',fontsize*.6, 'FontWeight','bold')
text(SEDMEXtime(2)+hours(4), -.3, 'SW', 'FontSize',fontsize*.6, 'FontWeight','bold')
text(SEDMEXtime(2)-hours(9), .1, 'onshore', 'FontSize',fontsize*.6, 'Color','r', 'FontWeight','bold')
text(SEDMEXtime(2)-hours(9), -.1, 'offshore', 'FontSize',fontsize*.6, 'Color','r', 'FontWeight','bold')

ax(6) = nexttile;
plot(TTmw.DateTime, TTmw.eta, 'Color','k', 'LineWidth',3); hold on
plot(TTmhw.DateTime, TTmhw.eta, 'Color','r', 'LineWidth',3)
plot(TTmlw.DateTime, TTmlw.eta, 'Color','k', 'LineWidth',3); hold off
ylabel('η (NAP+m)', 'FontSize',fontsize*.8)

yline(SEDMEX.MeanHW, '--', 'LineWidth',3)
yline(SEDMEX.MeanLW, '--', 'LineWidth',3)
% yline(SEDMEX.MeanHW, '--', 'MHW', 'LineWidth',3, 'FontSize',fontsize*.8)
% yline(SEDMEX.MeanLW, '--', 'MLW', 'LineWidth',3, 'FontSize',fontsize*.8, 'LabelVerticalAlignment','bottom')
text(SEDMEXtime(2)+hours(4),SEDMEX.MeanHW, 'MHW', 'FontSize',fontsize*.6, 'FontWeight','bold')
text(SEDMEXtime(2)+hours(4),SEDMEX.MeanLW, 'MLW', 'FontSize',fontsize*.6, 'FontWeight','bold')

xlim([ax(3) ax(4) ax(5) ax(6)], [SEDMEXtime(1), SEDMEXtime(2)])
xticks([ax(3) ax(4) ax(5) ax(6)], SEDMEXtime(1):days(2):SEDMEXtime(2))
xticklabels([ax(3) ax(4) ax(5)], [])
xtickformat('MMM dd')
xtickangle(45)
grid([ax(3) ax(4) ax(5) ax(6)],'on')
linkaxes([ax(3) ax(4) ax(5) ax(6)], 'x')
zoom xon

% yticks(ax(4),'manual') % exporting figure changes yticks...
ylim(ax(4), [45, 225])
yticks(ax(4), 90:45:180)
ylim(ax(5), [-.4, .4])
yticks(ax(5), -.4:.2:.4)
% ylim(ax(6), [-1, 1.2])
yticks(ax(6), -1:.5:1.5)

% xlim([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6)], [SEDMEXtime(1), SEDMEXtime(2)])
% xticks([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6)], SEDMEXtime(1):days(2):SEDMEXtime(2))
% xticklabels([ax(1) ax(2) ax(3) ax(4) ax(5)], [])
% xtickformat('MMM dd')
% xtickangle(45)
% grid([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6)],'on')
% linkaxes([ax(1) ax(2) ax(3) ax(4) ax(5) ax(6)], 'x')
% zoom xon
% 
% ylim(ax(2), [0 360])
% ylim(ax(4), [45 225])
% yticks([ax(2), ax(4)], 0:90:360)

% Add figure annotations
annotation('textbox', [0.94, 0.87, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');
annotation('textbox', [0.94, 0.666, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');
annotation('textbox', [0.94, 0.464, 0.1, 0.1], 'String','(c)', 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');
annotation('textbox', [0.94, 0.261, 0.1, 0.1], 'String','(d)', 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');

% Add sediment sampling moments
for i = 3:6
    % xline(ax(i), XsamplingDates, '-', 'LineWidth',2, 'Color','b', 'HandleVisibility','off')
    % xline(ax(i), YsamplingDates, '--', 'LineWidth',2, 'Color','g', 'HandleVisibility','off')
    xregion(ax(i), stormy_period, 'HandleVisibility','off', 'FaceColor',[.7 .7 .7])
    xregion(ax(i), storm_period, 'HandleVisibility','off', 'FaceColor',[.4 .4 .4])
end


%% Additional calculations
h = ncread(dataPath, 'h'); % height control volume above bed [m]
h_mean = mean(h, 'omitmissing');

ratio = abs(TT_L2C10.ulm ./ TT_L2C10.ucm);
ratio_mean = mean(ratio, 'omitmissing');