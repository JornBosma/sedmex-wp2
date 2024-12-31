%% Initialisation
close all
clear
clc

[~, fontsize, cbf, ~, ~] = sedmex_init;
% fontsize = 26; % ultra-wide screen

folderPath = [filesep 'Volumes' filesep 'T7 Shield' filesep...
    'DataDescriptor' filesep 'grainsizes' filesep];

% sand scraper depth settings
depth_mm = [0, -2:-4:-14, -20:-6:-50]';


%% Load sediment data
dataPath{1} = [folderPath 'GS_20211007s.csv'];
dataPath{2} = [folderPath 'GS_20211015s.csv'];

opts = detectImportOptions(dataPath{1});
opts = setvaropts(opts,'Date_ddMMyyyy','InputFormat','dd/MM/yyyy');

GS_20211007 = readtable(dataPath{1}, opts);
GS_20211015 = readtable(dataPath{2}, opts);


%% Modify data arrays (Mean)

% 7 October
LW07 = GS_20211007{1:10, "Mean_mu"}/1000;  % convert to mm
LW07 = [LW07; LW07(end)];

MW07 = GS_20211007{11:20, "Mean_mu"}/1000;
MW07 = [MW07; MW07(end)];

HW07 = GS_20211007{21:30, "Mean_mu"}/1000;
HW07 = [HW07; HW07(end)];

% 15 October
LW15 = GS_20211015{1:10, "Mean_mu"}/1000;
LW15 = [LW15; LW15(end)];

MW15 = GS_20211015{11:20, "Mean_mu"}/1000;
MW15 = [MW15; MW15(end)];

HW15 = GS_20211015{21:30, "Mean_mu"}/1000;
HW15 = [HW15; HW15(end)];

bar15 = GS_20211015{31:40, "Mean_mu"}/1000;
bar15 = [bar15; bar15(end)];

runnel15 = GS_20211015{41:48, "Mean_mu"}/1000;
runnel15 = [runnel15; NaN; NaN; NaN];  % groundwater hit

% Calculate mean stats
Mg_mean = mean([HW07(1:end-1), MW07(1:end-1), LW07(1:end-1), HW15(1:end-1),...
    MW15(1:end-1), LW15(1:end-1), bar15(1:end-1), runnel15(1:end-1)], 'omitmissing');
Mg_std = std([HW07(1:end-1), MW07(1:end-1), LW07(1:end-1), HW15(1:end-1),...
    MW15(1:end-1), LW15(1:end-1), bar15(1:end-1), runnel15(1:end-1)], 'omitmissing');
Mg_rel = Mg_std./Mg_mean;


%% Visualisation: pcolor
% f1 = figure("Position",[860 1681 923 330]);
% tiledlayout(1,8, "TileSpacing","compact")

f = figure('Position',[740, 1775, 643, 518]);
tiledlayout(2,8, "TileSpacing","compact")

ax(1) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [HW07, HW07])
% title('HW7')
title('A1')
ylabel('depth (mm)')

ax(2) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [MW07, MW07])
% title('MW7')
title('A2')

ax(3) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [LW07, LW07])
% title('LW7')
title('A3')

ax(4) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [HW15, HW15])
% title('HW15')
title('B1')

ax(5) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [MW15, MW15])
% title('MW15')
title('B2')

ax(6) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [LW15, LW15])
% title('LW15')
title('B3')

ax(7) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [bar15, bar15])
title('Bar')

ax(8) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [runnel15, runnel15])
% title('runnel')
title('Run')
ax(8).Color = cbf.skyblue;
yticks({})

c = colorbar;
% c.Position = [0.9 0.11 0.013 0.8];
c.Label.String = 'M_{G} (mm)';
c.FontSize = fontsize*.8;

set(ax(2:8), 'yTickLabel', [])
set(ax, 'YDir', 'normal')
set(ax,'xTick', [])
% set(ax, 'Colormap', crameri('lajolla'), 'CLim',[0.4 1.6])
set(ax, 'Colormap', brewermap([],"YlOrRd"), 'CLim',[0.4 1.6])

% for i = 1:length(ax)
%     shading(ax(i), 'interp')
% end


%% Modify data arrays (Sorting)

% 7 October
LW07 = GS_20211007{1:10, "Sorting"};
LW07 = [LW07; LW07(end)];

MW07 = GS_20211007{11:20, "Sorting"};
MW07 = [MW07; MW07(end)];

HW07 = GS_20211007{21:30, "Sorting"};
HW07 = [HW07; HW07(end)];

% 15 October
LW15 = GS_20211015{1:10, "Sorting"};
LW15 = [LW15; LW15(end)];

MW15 = GS_20211015{11:20, "Sorting"};
MW15 = [MW15; MW15(end)];

HW15 = GS_20211015{21:30, "Sorting"};
HW15 = [HW15; HW15(end)];

bar15 = GS_20211015{31:40, "Sorting"};
bar15 = [bar15; bar15(end)];

runnel15 = GS_20211015{41:48, "Sorting"};
runnel15 = [runnel15; NaN; NaN; NaN];  % groundwater hit

% Calculate mean stats
Sg_mean = mean([HW07(1:end-1), MW07(1:end-1), LW07(1:end-1), HW15(1:end-1),...
    MW15(1:end-1), LW15(1:end-1), bar15(1:end-1), runnel15(1:end-1)], 'omitmissing');
Sg_std = std([HW07(1:end-1), MW07(1:end-1), LW07(1:end-1), HW15(1:end-1),...
    MW15(1:end-1), LW15(1:end-1), bar15(1:end-1), runnel15(1:end-1)], 'omitmissing');
Sg_rel = Sg_std./Sg_mean;


%% Visualisation: pcolor
% f2 = figure("Position",[860 1271 923 330]);
% tiledlayout(1,8, "TileSpacing","compact")

ax(1) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [HW07, HW07])
% title('HW7')
ylabel('depth (mm)')

ax(2) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [MW07, MW07])
% title('MW7')

ax(3) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [LW07, LW07])
% title('LW7')

ax(4) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [HW15, HW15])
% title('HW15')

ax(5) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [MW15, MW15])
% title('MW15')

ax(6) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [LW15, LW15])
% title('LW15')

ax(7) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [bar15, bar15])
% title('bar')

ax(8) = nexttile;
pcolor([1, 2], [depth_mm, depth_mm], [runnel15, runnel15])
% title('runnel')
ax(8).Color = cbf.skyblue;
yticks({})

c = colorbar;
% c.Position = [0.9 0.11 0.013 0.8];
c.Label.String = 'σ_{G}';
c.FontSize = fontsize*.8;

set(ax(2:8), 'yTickLabel', [])
set(ax, 'YDir', 'normal')
set(ax,'xTick', [])
% set(ax, 'Colormap', crameri('imola'), 'CLim',[1 3])
set(ax, 'Colormap', brewermap([],"BuPu"), 'CLim',[1.5 3])

% for i = 1:length(ax)
%     shading(ax(i), 'interp')
% end

% Add figure annotations
annotation('textbox', [0.02, 0.88, 0.1, 0.1], 'String','(a)', 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');
annotation('textbox', [0.02, 0.43, 0.1, 0.1], 'String','(b)', 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','bold');

% annotation('textbox', [0.2, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(1)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.29, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(2)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.372, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(3)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.462, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(4)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.55, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(5)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.64, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(6)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.73, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(7)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.82, 0.65, 0.1, 0.1], 'String', sprintf('%.2f', Mg_mean(8)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');

% annotation('textbox', [0.2, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(1)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.29, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(2)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.372, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(3)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.462, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(4)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.55, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(5)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.64, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(6)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.73, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(7)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');
% annotation('textbox', [0.82, 0.21, 0.1, 0.1], 'String', sprintf('%.2f', Sg_mean(8)), 'EdgeColor','none', 'FontSize',fontsize*.8, 'FontWeight','normal', 'Rotation',90, 'Color','w');