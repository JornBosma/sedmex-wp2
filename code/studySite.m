%% Initialisation
close all
clear
clc

[~, fontsize, cbf, PHZ, ~] = sedmex_init;
% fontsize = 30; % ultra-wide screen

% Geographic limits NW Europe
lat1 = [44.1, 61.0];   % Latitude limits
lon1 = [-13.6, 22.7];  % Longitude limits

% Geographic limits NW Netherlands
lat2 = [52.6033, 53.4716];
lon2 = [4.0087, 5.6414];

% Geographic limits Prins Hendrikzanddijk
% lat3 = [PHZ.xLimWGS(1), PHZ.xLimWGS(2)];
% lon3 = [PHZ.yLimWGS(1), PHZ.yLimWGS(2)];
lat3 = [53.0006, 53.0323];
lon3 = [4.7812, 4.8385];


%% NW Europe
figure

% Load basemap data
geobasemap satellite

% Coordinates for the red edge (rectangular boundary)
edgeLat = [lat2(1), lat2(1), lat2(2), lat2(2), lat2(1)];
edgeLon = [lon2(1), lon2(2), lon2(2), lon2(1), lon2(1)];

% Plot a red edge
geoplot(edgeLat, edgeLon, 'r', 'LineWidth', 3);

% Set the limits
geolimits(lat1, lon1)

% Define north arrow parameters
arrowLength = 0.06; % Length of the arrow as a fraction of the figure height
arrowWidth = 0.2; % Width of the arrow as a fraction of the figure height
arrowPosition = [0.22, 0.8]; % Position of the arrow [x, y] in normalized figure coordinates

% Add north arrow
annotation('textarrow', ...
    'X', [arrowPosition(1), arrowPosition(1)], ...
    'Y', [arrowPosition(2), arrowPosition(2) + arrowLength], ...
    'String', 'N', ...
    'HeadWidth', arrowWidth * 100, ...
    'HeadLength', arrowLength * 100, ...
    'Color', 'w', ...
    'FontSize', fontsize, ...
    'LineWidth', 3);

% Set the figure properties
gp = gca;
gp.FontSize = fontsize;
gp.Scalebar.Visible = "on";
% gp.LatitudeLabel.String = '';
% gp.LongitudeLabel.String = '';


%% NW Netherlands
figure

% Load basemap data
geobasemap satellite

% Coordinates for the red edge (rectangular boundary)
edgeLat = [lat3(1), lat3(1), lat3(2), lat3(2), lat3(1)];
edgeLon = [lon3(1), lon3(2), lon3(2), lon3(1), lon3(1)];

% Plot a red edge
geoplot(edgeLat, edgeLon, 'r', 'LineWidth', 3);
hold on

% Add a red edge around the geoplot
% Coordinates for the red edge (rectangular boundary)
offset = 0.04;
edgeLat = [lat2(1), lat2(1), lat2(2), lat2(2), lat2(1)];
edgeLon = [lon2(1)-offset, lon2(2)+offset, lon2(2)+offset, lon2(1)-offset, lon2(1)-offset];

% Plot the red edge
geoplot(edgeLat, edgeLon, 'r', 'LineWidth', 5);
hold off

% Set the limits
geolimits(lat2, lon2)

% Define north arrow parameters
arrowLength = 0.06; % Length of the arrow as a fraction of the figure height
arrowWidth = 0.2; % Width of the arrow as a fraction of the figure height
arrowPosition = [0.28, 0.8]; % Position of the arrow [x, y] in normalized figure coordinates

% Add north arrow
annotation('textarrow', ...
    'X', [arrowPosition(1), arrowPosition(1)], ...
    'Y', [arrowPosition(2), arrowPosition(2) + arrowLength], ...
    'String', 'N', ...
    'HeadWidth', arrowWidth * 100, ...
    'HeadLength', arrowLength * 100, ...
    'Color', 'w', ...
    'FontSize', fontsize, ...
    'LineWidth', 3);

% Set the figure properties
gp = gca;
gp.FontSize = fontsize;
gp.Scalebar.Visible = "on";
% gp.LatitudeLabel.String = '';
% gp.LongitudeLabel.String = '';


%% Prins Hendrikzanddijk
figure

% Load basemap data
geobasemap satellite

% Add a red edge around the geoplot
% Coordinates for the red edge (rectangular boundary)
offset = 0;
edgeLat = [lat3(1), lat3(1), lat3(2), lat3(2), lat3(1)];
edgeLon = [lon3(1)-offset, lon3(2)+offset, lon3(2)+offset, lon3(1)-offset, lon3(1)-offset];

% Plot the red edge
geoplot(edgeLat, edgeLon, 'r', 'LineWidth', 5);

% Set the limits
geolimits(lat3, lon3)

% Define north arrow parameters
arrowLength = 0.06; % Length of the arrow as a fraction of the figure height
arrowWidth = 0.2; % Width of the arrow as a fraction of the figure height
arrowPosition = [0.34, 0.8]; % Position of the arrow [x, y] in normalized figure coordinates

% Add north arrow
annotation('textarrow', ...
    'X', [arrowPosition(1), arrowPosition(1)], ...
    'Y', [arrowPosition(2), arrowPosition(2) + arrowLength], ...
    'String', 'N', ...
    'HeadWidth', arrowWidth * 100, ...
    'HeadLength', arrowLength * 100, ...
    'Color', 'w', ...
    'FontSize', fontsize, ...
    'LineWidth', 3);

% Set the figure properties
gp = gca;
gp.FontSize = fontsize;
gp.Scalebar.Visible = "on";
% gp.LatitudeLabel.String = '';
% gp.LongitudeLabel.String = '';


%% Compilation
offset2 =  0.12;
offset3 =  0.005;

f1 = figure('Position',[615, 1751, 1844, 542]);
tiledlayout(1,3, 'TileSpacing','compact')

nexttile
% Coordinates for the red edge (NW Europe)
edgeLat = [lat2(1), lat2(1), lat2(2), lat2(2), lat2(1)];
edgeLon = [lon2(1)-offset2, lon2(2)+offset2, lon2(2)+offset2, lon2(1)-offset2, lon2(1)-offset2];

% Plot a red edge
geoplot(edgeLat, edgeLon, 'r', 'LineWidth', 4);

% Set the limits
geolimits(lat1, lon1)

% Add annotations
text(lat1(1)+12, lon1(1)+13, '\it{North Sea}', 'FontSize',fontsize*.8, 'Color','w', 'Rotation',0)

% Define north arrow parameters
arrowLength = 0.06; % Length of the arrow as a fraction of the figure height
arrowWidth = 0.2; % Width of the arrow as a fraction of the figure height
arrowPosition = [0.08, 0.8]; % Position of the arrow [x, y] in normalized figure coordinates

% Add north arrow
annotation('textarrow', ...
    'X', [arrowPosition(1), arrowPosition(1)], ...
    'Y', [arrowPosition(2), arrowPosition(2) + arrowLength], ...
    'String', 'N', ...
    'HeadWidth', arrowWidth * 100, ...
    'HeadLength', arrowLength * 100, ...
    'Color', 'w', ...
    'FontSize', fontsize, ...
    'LineWidth', 3);

% Set the figure properties
gp = gca;
gp.FontSize = fontsize;
gp.Scalebar.Visible = "on";
% gp.LatitudeLabel.String = '';
gp.LongitudeLabel.String = '';
gp.LatitudeAxis.TickLabelRotation = 90;
gp.GridColor = 'w';
gp.GridAlpha = .6;
gp.Scalebar.BackgroundAlpha = .8;

% Load basemap data
geobasemap satellite


nexttile
% Coordinates for the red edge (PHZ)
edgeLat = [lat3(1), lat3(1), lat3(2), lat3(2), lat3(1)];
edgeLon = [lon3(1)-offset3, lon3(2)+offset3, lon3(2)+offset3, lon3(1)-offset3, lon3(1)-offset3];

% Plot a red edge
geoplot(edgeLat, edgeLon, 'Color', [1 0 1], 'LineWidth', 3);
hold on

% Add a red edge around the geoplot
% Coordinates for the red edge (NW Europe)
edgeLat = [lat2(1), lat2(1), lat2(2), lat2(2), lat2(1)];
edgeLon = [lon2(1)-offset2, lon2(2)+offset2, lon2(2)+offset2, lon2(1)-offset2, lon2(1)-offset2];

% Plot the red edge
geoplot(edgeLat, edgeLon, 'Color', [1 0 0 .5], 'LineWidth', 10);
hold off

% Set the limits
geolimits(lat2, lon2)

% Add annotations
text(lat2(1)+0.5, lon2(1)+0.06, '\it{North Sea}', 'FontSize',fontsize*.8, 'Color','w', 'Rotation',0)
text(lat2(1)+0.4, lon2(1)+0.94, '\it{Wadden Sea}', 'FontSize',fontsize*.8, 'Color','w', 'Rotation',45)

% Define north arrow parameters
arrowLength = 0.06; % Length of the arrow as a fraction of the figure height
arrowWidth = 0.2; % Width of the arrow as a fraction of the figure height
arrowPosition = [0.39, 0.8]; % Position of the arrow [x, y] in normalized figure coordinates

% Add north arrow
annotation('textarrow', ...
    'X', [arrowPosition(1), arrowPosition(1)], ...
    'Y', [arrowPosition(2), arrowPosition(2) + arrowLength], ...
    'String', 'N', ...
    'HeadWidth', arrowWidth * 100, ...
    'HeadLength', arrowLength * 100, ...
    'Color', 'w', ...
    'FontSize', fontsize, ...
    'LineWidth', 3);

% Set the figure properties
gp = gca;
gp.FontSize = fontsize;
gp.Scalebar.Visible = "on";
gp.LatitudeLabel.String = '';
% gp.LongitudeLabel.String = '';
gp.LatitudeAxis.TickLabelRotation = 90;
gp.GridColor = 'w';
gp.GridAlpha = .6;
gp.Scalebar.BackgroundAlpha = .8;

% Load basemap data
geobasemap satellite


nexttile
% Add a red edge around the geoplot
% Coordinates for the red edge (PHZ)
edgeLat = [lat3(1), lat3(1), lat3(2), lat3(2), lat3(1)];
edgeLon = [lon3(1)-offset3, lon3(2)+offset3, lon3(2)+offset3, lon3(1)-offset3, lon3(1)-offset3];

% Plot the red edge
geoplot(edgeLat, edgeLon, 'Color', [1 0 1 .5], 'LineWidth', 10);

% Set the limits
geolimits(lat3, lon3)

% Add annotations
text(lat3(1)+0.01, lon3(1)+0.011, 'Prins Hendrikzanddijk', 'FontSize',fontsize, 'Color','w', 'Rotation',53.5)
text(lat3(1)+0.006, lon3(1)+0.035, '\it{Texelstroom}', 'FontSize',fontsize*.8, 'Color','w', 'Rotation',45)

% Define north arrow parameters
arrowLength = 0.06; % Length of the arrow as a fraction of the figure height
arrowWidth = 0.2; % Width of the arrow as a fraction of the figure height
arrowPosition = [0.7, 0.8]; % Position of the arrow [x, y] in normalized figure coordinates

% Add north arrow
annotation('textarrow', ...
    'X', [arrowPosition(1), arrowPosition(1)], ...
    'Y', [arrowPosition(2), arrowPosition(2) + arrowLength], ...
    'String', 'N', ...
    'HeadWidth', arrowWidth * 100, ...
    'HeadLength', arrowLength * 100, ...
    'Color', 'w', ...
    'FontSize', fontsize, ...
    'LineWidth', 3);

% Set the figure properties
gp = gca;
gp.FontSize = fontsize;
gp.Scalebar.Visible = "on";
gp.LatitudeLabel.String = '';
gp.LongitudeLabel.String = '';
gp.LatitudeAxis.TickLabelRotation = 79;
gp.GridColor = 'w';
gp.GridAlpha = .6;
gp.Scalebar.BackgroundAlpha = .8;

% Load basemap data
geobasemap satellite
