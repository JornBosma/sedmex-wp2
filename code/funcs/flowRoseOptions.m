function Options = flowRoseOptions(fontsize, tilenum)

% mycolormap = crameri('-roma');
mycolormap = brewermap([],"-RdYlBu");
% mycolormap = colormap(hsv(6));
% mycolormap = [0 0 1; 0 1 1; 0 1 0; .5 .2 .1; 1 1 0; 1 0 0];

Options.axes = gca;

Options.AngleNorth     = 90;
Options.AngleEast      = 0;
% Options.Labels         = {'N','NE','E','SE','S','SW','W','NW'};
if tilenum == 1
    Options.LegendType     = 2;
else
    Options.LegendType     = 0;
end
Options.LegendPosition = 'westoutside';
Options.FreqLabelAngle = 120;
Options.MaxFrequency   = 40;
Options.nFreq          = 4;
Options.nDirections    = 24;
Options.min_radius     = .1;
Options.vWinds         = linspace(0, .5, 6);
% Options.nSpeeds        = 6;

Options.TextFontname   = 'Arial';
Options.TitleString    = [];
Options.LabLegend      = 'U (m s^{-1})';
Options.LegendVariable = 'u';
Options.FigColor       = 'w';
Options.CMap           = mycolormap;
Options.Gap            = .2;
Options.EdgeColor      = 'k';
Options.EdgeWidth      = 2;

% Options.FrequencyFontColor = 'r';
% Options.FrequencyFontWeight = 'bold';
% Options.FrequencyFontName = 'Comic Sans MS';
Options.FrequencyFontSize = fontsize*.6;
Options.FrequencyFontAngle = 'italic';

% Options.AxesFontColor = 'b';
Options.AxesFontWeight = 'bold';
% Options.AxesFontName = 'Rockwell Extra Bold';
Options.AxesFontSize = fontsize*.8;
Options.AxesFontAngle = 'normal';

% Options.TitleColor = [1 0.7 0.7];
Options.TitleFontSize = fontsize;
Options.TitleFontWeight = 'normal';
% Options.TitleFontName = 'Jokerman';

% Options.LegendColor = [1 0.4 0];
Options.LegendFontSize = fontsize;
Options.LegendFontWeight = 'normal';
% Options.LegendFontName = 'Calibri';
Options.LegendFontAngle = 'italic';

% Options.LegendBarColor = [0 0.8 0];
Options.LegendBarFontSize = fontsize;
Options.LegendBarFontWeight = 'demi';
Options.LegendBarFontAngle = 'normal';
% Options.LegendBarFontName = 'Gill Sans Ultra Bold Condensed';

end