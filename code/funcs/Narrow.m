function [] = Narrow(fontsize)

% Define starting point and arrow length
start_point = [117290, 560606];
arrow_length = 200;

% Calculate endpoint of the arrow
% end_point = start_point + [0, arrow_length];

% Plot arrow
quiver(start_point(1)-5, start_point(2), 0, arrow_length, 'r',...
    'LineWidth',2, 'MaxHeadSize',3, 'HandleVisibility','off');
quiver(start_point(1), start_point(2), 0, arrow_length, 'k',...
    'LineWidth',2, 'MaxHeadSize',3, 'HandleVisibility','off');

% Add the letter 'N' at the tip of the arrow
text(start_point(1)+80, start_point(2) + arrow_length, 'N',...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom',...
    'FontSize',fontsize, 'FontWeight','bold', 'Color','r');
text(start_point(1)+85, start_point(2) + arrow_length, 'N',...
    'HorizontalAlignment','center', 'VerticalAlignment','bottom',...
    'FontSize',fontsize, 'FontWeight','bold', 'Color','k');

end