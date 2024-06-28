function plot_north_arrow()
    % Define length of arrow
    arrow_length = 1;
    
    % Plot North arrow
    figure;
    hold on;
    
    % Plot North line
    plot([0, 0], [0, arrow_length], 'k', 'LineWidth', 2);
    
    % Plot arrowhead
    arrowhead_size = 0.1;
    arrow_x = [0, -arrowhead_size, 0, arrowhead_size, 0];
    arrow_y = [arrow_length, arrow_length*(1 - arrowhead_size), arrow_length - 0.1, arrow_length*(1 - arrowhead_size), arrow_length];
    plot(arrow_x, arrow_y, 'k', 'LineWidth', 2);
    
    % Plot text for North
    text(0, arrow_length, 'N', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    
    % Plot arrow pointing 45 degrees
    arrow_angle = deg2rad(45);
    arrow_end_x = arrow_length * sin(arrow_angle);
    arrow_end_y = arrow_length * cos(arrow_angle);
    plot([0, arrow_end_x], [0, arrow_end_y], 'r', 'LineWidth', 2);
    
    % Plot text for 45 degrees
    text(arrow_end_x, arrow_end_y, '45°', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
    
    % Set aspect ratio to equal
    axis equal;
    
    hold off;
end
