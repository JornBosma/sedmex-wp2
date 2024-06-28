function fig = figureRH()
    % Get screen sizes
    screenSizes = get(0, 'MonitorPositions');

    % Find the row index with the largest value in the third column
    [~, max_row_index] = max(screenSizes(:, 3));

    % Assign the row with the largest value in the third position to a new variable
    screenSize = screenSizes(max_row_index, :);
    
    % Calculate position for the new figure
    figureWidth = screenSize(3) / 2;
    figureHeight = screenSize(4);
    figurePosition = [figureWidth+screenSize(1), 0, figureWidth, figureHeight];
    
    % Create new figure with custom position
    fig = figure('Position', figurePosition);
end
