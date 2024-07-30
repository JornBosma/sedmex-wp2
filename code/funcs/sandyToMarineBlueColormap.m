function cmap = sandyToMarineBlueColormap(n, reverse)
    % n is the number of levels in the colormap.
    % reverse is a boolean that controls whether the colormap is reversed.
    if nargin < 1
       n = 256; % Default number of levels
    end
    if nargin < 2
       reverse = false; % Default is not reversed
    end

    % Define RGB for darker sandy tone, white, and darker marine blue
    darkSandy = [0.68, 0.62, 0.45]; % Darker sandy tone
    white = [1, 1, 1];
    darkMarineBlue = [0.00, 0.45, 0.63]; % Darker marine blue

    % Create two gradients: one from dark sandy to white, and one from white to dark marine blue
    half_n = round(n / 2);
    gradient1 = [linspace(darkSandy(1), white(1), half_n)', linspace(darkSandy(2), white(2), half_n)', linspace(darkSandy(3), white(3), half_n)'];
    gradient2 = [linspace(white(1), darkMarineBlue(1), half_n)', linspace(white(2), darkMarineBlue(2), half_n)', linspace(white(3), darkMarineBlue(3), half_n)'];

    % Combine the two gradients
    cmap = [gradient1; gradient2];

    % If n is odd, repeat the middle row (white) to make the size of the colormap equal to n
    if mod(n, 2) == 1
        cmap = [cmap(1:half_n, :); white; cmap(half_n+1:end, :)];
    end

    % Reverse the colormap if required
    if reverse
        cmap = flipud(cmap);
    end
end
