function q = compute_transport_rate(phi, rho, rho_s, D, g)
    % INPUT
    %     phi       nondimensional transport predictor [-]
    %     rho       water density [kg/m^3]
    %     rho_s     sediment density [kg/m^3]
    %     D         representative grain size [m]
    %     g         gravitational acceleration [m/s^2]
    
    % OUTPUT
    %     q         sediment transport rate [m^3/m(width) /s] -> [m^2/s]
    %               (excluding pore space)

    % relative submerged density
    R = (rho_s - rho) / rho;

    % dimensional transport rate
    q = phi .* sqrt(R * g) * D^(1.5);

end
