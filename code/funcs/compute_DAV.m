function [u_c, u_star] = compute_DAV(u_z, z, k_sc, h)
    % Estimates the depth-averaged velocity assuming the law of the wall up
    % to the water surface (Kleinhans and Grasmeijer, 2006)

    % INPUT
    %     u_z       measured flow velocity [m/s]
    %     z         height above the bed [m]
    %     k_sc      Nikuradse current-related roughness [m]
    %     h         water depth [m]
    
    % OUTPUT
    %     u_c       depth-averaged velocity [m/s]
    %     u_star    shear velocity [m/s]
    
    % Constants
    kappa = 0.41;  % von Karman constant

    % Ensure input arrays are the same size
    if ~isequal(size(u_z), size(z), size(h))
        error('All input arrays must be of the same size');
    end

    % Calculate shear velocity (u_star) element-wise
    u_star = u_z .* (kappa ./ log(z ./ k_sc));

    % Calculate depth-average velocity (U) at depth 0.368h element-wise
    u_c = u_star .* (1 ./ kappa) .* log((0.368 .* h) ./ k_sc);

end