function [theta_cr, tau_cr, Dstar] = crit_Egiazaroff(D, D50, rho_s, rho_w, g)
    % Function to compute the critical bed shear stress using the hiding function of Egiazaroff
    
    % Inputs:
    %   D       - diameter of the sediment fraction (m)
    %   D50     - median or mean (Dm) diameter of the sediment mixture (m)
    %   rho_s   - density of sediment particles (kg/m^3)
    %   rho_w   - density of water (kg/m^3)
    %   g       - acceleration due to gravity (m/s^2)
    
    % Output:
    %   tau_cr   - critical bed shear stress (Pa)
    %   theta_cr - critical Shields number (-)
    %   Dstar    - dimensionless grain size [-]

    % constants
    nu = 1.3e-6;      % kinematic viscosity of water [m^2/s] (psu = 35 ppt, T = 12°C)
    
    s = rho_s ./ rho_w;  % relative density
    Dstar = D50 .* ((s-1) .* g ./ nu.^2).^(1/3);  % Bonneville nondimensional grain size
    
    % Soulsby (1997), no hiding-exposure
    theta_cr_noHE = 0.30 ./ (1+1.2 .* Dstar) + 0.055 .* (1-exp(-0.02 .* Dstar));

    % Hiding function
    xi = (log10(19) ./ log10(19 .* (D./D50))).^2;

    % Non-dimensional critical shear stress (Shields parameter)
    theta_cr = xi .* theta_cr_noHE;

    % Critical bed shear stress
    tau_cr = theta_cr .* (rho_s - rho_w) .* g .* D;

end
