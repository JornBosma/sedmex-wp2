function [theta_c, theta_w, theta_cw] = compute_Shields(tau_c, tau_w, tau_cw, rho_s, rho_w, D, g)
    
    % INPUT
    %     tau_c     current-induced bed shear stress [Pa]
    %     tau_w     wave-induced bed shear stress [Pa]
    %     tau_cw    total bed shear stress [Pa]
    %     rho_s     sediment density [kg/m^3]
    %     rho_w     water density [kg/m^3]
    %     D         grain diameter [m]
    %     g         gravitational acceleration [m/s^2]
    
    % OUTPUT
    %     theta_c   current-induced Shields number [-]
    %     theta_w   wave-induced Shields number [-]
    %     theta_cw  total Shields number [-]
    
    % non-dimensional Shields numbers
    theta_c = tau_c ./ ((rho_s - rho_w) * g * D);
    theta_w = tau_w ./ ((rho_s - rho_w) * g * D);
    theta_cw = tau_cw ./ ((rho_s - rho_w) * g * D);

end
