function [tau_c, tau_w, tau_cw] = compute_BSS_linearTheory(u_c, H, k, T, h, rho_w, phi_c, phi_w, D50, g)
    % Largely based on Zhu et al. (2016)

    % INPUT
    %     u_c       depth-averaged flow velocity [m/s]
    %     H         wave height [m]
    %     k         wave number [m^-1]
    %     T         wave period [s]
    %     h         water depth [m]
    %     rho_w     water density [kg/m^3]
    %     phi_c     current direction [°N]
    %     phi_w     wave propagation direction [°N]
    %     D50       median grain diameter [m]
    %     g         gravitational acceleration [m/s^2]
    
    % OUTPUT
    %     tau_c     current-induced bed shear stress [Pa]
    %     tau_w     wave-induced bed shear stress [Pa]
    %     tau_cw    total bed shear stress [Pa]

    k_s = 2.5*D50;                              % Nikuradse roughness (grain-related)
    C = 18 * log10(12 * h / k_s);               % Chézy coefficient
    tau_c = rho_w * g * u_c.^2 ./ C.^2;         % current-induced bed-shear stress
    
    % Linear wave theory
    Udelta = (pi * H) ./ (T .* sinh(k .* h));   % near-bed (peak) orbital velocity
    Adelta = H ./ (2 .* sinh(k .* h));          % near-bed (peak) orbital excursion
    r = Adelta ./ k_s;                          % relative roughness

    fw = exp(5.213 .* r.^(-0.194) - 5.977);     % wave friction factor
    tau_w = (1/2) .* rho_w .* fw .* Udelta.^2;  % wave-induced bed-shear stress
    
    phi_cw = wrapTo360(phi_w - phi_c);  % angle between current and wave propagation direction [°]

    % combined current- and wave-induced (total) bed-shear stress (Grant & Madsen, 1979)
    tau_cw = sqrt((tau_w + tau_c .* abs(cosd(phi_cw))).^2 + (tau_c .* abs(sind(phi_cw))).^2);

end
