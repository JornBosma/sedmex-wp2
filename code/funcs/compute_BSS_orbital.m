function [tau_c, tau_w, tau_cw] = compute_BSS_orbital(u_c, h, Urms, Trms, rho_w, phi_c, phi_w, D90, g)
    % Based on Kleinhans and Grasmeijer (2006), assuming hydraulic rough
    % conditions (Re* >= 11.63)

    % INPUT
    %     u_c       depth-averaged flow velocity [m/s]
    %     h         water depth [m]
    %     Urms      root-mean-square orbital velocity amplitude [m/s]
    %     T         wave period [s]
    %     rho_w     water density [kg/m^3]
    %     phi_c     current direction [°N]
    %     phi_w     wave propagation direction [°N]
    %     D90       roughness diameter [m]
    %     g         gravitational acceleration [m/s^2]
    
    % OUTPUT
    %     tau_c     current-induced bed-shear stress [Pa] (or [N/m^2])
    %     tau_w     wave-induced bed-shear stress [Pa]
    %     tau_cw    total bed-shear stress [Pa]
    
    k_s = D90;                               % Nikuradse roughness (grain-related)
    C = 18 * log10(12 * h / k_s);            % Chézy coefficient
    tau_c = rho_w * g * u_c.^2 ./ C.^2;      % current-induced bed-shear stress
    
    w = (2*pi) ./ Trms;                      % angular frequency [rad/s]
    Arms = Urms ./ w;                        % near-bed orbital velocity amplitude
    r = Arms / k_s;                          % relative roughness

    % Bosboom & Stive (2022) p. 198
    % if r < 1.59
    %     fw = 0.30;
    % else    
        fw = exp(5.213 * r.^(-0.194) - 5.977);   % wave friction factor
    % end

    tau_w = (1/2) .* rho_w .* fw .* Urms.^2; % wave-induced bed-shear stress
    
    phi_cw = wrapTo360(phi_w - phi_c);       % angle between current and wave propagation direction [°]

    % combined current- and wave-induced (total) bed-shear stress (Grant & Madsen, 1979)
    tau_cw = sqrt((tau_w + tau_c .* abs(cosd(phi_cw))).^2 + (tau_c .* abs(sind(phi_cw))).^2);

end
