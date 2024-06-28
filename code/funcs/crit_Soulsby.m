function [theta_cr, tau_cr, Dstar] = crit_Soulsby(D, D50, rho_s, rho_w, g)
    
    % INPUT
    %     D          representative grain size [m]
    %     D50        median or mean (Dm) grain size [m]
    %     rho_s      sediment density [kg/m^3]
    %     rho_w      water density [kg/m^3]
    %     g          gravitational acceleration [m/s^2]

    % OUTPUT
    %     theta_cr   critical Shields number [-]
    %     tau_cr     critical bed shear stress [Pa]
    %     Dstar      dimensionless grain size [-]
    
    % constants
    nu = 1.3e-6;      % kinematic viscosity of water [m^2/s] (psu = 35 ppt, T = 12°C)
    
    s = rho_s ./ rho_w;  % relative density
    Dstar = D50 * ((s-1) * g / nu^2)^(1/3);  % Bonneville nondimensional grain size
    
    % % Brownlie (1981), also in 'The Civil Engineering Handbook' from Chen (1995)
    % Reg = sqrt(9.81 * (s - 1) * D.^3) / nu;      % grain Reynolds number
    % Y = Reg.^(-0.6);
    % theta_cr = 0.22 * Y + 0.06 * 10.^(-7.7 .* Y); % critical Shields number
    
    % % Soulsby & Whitehouse (1997)
    % theta_cr = 0.24 ./ Dstar + 0.055 * (1-exp(-0.02 * Dstar));
    
    % Soulsby (1997), no hiding-exposure
    theta_cr = 0.30 ./ (1+1.2 .* Dstar) + 0.055 .* (1-exp(-0.02 .* Dstar));
    
    % % Hanson & Camenen (2007)
    % theta_cr = 0.08 * (1-exp(-15./Dstar - 0.02.*Dstar));
    
    % dimensionalise threshold
    tau_cr = theta_cr .* (rho_s - rho_w) .* g .* D;

end
