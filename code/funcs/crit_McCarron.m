function [theta_cr, tau_cr, Dstar] = crit_McCarron(D, D50, fg, rho_s, rho_w, g)

% INPUT
%     D          grain-size fraction [m]
%     D50        median or mean (Dm) grain size [m]
%     fg         gravel fraction (2-64 mm) [-]
%     rho_s      sediment density [kg/m^3]
%     rho_w      water density [kg/m^3]
%     g          gravitational acceleration [m/s^2]

% OUTPUT
%     theta_cr   critical Shields number incl. HE [-]
%     tau_cr     critical bed shear stress [Pa]
%     Dstar      dimensionless grain size [-]

% constants
nu = 1.3e-6;      % kinematic viscosity of water [m^2/s] (psu = 35 ppt, T = 12°C)
n = 1;            % sediment classes have equal density

s = rho_s ./ rho_w;  % relative density
Dstar = D50 .* ((s-1) .* g ./ nu.^2).^(1/3);  % Bonneville nondimensional grain size

% Soulsby (1997), no hiding-exposure
theta_cr_noHE = 0.30 ./ (1+1.2 .* Dstar) + 0.055 .* (1-exp(-0.02 .* Dstar));

% McCarron et al. (2019)
y_s = 0.68;  % empirically derived
y_g = 0.86;  % empirically derived

y_s = y_s + (y_g - y_s) .* fg.^1.73;
xi = n .* (D ./ D50).^(-y_s);
theta_cr = xi .* theta_cr_noHE;

% dimensionalise threshold
tau_cr = theta_cr .* (rho_s-rho_w) .* g .* D;

end
