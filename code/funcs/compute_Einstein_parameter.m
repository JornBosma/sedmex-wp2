function phi_b = compute_Einstein_parameter(theta, theta_cr, alpha, beta)
    % INPUT
    %     theta     Shields number [-]
    %     theta_cr  critical Shields number [-]
    %     alpha     constant
    %     beta      constant
    
    % OUTPUT
    %     phi_b     nondimensional bedload predictor [-]

    % Ensure theta_cr is compatible with theta for element-wise operations
    if isscalar(theta_cr)
        theta_cr = theta_cr * ones(size(theta));
    end
    
    theta_diff = theta - theta_cr;

    % Compute phi_b using element-wise operations
    phi_b = alpha * (theta_diff).^beta;
    
    % Set phi_b to 0 where theta_diff is less than or equal to 0
    phi_b(theta_diff <= 0) = 0;
end
