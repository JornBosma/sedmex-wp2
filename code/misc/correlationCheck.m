% Assuming A and B are your time series vectors

% Step 1: Create logical masks where values are not NaN
valid_A = ~isnan(u_z);
valid_B = ~isnan(phi_w);

% Step 2: Combine the masks to find where both time series have non-NaN values
valid_indices = valid_A & valid_B;

% Step 3: Extract the values from A and B where both are valid
A = u_z(valid_indices);
B = phi_w(valid_indices);

% Optionally, you can also extract the corresponding time steps if you have them
% t_filtered = t(valid_indices); % Assuming you have a time vector 't'

%% 1. Correlation Analysis

% Assuming A is the flow velocity time series and B is the wave height time series
correlation_coefficient = corr(A, B);

fprintf('The correlation coefficient between A and B is: %.3f\n', correlation_coefficient);

%% 2. Linear Regression

% Perform linear regression: A = a*B + b
coefficients = polyfit(B, A, 1);  % 1 for linear fit
slope = coefficients(1);
intercept = coefficients(2);

% Display the result
fprintf('Linear regression model: A = %.3f*B + %.3f\n', slope, intercept);

% To visualize the fit:
A_predicted = polyval(coefficients, B);
figure;
plot(B, A, 'o', B, A_predicted, '-');
xlabel('Wave Height (B)');
ylabel('Flow Velocity (A)');
legend('Measured A', 'Fitted A');
title('Linear Regression of A vs B');

%% 3. Cross-Spectral Analysis

% Perform cross-spectral analysis
[cxy, f] = mscohere(B, A, [], [], [], 1);  % mscohere computes magnitude-squared coherence

% Plot the coherence to see the frequency dependence
figure;
plot(f, cxy);
xlabel('Frequency (Hz)');
ylabel('Magnitude-Squared Coherence');
title('Coherence between Wave Height and Flow Velocity');

%% 4. Wave Filtering using Fourier Transform

% Fourier transform of A
n = length(A);
A_fft = fft(A);
f = (0:n-1)*(1/n);  % Frequency range

% Apply filter based on the wave frequency range
% Assuming the wave frequency is known (e.g., between f_min and f_max)
f_min = 0.2; % Example lower bound of wave frequency in Hz
f_max = 2;  % Example upper bound of wave frequency in Hz
filter = (f >= f_min & f <= f_max);

% Filter out components not associated with waves
A_filtered = A_fft;
A_filtered(~filter) = 0;  % Zero out non-wave components

% Inverse Fourier transform to get the filtered velocity
A_wave_component = ifft(A_filtered);

% Plot the original and filtered signals
figure;
plot(A, 'b');
hold on;
plot(real(A_wave_component), 'r');
xlabel('Time');
ylabel('Flow Velocity (A)');
legend('Original A', 'Wave-Influenced A');
title('Original vs Wave-Influenced Flow Velocity');

% After trying these methods, you could compare the results and assess how
% much the waves (represented by B) influence the flow velocity (A). If
% you find a high correlation or coherence, it suggests that the waves
% strongly influence the velocity.