% TRICAL_DEMO  Demonstrate the TRICAL magnetometer calibration filter.
%
%   Simulates a magnetometer with known hard-iron bias and soft-iron scale
%   distortion, runs the TRICAL UKF to recover the calibration parameters,
%   and plots convergence of the bias and calibrated field magnitude.
%
%   See also: TRICAL_INIT, TRICAL_ESTIMATE_UPDATE, TRICAL_CALIBRATE

clear; close all; clc;

%% -----------------------------------------------------------------------
%  Simulation parameters
% -----------------------------------------------------------------------
rng(0);                        % reproducible random numbers
N            = 1000;           % number of sensor readings

true_field_norm = 50.0;        % e.g. 50 µT Earth-strength field

% Hard-iron (additive) bias
true_bias = [5.0; -3.0; 2.0];

% Soft-iron (multiplicative) scale-error matrix D
% Calibrated measurement = (I + D) * (raw - bias)
% so raw = (I+D)^{-1} * true_field + bias
true_D = [0.10  0.05  0.02;
          0.05 -0.08  0.03;
          0.02  0.03  0.06];

sensor_noise_std = 0.5;        % measurement noise standard deviation

%% -----------------------------------------------------------------------
%  Generate synthetic magnetometer data
%
%  At each time step we pick a random orientation (random unit vector),
%  scale it to the true field norm, apply the sensor distortion model,
%  and add Gaussian noise.
% -----------------------------------------------------------------------
raw_measurements = zeros(3, N);
true_field_dirs  = zeros(3, N);   % undistorted, for reference

for k = 1:N
    % Random orientation unit vector
    dir  = randn(3, 1);
    dir  = dir / norm(dir);

    % Undistorted field in body frame
    B_true = true_field_norm * dir;

    % Distorted raw sensor output: raw = (I + D) * B_true + bias + noise
    raw_measurements(:, k) = (eye(3) + true_D) * B_true + true_bias ...
                             + sensor_noise_std * randn(3, 1);
    true_field_dirs(:, k)  = B_true;
end

%% -----------------------------------------------------------------------
%  Run TRICAL filter
% -----------------------------------------------------------------------
instance = TRICAL_init(true_field_norm, sensor_noise_std);

% Logging
bias_history  = zeros(3, N);
D_diag_history = zeros(3, N);
mag_history   = zeros(1, N);

for k = 1:N
    % For attitude-independent calibration, pass the raw measurement as
    % both the observation and the reference field direction.
    % (See filter.c: "the same vector can be supplied for measurement and field")
    instance = TRICAL_estimate_update(instance, raw_measurements(:, k), ...
                                      raw_measurements(:, k));

    % Log state
    bias_history(:, k)   = instance.state(1:3);
    D_diag_history(:, k) = [instance.state(4); instance.state(8); instance.state(12)];

    % Apply calibration and record output magnitude
    cal = TRICAL_calibrate(instance, raw_measurements(:, k));
    mag_history(k) = norm(cal);
end

%% -----------------------------------------------------------------------
%  Report final estimates
% -----------------------------------------------------------------------
fprintf('--- TRICAL Demo Results (%d measurements) ---\n\n', N);

fprintf('Hard-iron bias estimate:  [%8.4f  %8.4f  %8.4f]\n', instance.state(1:3));
fprintf('Hard-iron bias truth:     [%8.4f  %8.4f  %8.4f]\n\n', true_bias);

D_est = reshape(instance.state(4:12), 3, 3)';
fprintf('Soft-iron D estimate:\n');
disp(D_est);
fprintf('Soft-iron D truth:\n');
disp(true_D);

fprintf('Mean calibrated magnitude (last 100 samples): %.4f  (target: %.4f)\n', ...
    mean(mag_history(end-99:end)), true_field_norm);

%% -----------------------------------------------------------------------
%  Plots
% -----------------------------------------------------------------------
t = 1:N;

figure('Name', 'TRICAL Calibration Convergence');

subplot(3, 1, 1);
plot(t, bias_history(1,:), 'r', t, bias_history(2,:), 'g', t, bias_history(3,:), 'b');
yline(true_bias(1), 'r--'); yline(true_bias(2), 'g--'); yline(true_bias(3), 'b--');
xlabel('Measurement index'); ylabel('Bias estimate');
title('Hard-iron bias convergence');
legend('bx','by','bz','bx true','by true','bz true','Location','best');
grid on;

subplot(3, 1, 2);
plot(t, D_diag_history(1,:), 'r', ...
     t, D_diag_history(2,:), 'g', ...
     t, D_diag_history(3,:), 'b');
yline(true_D(1,1), 'r--'); yline(true_D(2,2), 'g--'); yline(true_D(3,3), 'b--');
xlabel('Measurement index'); ylabel('D diagonal estimate');
title('Soft-iron scale D (diagonal elements)');
legend('D11','D22','D33','D11 true','D22 true','D33 true','Location','best');
grid on;

subplot(3, 1, 3);
plot(t, mag_history, 'k');
yline(true_field_norm, 'r--', 'Target');
xlabel('Measurement index'); ylabel('||calibrated||');
title('Calibrated measurement magnitude');
grid on;
