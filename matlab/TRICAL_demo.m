% TRICAL_DEMO  Demonstrate the TRICAL magnetometer calibration filter.
%
%   Simulates a tri-axial magnetometer with known hard-iron bias and
%   soft-iron scale distortion, runs the TRICAL UKF to recover those
%   calibration parameters, and plots convergence.
%
%   Two use-cases are shown:
%     CASE 1 – Attitude-dependent:  the true field vector B_true is supplied
%              as the reference field (e.g. from a WMM model + attitude
%              estimate).  Recommended for fastest convergence.
%
%     CASE 2 – Attitude-independent:  the raw measurement itself is used as
%              the reference direction.  Requires field_norm ≈ 1 (normalize
%              measurements first) to avoid near-zero innovations at startup.
%
%   See also: TRICAL_INIT, TRICAL_ESTIMATE_UPDATE, TRICAL_CALIBRATE

clear; close all; clc;

%% -----------------------------------------------------------------------
%  Simulation parameters
% -----------------------------------------------------------------------
rng(0);
N = 2000;                       % number of sensor readings

true_field_norm = 50.0;         % e.g. 50 µT Earth-strength field

true_bias = [5.0; -3.0; 2.0];  % hard-iron bias

% Soft-iron / scale-error matrix D  (calibration model: cal = (I+D)*(raw-b))
true_D = [0.10  0.05  0.02;
          0.05 -0.08  0.03;
          0.02  0.03  0.06];

sensor_noise_std = 0.5;

%% -----------------------------------------------------------------------
%  Generate synthetic magnetometer data
% -----------------------------------------------------------------------
raw_meas  = zeros(3, N);
B_true_all = zeros(3, N);

for k = 1:N
    dir  = randn(3, 1);
    dir  = dir / norm(dir);
    B_true = true_field_norm * dir;

    % Distorted sensor output:  raw = (I+D)*B_true + bias + noise
    raw_meas(:, k)   = (eye(3) + true_D) * B_true + true_bias ...
                       + sensor_noise_std * randn(3, 1);
    B_true_all(:, k) = B_true;
end

%% =======================================================================
%  CASE 1 – Attitude-dependent calibration
%  Reference field = B_true  (known from WMM + attitude in a real system)
%  z = sqrt(|dot(cal, B_true)|)  →  equals field_norm when cal == B_true
% =======================================================================
inst1 = TRICAL_init(true_field_norm, sensor_noise_std);

bias_hist1  = zeros(3, N);
Ddiag_hist1 = zeros(3, N);
mag_hist1   = zeros(1, N);

for k = 1:N
    inst1 = TRICAL_estimate_update(inst1, raw_meas(:, k), B_true_all(:, k));

    bias_hist1(:, k)   = inst1.state(1:3);
    Ddiag_hist1(:, k)  = [inst1.state(4); inst1.state(8); inst1.state(12)];
    mag_hist1(k)       = norm(TRICAL_calibrate(inst1, raw_meas(:, k)));
end

%% =======================================================================
%  CASE 2 – Attitude-independent calibration
%  Normalize measurements to field_norm=1 so that at the initial state
%  (cal == raw_norm) the magnitude  ≈ 1 == field_norm  and innovations
%  are driven by the ellipsoidal spread, not the DC magnitude.
% =======================================================================
raw_norm = raw_meas / true_field_norm;   % scale to unit sphere

inst2 = TRICAL_init(1.0, sensor_noise_std / true_field_norm);

bias_hist2  = zeros(3, N);
mag_hist2   = zeros(1, N);

for k = 1:N
    % Pass the same (normalized) vector as both measurement and reference
    inst2 = TRICAL_estimate_update(inst2, raw_norm(:, k), raw_norm(:, k));

    bias_hist2(:, k)  = inst2.state(1:3) * true_field_norm;  % rescale for display
    cal_norm          = TRICAL_calibrate(inst2, raw_norm(:, k));
    mag_hist2(k)      = norm(cal_norm) * true_field_norm;
end

%% -----------------------------------------------------------------------
%  Report results
% -----------------------------------------------------------------------
fprintf('=== CASE 1 – Attitude-dependent (%d measurements) ===\n', N);
fprintf('  Bias estimate : [%7.4f  %7.4f  %7.4f]\n', inst1.state(1:3));
fprintf('  Bias truth    : [%7.4f  %7.4f  %7.4f]\n\n', true_bias);
D1 = reshape(inst1.state(4:12), 3, 3)';
fprintf('  D estimate:\n'); disp(D1);
fprintf('  D truth:\n');    disp(true_D);
fprintf('  Mean |cal| last 200: %.4f  (target %.4f)\n\n', ...
    mean(mag_hist1(end-199:end)), true_field_norm);

fprintf('=== CASE 2 – Attitude-independent (%d measurements) ===\n', N);
fprintf('  Bias estimate (rescaled): [%7.4f  %7.4f  %7.4f]\n', bias_hist2(:,end)');
fprintf('  Bias truth               : [%7.4f  %7.4f  %7.4f]\n\n', true_bias);
fprintf('  Mean |cal| last 200: %.4f  (target %.4f)\n\n', ...
    mean(mag_hist2(end-199:end)), true_field_norm);

%% -----------------------------------------------------------------------
%  Plots
% -----------------------------------------------------------------------
t = 1:N;

figure('Name','TRICAL – Case 1: Attitude-dependent','NumberTitle','off');

subplot(3,1,1);
plot(t, bias_hist1(1,:),'r', t, bias_hist1(2,:),'g', t, bias_hist1(3,:),'b');
yline(true_bias(1),'r--'); yline(true_bias(2),'g--'); yline(true_bias(3),'b--');
xlabel('Measurement index'); ylabel('Bias estimate');
title('Case 1 – Hard-iron bias convergence'); legend('bx','by','bz'); grid on;

subplot(3,1,2);
plot(t, Ddiag_hist1(1,:),'r', t, Ddiag_hist1(2,:),'g', t, Ddiag_hist1(3,:),'b');
yline(true_D(1,1),'r--'); yline(true_D(2,2),'g--'); yline(true_D(3,3),'b--');
xlabel('Measurement index'); ylabel('D diag estimate');
title('Case 1 – Soft-iron D (diagonal)'); legend('D11','D22','D33'); grid on;

subplot(3,1,3);
plot(t, mag_hist1,'k');
yline(true_field_norm,'r--','Target');
xlabel('Measurement index'); ylabel('||calibrated||');
title('Case 1 – Calibrated magnitude'); grid on;

figure('Name','TRICAL – Case 2: Attitude-independent','NumberTitle','off');

subplot(2,1,1);
plot(t, bias_hist2(1,:),'r', t, bias_hist2(2,:),'g', t, bias_hist2(3,:),'b');
yline(true_bias(1),'r--'); yline(true_bias(2),'g--'); yline(true_bias(3),'b--');
xlabel('Measurement index'); ylabel('Bias (rescaled)');
title('Case 2 – Hard-iron bias convergence (attitude-independent)');
legend('bx','by','bz'); grid on;

subplot(2,1,2);
plot(t, mag_hist2,'k');
yline(true_field_norm,'r--','Target');
xlabel('Measurement index'); ylabel('||calibrated||');
title('Case 2 – Calibrated magnitude'); grid on;
