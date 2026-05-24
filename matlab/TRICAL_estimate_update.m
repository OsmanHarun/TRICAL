function instance = TRICAL_estimate_update(instance, measurement, field)
% TRICAL_ESTIMATE_UPDATE  Run one UKF update step for TRICAL calibration.
%
%   instance = TRICAL_ESTIMATE_UPDATE(instance, measurement, field)
%
%   INPUTS:
%     instance    - Struct from TRICAL_INIT (or a previous call here)
%     measurement - 3x1 raw sensor reading in the same units as field_norm
%     field       - 3x1 reference field direction vector (need not be a unit
%                   vector; only its direction matters for the scalar
%                   measurement projection). For attitude-independent mode,
%                   pass the same raw measurement vector as field.
%
%   OUTPUT:
%     instance    - Updated struct with new .state, .P, .measurement_count
%
%   ALGORITHM OVERVIEW:
%     Implements the Scaled Unscented Kalman Filter from:
%       "Real-Time Attitude-Independent Three-Axis Magnetometer Calibration"
%       Alonso & Shuster, J. Astronaut. Sci., 2002
%       http://www.acsu.buffalo.edu/~johnc/mag_cal05.pdf
%
%     UKF parameters (alpha^2=1, beta=0, kappa=1  =>  lambda=1):
%       n               = 12 (state dimension)
%       dim_plus_lambda = 13
%       WM0 = WC0       = 1/13
%       WMI = WCI       = 1/26
%       Number of sigma points = 2*12 + 1 = 25
%
%     The scalar measurement model reduces each sigma-point state to a
%     single number:
%       z = sqrt(|dot(calibrated_B, field)|)
%     which is compared against field_norm (the expected value when fully
%     calibrated).
%
%   See also: TRICAL_INIT, TRICAL_CALIBRATE

% -------------------------------------------------------------------------
% UKF scaled parameters (fixed constants matching the C implementation)
% -------------------------------------------------------------------------
n               = 12;
dim_plus_lambda = 13.0;   % alpha^2 * (n + kappa) = 1 * (12 + 1)
WM0             = 1.0 / 13.0;
WC0             = 1.0 / 13.0;  % WM0 + (1 - alpha^2 + beta) = 1/13 + 0
WMI             = 1.0 / 26.0;  % 1 / (2 * dim_plus_lambda)
WCI             = 1.0 / 26.0;

state      = instance.state(:);   % ensure 12x1 column vector
P          = instance.P;
field_norm = instance.field_norm;
meas_noise = instance.measurement_noise;

measurement = measurement(:);
field       = field(:);

% -------------------------------------------------------------------------
% Step 1 – Cholesky decomposition of the scaled covariance
%   Find lower-triangular L  s.t.  L * L' = P * dim_plus_lambda
% -------------------------------------------------------------------------
L = chol(P * dim_plus_lambda, 'lower');

% -------------------------------------------------------------------------
% Step 2 – Generate 25 sigma points and compute scalar measurement estimates
%   z(1)        : central sigma point (state itself)
%   z(2..13)    : state + i-th column of L  (positive perturbations)
%   z(14..25)   : state - i-th column of L  (negative perturbations)
% -------------------------------------------------------------------------
z = zeros(2 * n + 1, 1);
z(1) = measure_reduce(state, measurement, field);

z_pair_sum = 0.0;
for i = 1:n
    z_pos        = measure_reduce(state + L(:, i), measurement, field);
    z_neg        = measure_reduce(state - L(:, i), measurement, field);
    z(i + 1)     = z_pos;
    z(i + 1 + n) = z_neg;
    z_pair_sum   = z_pair_sum + z_pos + z_neg;
end

% -------------------------------------------------------------------------
% Step 3 – Weighted mean of measurement estimates
%   z_mean = WM0 * z(1) + WMI * sum(z(2:25))
% -------------------------------------------------------------------------
z_mean = WM0 * z(1) + WMI * z_pair_sum;

% -------------------------------------------------------------------------
% Step 4 – Convert to deviations; compute scalar measurement covariance Pzz
%
%   NOTE: The original C implementation sums the squared deviations WITHOUT
%   per-point weights (i.e. unweighted sum).  This matches filter.c verbatim.
% -------------------------------------------------------------------------
dz  = z - z_mean;
Pzz = sum(dz .* dz) + meas_noise^2;

% -------------------------------------------------------------------------
% Step 5 – Cross-correlation vector Pxz (12x1)
%
%   Inner accumulation (unscaled), then apply WCI / WC0 weights.
%   Regenerates sigma-point positions on-the-fly (no extra storage needed).
%
%   Mathematical note: due to the zero-mean property of deviations,
%   the extra "state * sum(dz_pairs)" term vanishes and the result equals
%   the standard UKF cross-correlation formula:
%       Pxz = WCI * sum_i (dz_{+i} - dz_{-i}) * L(:,i)
% -------------------------------------------------------------------------
Pxz = zeros(n, 1);
for i = 1:n
    Pxz = Pxz ...
        + dz(i + 1)     * (state + L(:, i)) ...
        + dz(i + 1 + n) * (state - L(:, i));
end
Pxz = WCI * Pxz + WC0 * dz(1) * state;

% -------------------------------------------------------------------------
% Step 6 – Innovation: expected field norm vs estimated measurement mean
% -------------------------------------------------------------------------
innovation = field_norm - z_mean;

% -------------------------------------------------------------------------
% Step 7 – State update
%   Kalman gain K = Pxz / Pzz  (scalar denominator – no matrix inverse)
% -------------------------------------------------------------------------
K = Pxz / Pzz;
instance.state = state + K * innovation;

% -------------------------------------------------------------------------
% Step 8 – Covariance update
%   Standard form: P -= K * Pzz * K'
%   Since K = Pxz/Pzz:  K*Pzz*K' = (Pxz/Pzz)*Pzz*(Pxz/Pzz)' = Pxz*Pxz'/Pzz
%
%   The original C code drops the /Pzz (see filter.c comment), which is a
%   latent bug: it only works when field_norm ≈ 1 so Pzz ≈ O(1).  For any
%   other scale Pxz*Pxz' >> P and P loses positive-definiteness immediately.
% -------------------------------------------------------------------------
instance.P = P - Pxz * Pxz' / Pzz;

instance.measurement_count = instance.measurement_count + 1;
end

% =========================================================================
%  Local helper functions
% =========================================================================

function z = measure_reduce(state, measurement, field)
% Scalar measurement: sqrt(|dot(calibrated_B, field)|)
    cal = measure_calibrate(state, measurement);
    z   = sqrt(abs(dot(cal, field)));
end

function cal = measure_calibrate(state, measurement)
% Apply hard-iron bias and soft-iron scale correction.
%   cal = (I3 + D) * (measurement - bias)
%
%   state(1:3)  = bias b
%   state(4:12) = D matrix elements stored row-major in C convention.
%   reshape(...,3,3)' converts from C row-major to MATLAB column-major.
    b   = state(1:3);
    D   = reshape(state(4:12), 3, 3)';
    v   = measurement - b;
    cal = (eye(3) + D) * v;
end
