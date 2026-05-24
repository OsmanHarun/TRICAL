function instance = TRICAL_init(field_norm, measurement_noise)
% TRICAL_INIT  Initialize a TRICAL magnetometer calibration filter instance.
%
%   instance = TRICAL_INIT()
%   instance = TRICAL_INIT(field_norm)
%   instance = TRICAL_INIT(field_norm, measurement_noise)
%
%   INPUTS:
%     field_norm        - Expected magnitude of the calibrated field
%                         (e.g. output of WMM for a magnetometer). Default: 1.0
%     measurement_noise - Sensor noise standard deviation. Default: 1e-6
%
%   OUTPUT:
%     instance - Struct with the following fields:
%       .field_norm        - Expected field magnitude
%       .measurement_noise - Measurement noise std dev
%       .state             - 12x1 state vector [bias(3); scale_D_rowmajor(9)]
%       .P                 - 12x12 state covariance matrix
%       .measurement_count - Number of measurements processed so far
%
%   STATE LAYOUT:
%     state(1:3)  = hard-iron bias vector [bx; by; bz]
%     state(4:12) = soft-iron scale matrix D, stored row-major:
%                   [D11 D12 D13  D21 D22 D23  D31 D32 D33]
%
%   CALIBRATION MODEL:
%     calibrated = (I3 + D) * (raw_measurement - bias)
%
%   Pass the returned struct to TRICAL_ESTIMATE_UPDATE for each sensor
%   reading, then optionally retrieve calibrated output with TRICAL_CALIBRATE.
%
%   See also: TRICAL_ESTIMATE_UPDATE, TRICAL_CALIBRATE, TRICAL_DEMO

if nargin < 1 || isempty(field_norm)
    field_norm = 1.0;
end
if nargin < 2 || isempty(measurement_noise)
    measurement_noise = 1e-6;
end

instance.field_norm        = field_norm;
instance.measurement_noise = measurement_noise;
instance.state             = zeros(12, 1);
instance.P                 = 0.01 * eye(12);
instance.measurement_count = 0;
end
