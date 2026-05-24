function calibrated = TRICAL_calibrate(instance, measurement)
% TRICAL_CALIBRATE  Apply the current calibration estimate to a raw measurement.
%
%   calibrated = TRICAL_CALIBRATE(instance, measurement)
%
%   INPUTS:
%     instance    - Struct from TRICAL_INIT / TRICAL_ESTIMATE_UPDATE
%     measurement - 3x1 raw sensor reading (same units as field_norm)
%
%   OUTPUT:
%     calibrated  - 3x1 calibrated reading
%
%   Implements:
%     calibrated = (I3 + D) * (measurement - bias)
%
%   where bias = instance.state(1:3) and D is the 3x3 scale-error matrix
%   stored row-major in instance.state(4:12).
%
%   WARNING: Always pass RAW measurements to TRICAL_ESTIMATE_UPDATE.
%   Do NOT feed calibrated output back into the filter.
%
%   See also: TRICAL_INIT, TRICAL_ESTIMATE_UPDATE

state = instance.state(:);
b     = state(1:3);
D     = reshape(state(4:12), 3, 3)';   % C row-major -> MATLAB column-major
v     = measurement(:) - b;
calibrated = (eye(3) + D) * v;
end
