function [waist_x, waist_y] = find_waistxy(E_in, N, x, y)
    % Calculate the intensity of the beam

    intensity = abs(E_in).^2;  % Intensity of the beam

% Find the maximum intensity and normalize
I_max = max(intensity(:));
intensity = intensity / I_max;

%% Take center row and column
center_idx = floor(N/2);
Ix = intensity(center_idx, :); % Horizontal slice (along x)
Iy = intensity(:, center_idx); % Vertical slice (along y)

%% Only use positive side (from center to end)
Ix_pos = Ix(center_idx:end);
Iy_pos = Iy(center_idx:end);
xx_pos = x(center_idx:end);
yy_pos = y(center_idx:end);

% Threshold for 1/e^2
threshold = exp(-2); % ≈ 0.135

% Find where intensity drops below 1/e^2
idx_x = find(Ix_pos <= threshold, 1);
idx_y = find(Iy_pos <= threshold, 1);

% Compute waist distances
waist_x = xx_pos(idx_x);
waist_y = yy_pos(idx_y);
end