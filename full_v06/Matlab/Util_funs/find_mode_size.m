function [sig_x2, sig_y2] = find_mode_size(Eint_2D, xgrid, ygrid)
%% Function for Finding the Second Moment Mode(Spot)-Size


% Determine Grid
dx = xgrid(2)-xgrid(1);
dy = ygrid(2)-ygrid(1);

% Integrate Intensity
norm = sum(Eint_2D, [1 2]) *dx*dy;

x2 = reshape(xgrid.^2, [], 1,1);
y2 = reshape(ygrid.^2, 1, [],1);

sig_x2 = sum(Eint_2D.*x2,[1 2] ) *dx*dy./norm;
sig_y2 = sum(Eint_2D.*y2, [1 2]) *dx*dy./norm;

end