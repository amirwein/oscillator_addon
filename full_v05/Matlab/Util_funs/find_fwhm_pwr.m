function [tau, peakPos, pwr_max] = find_fwhm_pwr(pwr_s, rad_s, rad_e, delt)
%% Function to find the Pulse Duration from its Profile

time_ax = (rad_e-delt/2:-delt:rad_s+delt/2);

xq = linspace(min(time_ax), max(time_ax), 10000);

yq = interp1(time_ax, pwr_s, xq, 'pchip');

[ymax, peakIndex] = max(yq);
peakPos = xq(peakIndex);
halfMax = ymax/2;

above = yq >= halfMax;
i1 = find(diff(above) == 1) + 1;
i2 = find(diff(above) == -1);

filt_l = min( abs(peakPos - xq(i1)) );
filt_r = min( abs(peakPos - xq(i2)) );

tau = filt_l + filt_r;
pwr_max = ymax;

end
