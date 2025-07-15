function [width, peakPos] = my_fwhm(x, y)
    % Normalize the signal
    y = y / max(y);

    % Find the peak position
    [~, peakIndex] = max(y);
    peakPos = x(peakIndex);
    
    % Find where the signal crosses half maximum
    halfMax = 0.5;
    indices = find(y >= halfMax);
    
    % get the indices
    i1 = find(y(1:indices(1)) < halfMax, 1, 'last');
    i2 = find(y(indices(end):end) < halfMax, 1, 'first')+ indices(end) - 1;

    % insert the indices
    x1 = x(i1);%interp1(y(i1:i1+1), x(i1:i1+1), halfMax);
    x2 = x(i2);%interp1(y(i2-1:i2), x(i2-1:i2), halfMax);

    width = abs(x2 - x1);
end