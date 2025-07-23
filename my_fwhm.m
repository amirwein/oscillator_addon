function [width, peakPos] = my_fwhm(x, y)
    %% A function to return the optical pulse width and Peak position
    % x is the time axis
    % y is the pulse data
    
    % Normalize the signal
    y = y / max(y);

    % Find the peak position
    [~, peakIndex] = max(y);
    peakPos = x(peakIndex);
    
    % Find where the signal crosses half maximum
    halfMax = 0.5;
    indices = find(y >= halfMax);
    
    % get the indices from edges toward center
    i1 = find(y(1:indices(1)) < halfMax, 1, 'last');
    i2 = find(y(indices(end):end) < halfMax, 1, 'first')+ indices(end) - 1;

    % insert the indices
    x1 = x(i1);
    x2 = x(i2);
    
    % return pulse width
    width = abs(x2 - x1);
end
