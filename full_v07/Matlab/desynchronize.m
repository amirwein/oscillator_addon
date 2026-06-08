function [ar,ai]= desynchronize(arn1,ain1,time_axs,my_time,xlngth,ylngth)
%% Function to calculate the Radiation pulse in time and frequency domain

sz = size(arn1);
ar = zeros(sz);
ai = zeros(sz);

    for r=1:1:xlngth
        for c=1:1:ylngth
            testr = squeeze( arn1(r,c,:) );
            arnn1 = interp1(time_axs,testr,my_time, 'linear');
            TFr = isnan(arnn1);
            arnn1(TFr) = 0;

            testi = squeeze( ain1(r,c,:) );
            ainn1 = interp1(time_axs,testi,my_time, 'linear');
            TFi = isnan(ainn1);
            ainn1(TFi) = 0;
            
            ar(r,c,:) = arnn1;
            ai(r,c,:) = ainn1;

        end
    end
end