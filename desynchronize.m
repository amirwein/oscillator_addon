function [ar,ai]= desynchronize(arn1,ain1,time_axs,my_time,xlngth,ylngth, PWR_C)
%% Function to implement desynchronization 
% walk through the matrix
% use linear interpolation to implement the effect
% Multiply the field by the power retained in the cavity

 for r=1:1:xlngth
        for c=1:1:ylngth
            testr=squeeze(squeeze(arn1(r,c,:)));
            arnn1=interp1(time_axs,testr,my_time, 'linear', 'extrap');
        
            testi=squeeze(squeeze(ain1(r,c,:)));
            ainn1=interp1(time_axs,testi,my_time, 'linear', 'extrap');
            ar(r,c,:)=arnn1*PWR_C^0.5;
            ai(r,c,:)=ainn1*PWR_C^0.5;  
        end
    end
end
