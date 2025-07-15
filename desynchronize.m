function [ar,ai]= desynchronize(arn1,ain1,time_axs,my_time,xlngth,ylngth, LOSS_C)
%% Function to calculate the Radiation pulse in time and frequency domain

 for r=1:1:xlngth
        for c=1:1:ylngth
            testr=squeeze(squeeze(arn1(r,c,:)));
            arnn1=interp1(time_axs,testr,my_time, 'linear', 'extrap');
        
            testi=squeeze(squeeze(ain1(r,c,:)));
            ainn1=interp1(time_axs,testi,my_time, 'linear', 'extrap');
            ar(r,c,:)=arnn1*LOSS_C^0.5;
            ai(r,c,:)=ainn1*LOSS_C^0.5;
            % if (r<20) && (c<20)
            %     figure(4)
            %     title([' ',num2str(r),' ',num2str(c)])
            %     plot(squeeze(squeeze(ar1(r,c,:))))
            %     hold on
            %     plot(squeeze(squeeze(arn1(r,c,:))))
            %     %pause(1)
            %     hold off
            % end
        end
    end
end