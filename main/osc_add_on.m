%% MAin file for FEL Oscillator add on for GINGER3d
clear all;
close all;

addpath("ginger","tst");
%% Some configurations
SAV=0
sav_dir='/Users/amirweinberg/ginger_int/superrad/';
nam="025ps480c";
%sav_dir1='/Users/amirweinberg/ginger_int/plslngt/';
TSTIN=1


%% General parameters
c= 299792458;                                                              % speed of light

%% UHM Beam parameters
energy  = 40;        						                               % electron energy [MeV]
eSpread = 0.5e-2;       					                               % fractional rms energy spread 
emitN   = 8e-6;           					                               % normalized transverse emittance [mm-mrad]
currentMax = 480;					                                       % peak current [Ampere]
betaX_twiss = 1.4;                                                         % twiss parameters
betaY_twiss = 0.24;
alphaX_twiss = 0.4714;
alphaY_twiss = 0.0;
bunch_sigma = 0.125e-12;                                                   % e-beam bunch rms Gaussian pulse width

%% UHM FEL parameters
unduPeriod = 0.023;     						                           % undulator period [meter]
unduK = 1.2/sqrt(2);    						    	                   % RMS undulator parameter, (K) 
NumPeriod = 47;                                                            % Number of undulator periods
unduL = unduPeriod*NumPeriod;  						                       % length of undulator [meter]
radWavelength = 3.2281e-06;                                                % seed wavelength? [meter]
fcar =c/radWavelength;                                                     % carrier frequency
slipp = radWavelength*NumPeriod                                            % Slippage length 

%% UHM Oscillator parameters
CavityLength= 2.0469;                                                      % Cavity length
MirrorRadiusCurv=1.3;                                                      % Mirrors Radius of Curvatures
LOSS_C = 0.93;                                                             % CAVITY Reflectivity    
npass=200;                                                                   % Number of passes in the cavity
dd=0.007;%s[0.007 0.025 0.05];                                             % desynchronization values                                                   % Desynchronization value         
if TSTIN
    npass=5;
end


%% GINGER3D simulation paramteres
ebunch_slices = 96;                                                        % Number of electron bunch slices
rad_slices = 192;                                                          % Number of Radiaion slices
intrct_slices = 15;                                                        % Totall # of Photon slices each elec beam slice will interact 
r_grid = 0.080e-3;                                                         % resolution in x,y
n_grid = 129;                                                              % Number of radial grid zones



inp_struc.energy      = energy;
inp_struc.eSpread     = eSpread;
inp_struc.emitx       = emitN;
inp_struc.emity       = emitN;
inp_struc.currentMax  = currentMax;
inp_struc.betaX       = betaX_twiss
inp_struc.betaY       = betaY_twiss;
inp_struc.alphaX      = alphaX_twiss;
inp_struc.alphaY      = alphaY_twiss;
inp_struc.RMSwidth    = bunch_sigma;

inp_struc.unduPeriod  = unduPeriod;
inp_struc.NumPeriod   = NumPeriod;
inp_struc.unduL       = unduL;
inp_struc.RMSunduK    = unduK;
inp_struc.radWvlength = radWavelength;

inp_struc.eslices     = ebunch_slices;
inp_struc.rslices     = rad_slices;
inp_struc.intslices   = intrct_slices;

inp_struc.rgrid       = r_grid;
inp_struc.ngrid       = n_grid;


myvid=0;                                                                   % Whether to Save a video file

for desyn=dd
desyn_shft=desyn*(47*radWavelength);                                       % Calculate cavity shift length in meters 
filnme="video_pulse"+num2str(desyn)+".mp4"
if (myvid)
v =  VideoWriter(filnme, 'MPEG-4');% Use 'MPEG-4' for .mp4 format
v.FrameRate = 10; % Set frames per second
open(v);
end

                                                                  
mycolors=parula(npass);
my_pls=zeros(1,npass);


%% First pass in undulator
if TSTIN
    [ar1, ai1,delt,i0,dx1,rad_s,rad_e,ebm_s,ebm_e] = input_tst(1,1,1,inp_struc);                %TEST withought GINGER 5 passes
else
    [ar1, ai1,delt,i0,dx1,rad_s,rad_e,ebm_s,ebm_e] = input_EE(1,1,1,inp_struc);                % Run GINGER3d for first pass 
end


xlngth=length(ar1(:,1,1));
ylngth=length(ar1(1,:,1));
tlngth=length(ar1(1,1,:));


[arn1,ain1]=inser_lens2(ar1,ai1,radWavelength,CavityLength,unduL, ...
    MirrorRadiusCurv/2,dx1,desyn_shft);                                    % implement cavity 

dshifts_s=desyn_shft/c;                                                    %Desynchronization value in time

time_axs=(rad_e-delt:-delt:rad_s);                                         % Generate time axis
f_axs=(-0:(length(time_axs)-1))/((time_axs(1)-time_axs(end)));             % Generate frequencey axis 
f_axs=f_axs-median(f_axs)+fcar;                                            % Shift frequency axis to carrier frequency

my_time=time_axs+dshifts_s;

[ar1,ai1]= desynchronize(arn1,ain1,time_axs,my_time, ...
    xlngth,ylngth, LOSS_C);                                                % implement desynchronization 



[power_z,power_s]=my_funct(ar1,ai1,i0,dx1);                                % Calculate radiation pulse in time and frequency domain 
my_pls(1)=max(power_s);                                                    % Save radiation pulse peak power

disp(['pulse1 pass=',num2str(1),' pulse peak power=',num2str(my_pls(1)),'[W]'])     % Display input power



%% Second pass and on in undulator
disp_data1=zeros(length(ar1(1,1,:)), 25);                                   %initialize Plotting paramters
disp_data_labels1=strings(1,25);
disp_data2=zeros(length(ar1(1,1,:)), 25);
disp_data_labels2=strings(1,25);
disp_data3=zeros(length(ar1(1,1,:)), 25);                                   %Plostting paramter
disp_data_labels3=strings(1,25);
disp_data4=zeros(length(ar1(1,1,:)), 25);
disp_data_labels4=strings(1,25);

width=zeros(1,npass);
peakPos=zeros(1,npass);

i=1;
k=1;
l=1;
m=1;
wst1=0;

wst2=0;

for j=2:npass           
    if TSTIN
        [arn, ain,delt,i0,dx,rad_s1,rad_e1,ebm_s,ebm_e] = input_tst(j,ar1,ai1,inp_struc);                    %TEST withought GINGER 5 passes
        pause(1);
    else
        [arn, ain,delt,i0,dx,rad_s1,rad_e1,ebm_s,ebm_e] = input_EE(j,ar1,ai1,inp_struc);                     % Run GINGER3d for j pass
    end 
    
    [arn1,ain1,wst1,wst2]=inser_lens2(arn,ain,radWavelength,CavityLength,unduL, ...
    MirrorRadiusCurv/2,dx,desyn_shft);                                     % implement cavity 

    Zray(wst1,wst2, CavityLength, unduL, desyn);                           % Return Z rayleigh  
    
    [ar1,ai1]= desynchronize(arn1,ain1,time_axs,my_time, ...
        xlngth,ylngth, LOSS_C);                                            % implement desynchronization 
    
    [power_z,power_s,angl]=my_funct(arn,ain,i0,dx);                        % Calculate radiation pulse 
    
    my_pls(j)=max(power_s);                                                % Save radiation pulse peak power
   
    fig23=figure(23);
    
    % Display radiation pulse 
    subplot(1,2,1);%hold on;
    yyaxis left
    plot(time_axs, power_s,'Color',mycolors(j,:),'LineWidth',1);
    title(['pass=',num2str(j)]);
    %legend(['desync=',num2str(desyn)]);
    xline(0-ebm_e+slipp/c,'--','e-beam head')
    xline(0+ebm_e+slipp/c,'--','e-beam tail')
%    yyaxis right;plot(s,unwrap(angle(Et))/pi,'Color',mycolors(1,:),'LineWidth',1.5,'LineStyle','--');
    xlabel('time [s]')
    ylabel('power (W)')
    yyaxis right
    plot(time_axs, angl)
    ylabel('Phase [rad]')
    %ylim([-5*pi pi])
    % ylabel('power (GW)')
    %savefig(fig23, sav_dir1+"pls"+num2str(desyn)+"pass"+num2str(j)+".fig")
    
    if (myvid)
    frame = getframe(gcf);  % Or getframe(gca) for axes only
    writeVideo(v, frame);
    end

    %fig25=figure(25);
    subplot(1,2,2);
    plot(f_axs,power_z,'.-','LineWidth',1.5);
    xlabel('f (Hz)')
    ylabel('a.u.')
    title(['pass=',num2str(j)]);
    drawnow
    disp(['pulse pass=', num2str(j),', pulse peak power=',num2str(my_pls(j)),'[W]'])
    %savefig(fig25, sav_dir1+"spct"+num2str(desyn)+"pass"+num2str(j)+".fig")

    [width(j), peakPos(j)] = my_fwhm(time_axs, power_s);
    
    fprintf('FWHM: %.3f ps\n', width(j)* 1e12);
    fprintf('Peak position: %.3f ps\n', peakPos(j)*1e12);

    
    if (j>=101) && (j<=125) 
        disp_data_labels1(i) = ["pulse #" + num2str(j)];
        disp_data1(:,i)=power_s;
        i=i+1;
    elseif (j>=126) && (j<=150)
        disp_data_labels2(k) = ["pulse #" + num2str(j)];
        disp_data2(:,k)=power_s;
        k=k+1;
    elseif (j>=151) && (j<=175) 
        disp_data_labels3(l) = ["pulse #" + num2str(j)];
        disp_data3(:,l)=power_s;
        l=l+1;
    elseif (j>=176) && (j<=200)
        disp_data_labels4(m) = ["pulse #" + num2str(j)];
        disp_data4(:,m)=power_s;
        m=m+1;
    end
end


%% Display stacked plots
fig3=figure(3)
stackedplot(time_axs, disp_data1(:,1:25), 'DisplayLabels', disp_data_labels1(1:25));
title('Stacked Plot of 25 radiation pulses 101-125');
xlabel('Time');
fig4=figure(4)
stackedplot(time_axs, disp_data2(:,1:25), 'DisplayLabels', disp_data_labels2(1:25));
title('Stacked Plot of 25 radiation pulses 126-150');
xlabel('Time');
fig8=figure(8)
stackedplot(time_axs, disp_data3(:,1:25), 'DisplayLabels', disp_data_labels3(1:25));
title('Stacked Plot of 25 radiation pulses 151-175');
xlabel('Time');
fig9=figure(9)
stackedplot(time_axs, disp_data4(:,1:25), 'DisplayLabels', disp_data_labels4(1:25));
title('Stacked Plot of 25 radiation pulses 176-200');
xlabel('Time');

%% Display radiation pulse vs Oscillator pass number
fig1= figure(1)
plot(my_pls)
title('Radiation pulse peak power vs Oscillator pass # emittance', num2str(emitN))
ylabel('W')
xlabel('pass number')


fig6=figure(6)
plot(width*1e12);
title('\Delta\tau - FWHM');
ylabel('pSec')
xlabel('pass number');

fig7=figure(7)
plot(peakPos*1e12);
title('Peak Position');
ylabel('pSec')
xlabel('pass number');


%% save figures
if (SAV)
    savefig(fig23, sav_dir+nam+"lstpls"+num2str(desyn)+".fig")
    savefig(fig3, sav_dir+nam+"101-125"+num2str(desyn)+".fig")
    savefig(fig4, sav_dir+nam+"126-150"+num2str(desyn)+".fig")   
    savefig(fig8, sav_dir+nam+"151-175"+num2str(desyn)+".fig")
    savefig(fig9, sav_dir+nam+"176-200"+num2str(desyn)+".fig")
    savefig(fig1, sav_dir+nam+"peakpass"+num2str(desyn)+".fig")
    savefig(fig6, sav_dir+nam+"FWHM"+num2str(desyn)+".fig")
    savefig(fig7, sav_dir+nam+"peakPos"+num2str(desyn)+".fig")
end

if (myvid)
close(v);
%MYPOWER(vr)=my_pls(end)
end

end

