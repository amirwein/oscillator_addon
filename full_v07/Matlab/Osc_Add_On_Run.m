function out = ...
    Osc_Add_On_Run(inp_struc,npass,cav_d,jitterRMS,d_ramp,EFF,path_struc)
%% Function to run FEL Oscillator add on for GINGER3d
% Replaces Osc_Add_On script, but has identical functionality
% Written for parallellable usage


if length(d_ramp) < npass
    d_ramp(end+1:npass) = 0;
end

%% General params
c = 299792458;
m_e = 0.511;                                     % electron rest mass [MeV]

%% UHM Oscillator parameters
CavityLength= 2.0469;                                       % Cavity length
MirrorRadiusCurv=1.3;                        % Mirrors Radius of Curvatures
PWR_C = 0.93;                                         % CAVITY Reflectivity   


%% Localize variables
radWavelength = inp_struc.radWvlength;
slipp = radWavelength* inp_struc.NumPeriod;
unduL = inp_struc.unduL;
rslices = inp_struc.rslices;
fcar = c/radWavelength;


%% Setup Axes/First Run
% Simulation output
my_pls=zeros(1,npass);
my_nrg=zeros(1,npass);
my_tau=zeros(1,npass);
pwr_prfl = zeros(npass, rslices);
pwr_spct = zeros(npass, rslices);

% Run first pass GINGER3d for params

[ar1, ai1,delt,i0,dx1,rad_s,rad_e,ebm_s,ebm_e] = input_EE(1,1,1,...
                                                inp_struc,path_struc,EFF);   

% Number of gridpoints along each dimension
xlngth=length(ar1(:,1,1));
ylngth=length(ar1(1,:,1));

% Time axes
time_axs = (rad_e-delt/2:-delt:rad_s+delt/2);

% Frq axes
f_ax = (-0:(length(time_axs)-1))/((time_axs(1)-time_axs(end)));
f_axs= f_ax-median(f_ax)+fcar;

% Cavity desync in meters, seconds
desyn_shft = cav_d * slipp;
dshift_s = desyn_shft/c;

% Time frame of field entering undulator
my_time = time_axs+dshift_s;


%% Finish 1st pass, Start later passes
% Implement Cavity
[arn1,ain1]=inser_lens2(ar1,ai1,radWavelength,CavityLength,unduL,...
                        MirrorRadiusCurv/2,dx1,desyn_shft,0);

% Shift radiation field time-frame (due to cavity desync)
[ars1, ais1] = desynchronize(arn1,ain1,time_axs,my_time,xlngth,ylngth);

% Cavity loss
arL = ars1 * PWR_C^0.5;
aiL = ais1 * PWR_C^0.5;

% Calculate radiation pulse exiting/entering undulator
[power_z, power_s, angl] = find_pwr(ar1, ai1,i0,dx1);
[power_zL,power_sL,anglL] = find_pwr(arL,aiL,i0,dx1);

% Save peak radiation power/energy
%my_pls(1) = max(power_s);                      
my_nrg(1) = sum(power_s) * delt;
[my_tau(1), peakPos, my_pls(1)] = find_fwhm_pwr(power_s, rad_s,rad_e,delt);

% Save power profile/spectra
pwr_prfl(1,:) = power_s;
pwr_spct(1,:) = power_z;

% Display input power - don't in parallel
%disp(['pulse pass=', num2str(1), ' pulse peak power=',...
%        num2str(my_pls(1)), '[W]'])


for j=2:npass
    % Get shift for current pass
    my_shft = jitterRMS*randn() + d_ramp(j)*slipp/c;
    new_time = time_axs + my_shft;

    % Shift radiation to current frame
    [ars,ais] = desynchronize(arL,aiL,time_axs,new_time,xlngth,ylngth);

    % Call Ginger-3D (centered on e-beam) for pass j
    [ar, ai,delt,i0,dx,rad_s,rad_e,ebm_s,ebm_e]=input_EE(j, ars, ais, ...
                                                 inp_struc,path_struc,EFF);

    % Reverse timing-offset shift
    [arns, ains] = desynchronize(ar,ai, new_time, time_axs,xlngth,ylngth);

    % Implement Cavity
    [arn,ain] = inser_lens2(arns,ains,radWavelength,CavityLength,unduL,...
                            MirrorRadiusCurv/2,dx,desyn_shft,0);

    % Shift field to undulator entrance
    [arnp,ainp] = desynchronize(arn,ain, time_axs,my_time,xlngth,ylngth);

    % Cavity Loss
    arL = arnp * PWR_C^0.5;
    aiL = ainp * PWR_C^0.5;

    % Calculate radiation exiting/entering undulator
    [power_z, power_s, angl] = find_pwr(ar, ai, i0,dx);
    [power_zL,power_sL,anglL] = find_pwr(arL,aiL, i0,dx);

    % Save peak power/pulse energy
    my_pls(j) = max(power_s);
    my_nrg(j) = sum(power_s) * delt;
    [my_tau(j), peakPos, ~] = find_fwhm_pwr(power_s, rad_s,rad_e,delt);

    % Save power profile/spectra
    pwr_prfl(j,:) = power_s;
    pwr_spct(j,:) = power_z;

    % disp(['desyn=',num2str(cav_d),...
    %     ...' dramp end=',num2str(nr), ...
    %     ' ver=', num2str(run_n), ...
    %     ' pulse pass=', num2str(j), ', pulse peak power=',...
    %         num2str(max(power_s)), '[W]'])


end

%disp(['All ', num2str(npass), ' passes completed!'])

%% Output struct containing simulation data
%vs pass
out.numPass  = npass;
out.peak_pwr = my_pls;
out.pls_nrg  = my_nrg;
out.pls_tau  = my_tau;

%axes
out.time_ax = time_axs;
out.frq_ax  = f_axs;

%each pass
out.power_profile = pwr_prfl;
out.power_spectra = pwr_spct;




end

