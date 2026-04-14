%% Main file for FEL Oscillator add on for GINGER3d
% Allows for electron time-of-flight jitter, cavity desynchronization,
% and electron-bunch desynchronization

% Reset vars,figures before run
clear;
close all;

% Include directories containing functions
addpath("Util_funs","File_gen");

% Before continuing, go to the input_EE function in File_gen to point your
% user/input/output file directories, ginger executable, run name


%% Config toggles

% Live figs showing the transverse intensity profile
DISP=1;

% Efficient storage
EFF=0;


%% General params
c = 299792458;
m_e = 0.511;                                     % electron rest mass [MeV]

%% UHM Beam parameters
energy  = 40;        				       % electron kinestic energy [MeV]
eSpread = 0.5e-2;       					 % fractional rms energy spread 
emitN   = 8e-6;                 % normalized transverse emittance [mm-mrad]
eGamma  = 1 + energy/m_e;                % total relativistic e-beam energy
currentMax = 37.5;					                % peak current [Ampere]
betaX_twiss = 1.4;                                       % twiss parameters
betaY_twiss = 0.24;
alphaX_twiss = 0.4714;
alphaY_twiss = 0.0;
% e-beam bunch rms Gaussian pulse width
bunch_len = 2e-12;                            % fwhm - pulse duration [sec]
bunch_sigma = bunch_len/sqrt( log(256) );         % rms - pulse width [sec]
%bunch_sigma = 0.8493e-12;                      % 2ps pulse uses 0.8493e-12
%bunch_sigma = 0.215e-12;                      % 0.5ps pulse uses 0.215e-12


%% UHM FEL parameters
unduPeriod = 0.023;     						% undulator period  [meter]
unduK = 1.2/sqrt(2);    				     % RMS undulator parameter, (K) 
NumPeriod = 47;                               % Number of undulator periods
unduL = unduPeriod*NumPeriod;  		   % total length of undulator  [meter]
%radWavelength = 3.2281e-06;                 % radiation wavelength  [meter]
radWavelength = unduPeriod/2/eGamma^2 * (1 + unduK^2); 
fcar =c/radWavelength;                            % carrier frequency  [hZ]
slipp = radWavelength*NumPeriod;                 % Slippage length  [meter]


%% UHM Oscillator parameters
CavityLength= 2.0469;                                       % Cavity length
MirrorRadiusCurv=1.3;                        % Mirrors Radius of Curvatures
PWR_C = 0.93;                                         % CAVITY Reflectivity    
npass=400;                                 % Number of passes in the cavity

% Electron Jitter
jitterRMS = 0;%0.4e-12;

% Cavity Desync     % desync units = 2 (change in cavity length) / slippage
cav_d=0;

% Electron Desync (ramping)                       % Given in units [desync]
%dmax = 0;
%d_ramp = linspace(dmax, 1,npass);       % Linear Ramping desynchronization

%pass_i = 1:npass; pass_i = (pass_i-1)./399;
%d_ramp = 399*0.005/199 * (1 - pass_i);         % Dynamic desynchronization

d_ramp = zeros(1,npass);                    % No electron desynchronization


%% GINGER3D simulation paramteres
ebunch_slices = 192;                     % Number of electron bunch slices
%ebunch_slices = 96;                          % 0.5ps uses 96, 2ps uses 192
rad_slices = 256;                              % Number of Radiaion slices
%rad_slices = 192;                           % 0.5ps uses 192, 2ps uses 256

% Number of Photon slices each elec beam slice will interact with   %15
intrct_slices = 15;

r_grid = 0.080e-3;                    % Spatial (x,y) grid spacing  [meter]

% Number of grid points along x or y
n_grid = 129;                       % odd value ensures even # of divisions


% Create structural object containing input parameters
inp_struc.energy      = energy;
inp_struc.eSpread     = eSpread;
inp_struc.emitx       = emitN;
inp_struc.emity       = emitN;
inp_struc.currentMax  = currentMax;
inp_struc.betaX       = betaX_twiss;
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


%% Setup Axes/First Run

pass_axs = 1:npass;

my_pls=zeros(1,npass);

% Run first pass GINGER3d for params
[ar1, ai1,delt,i0,dx1,rad_s,rad_e,ebm_s,ebm_e] = input_EE(1,1,1,...
                                                            inp_struc,EFF);   

% Number of gridpoints along each dimension
xlngth=length(ar1(:,1,1));
ylngth=length(ar1(1,:,1));

% Cavity desync in meters, seconds
desyn_shft = cav_d * slipp;
dshift_s = desyn_shft/c;

% Time,Frequency axes
time_axs = (rad_e-delt/2:-delt:rad_s+delt/2);
f_axs = (-0:(length(time_axs)-1))/(time_axs(1)-time_axs(end));

% Shift frequency axis to carrier frequency
f_axs = f_axs-median(f_axs)+fcar;

% Time frame of field entering undulator
my_time = time_axs+dshift_s;


%% Finish 1st pass, Start later passes
% Implement Cavity
[arn1,ain1]=inser_lens2(ar1,ai1,radWavelength,CavityLength,unduL,...
                        MirrorRadiusCurv/2,dx1,desyn_shft,DISP);

% Shift radiation field time-frame (due to cavity desync)
[ars1, ais1] = desynchronize(arn1,ain1,time_axs,my_time,xlngth,ylngth);

% Cavity loss
arL = ars1 * PWR_C^0.5;
aiL = ais1 * PWR_C^0.5;

% Calculate radiation pulse entering undulator
[power_z,power_s,angl] = find_pwr(arL,aiL,i0,dx1);
my_pls(1) = max(power_s);                       % Save peak radiation power

% Display input power
disp(['pulse pass=', num2str(1), ' pulse peak power=',...
        num2str(my_pls(1)), '[W]'])

for j=2:npass
    % Get shift for current pass
    my_shft = jitterRMS*randn() + d_ramp(j)*slipp/c;
    new_time = time_axs + my_shft;

    % Shift radiation to current frame
    [ars,ais] = desynchronize(arL,aiL,time_axs,new_time,xlngth,ylngth);

    % Call Ginger-3D (centered on e-beam) for pass j
    [ar, ai,delt,i0,dx,rad_s,rad_e,ebm_s,ebm_e]=input_EE(j, ars, ais, ...
                                                            inp_struc,EFF);

    % Reverse timing-offset shift
    [arns, ains] = desynchronize(ar,ai, new_time, time_axs,xlngth,ylngth);

    % Implement Cavity
    [arn,ain] = inser_lens2(arns,ains,radWavelength,CavityLength,unduL,...
                            MirrorRadiusCurv/2,dx,desyn_shft,DISP);

    % Shift field to undulator entrance
    [arnp,ainp] = desynchronize(arn,ain, time_axs,my_time,xlngth,ylngth);

    % Cavity Loss
    arL = arnp * PWR_C^0.5;
    aiL = ainp * PWR_C^0.5;

    % Calculate radiation entering undulator
    [power_z,power_s,angl] = find_pwr(arL,aiL, i0,dx);
    my_pls(j) = max(power_s);

    disp(['pulse pass=', num2str(j), ', pulse peak power=',...
            num2str(max(power_s)), '[W]'])


end


disp(['All ', num2str(npass), ' passes completed!'])



