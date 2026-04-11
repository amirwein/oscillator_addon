function MK3_make_G3d_2in_file(filename,inp_struc)
%% Function file for making a GINGER3D "in" file for the UHM FEL Oscillator 2nd pass and on
%% Refer to GINGER3D manual for variable meaning

fileID = fopen(filename, 'w');
fprintf(fileID, 'input file for parameters corresponding to UHM IR FEL benchmark$$$\n');
fprintf(fileID, '&in_beams\n');
fprintf(fileID, '!>>> e-beam parameters:\n');
fprintf(fileID, 'emitx_mks = %4.1e\n',inp_struc.emitx);
fprintf(fileID, 'emity_mks = %4.1e\n',inp_struc.emity);
% 4-D transverse phase space Gaussian distribution
fprintf(fileID, 'jmg = %2d\n', 2);
fprintf(fileID, 'betaX_twiss = %4.4f\n',inp_struc.betaX);
fprintf(fileID, 'betaY_twiss = %4.4f\n', inp_struc.betaY);
fprintf(fileID, 'alphaX_twiss = %4.4f\n', inp_struc.alphaX);
fprintf(fileID, 'alphaY_twiss = %4.4f\n', inp_struc.alphaY);
fprintf(fileID, 'energy_MeV = %4.4f\n',inp_struc.energy);       % unit: MeV
fprintf(fileID, 'sigmaE_norm = %4.4f\n',inp_struc.eSpread);
fprintf(fileID, 'gam_load = ''gaussian''\n');
% ^Gaussian energy spread distribution

fprintf(fileID, 'current = %4.4f\n',inp_struc.currentMax);
% gaussian current shape 
fprintf(fileID, 'ebeam_pulse_shape=''gaussian''\n');
% RMS electron beam temporal pulse width
fprintf(fileID, 'ebeam_tbody=%1.4e\n',inp_struc.RMSwidth);
fprintf(fileID, 'ntestp = %4d\n',8192);               % # of macroparticles 
fprintf(fileID, 'nfold_sym = %4d\n',8);               	    % default value
fprintf(fileID, 'lshot = .t\n');   

% 
%!> shot noise *on*
fprintf(fileID, '\n!>>> undulator parameters:\n');
fprintf(fileID, 'wavelw = %4.4f\n',inp_struc.unduPeriod);    % unit: meters
fprintf(fileID, 'llinear = .t\n'); 	 	          % linearly-polarized undu
fprintf(fileID, 'idesign = %2d\n',1); 	           % constant strength undu
fprintf(fileID, 'aw0 = %4.4f\n',inp_struc.RMSunduK); 	   % rms value of K 
fprintf(fileID, 'rkxkw = %4.4f\n', 0.00);
fprintf(fileID, '\n!>>> input laser characteristics:\n');
fprintf(fileID, 'wavels = %1.4e\n',inp_struc.radWvlength);     	  % unit: m
fprintf(fileID, 'plaser = %1.1e\n',0);              	      % unit: watts
fprintf(fileID, '/END\n');

fprintf(fileID, '&in_sim_diag\n');
fprintf(fileID, 'l_write_fld=.t\n');
fprintf(fileID, 'l_write_rst=.t\n'); % write a restart file for phase space
%!window ≡ nphoton × dt slice = (nside / nsidep) × (Lw × λs /λw )

% how many e-beam slices you need
fprintf(fileID, 'nside = %4d\n',inp_struc.eslices);
%192 how many radiation slices, should be a power of two
fprintf(fileID, 'nslice_sim = %4d\n',inp_struc.rslices);
fprintf(fileID, 'lshort=.t\n');%
% 15 totall # of Photon slices each elec beam slice will interact 
fprintf(fileID, 'nsidep = %4d\n',inp_struc.intslices); 		                               
fprintf(fileID, ' iseeds(1:2)=20241810	12418367\n');



fprintf(fileID, '\n!>> simulation and diagnostic parameters:\n');
fprintf(fileID, 'lfred = .f\n');
fprintf(fileID, 'l_2Dxy=.t\n');
fprintf(fileID, 'zmxmeter = %4.4f\n',inp_struc.unduL);   	% unit = meters
fprintf(fileID, 'dzstep = %4.4f\n',inp_struc.unduPeriod);
fprintf(fileID, '!rmaxsim = %4.4f\n',0.025);    % grid boundary unit = cm!!
fprintf(fileID, 'dxg = %1.4e\n',inp_struc.rgrid);    	         % unit = m
fprintf(fileID, 'nxg = %4d\n',inp_struc.ngrid);
fprintf(fileID, '/END\n');
fclose(fileID);
end