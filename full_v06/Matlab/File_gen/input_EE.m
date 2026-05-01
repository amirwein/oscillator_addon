function [E1,E2,delt,i0,dx,rad_s,rad_e,ebm_s,ebm_e]= input_EE(reflect,ar,ai,inp_struc,EFF)
%% function for Generating GINGER-3D files


% Edit according to user
home_dir='/Users/levikf/';                                  %User Directory
my_dir=home_dir+"full_v06/Matlab/";                      %Working Directory
file_indir='File_gen/in_dir/';                        %Input file directory
file_outdir='out_dir/';                 %New directory to save ginger files

% Absolute path for output-file directory
out_path= string( strcat(my_dir, file_outdir) );

%Absolute path to Ginger program executable
exe_path= '/Users/amirweinberg/Documents/ginger/3d/xg3D-static-1.1.0-a5';  

%EFF = 0;                                        % Toggle efficient storage

% Fresh vs. modulated beam
if reflect==1                                % First pass through undulator
    % Clear old input/output files
    delete(strcat(file_indir,'G3d*'));
    delete(strcat(file_outdir,'*.hdf5'));

    % Create input file
    myinfile = strcat(file_indir,'G3din',num2str(reflect));
    in_path = strcat(my_dir, myinfile);
    MK3_make_G3d_in_file(in_path,inp_struc);

    % Name to run Ginger3D
    myoutfile = 'UHM_G3D1';
    myCommand = "COLUMNS=1000 "+exe_path+" -i "+myinfile+" -r "+myoutfile;


else                                                         % Later passes
    passnum = reflect;
    % Save fewer files, once every ## passes
     if EFF && (mod(reflect,20))
         reflect=2; 
     end
    % File used for input radiation field
    fld_file = file_outdir+"tmp_fld.hdf5";
    
    % Concatenate complex field in Ginger convention
    E3 = cat(1, ar, ai);
    
    % Write back to field input file
    h5write(fld_file,'/radiation/2D_xy_field',E3);

    % Name the in-file for 2+ passes in oscillator
    myinfile = strcat(file_indir,'G3d2in',num2str(reflect));
    in_path = strcat(my_dir, myinfile);
    MK3_make_G3d_2in_file(in_path,inp_struc,passnum);
    
    % Name of run
    myoutfile = strcat('UHM_G3D2',num2str(reflect));
    
    %Command for running GINGER3D
    myCommand = "COLUMNS=400 "+exe_path+" -i "+myinfile+...
        " -f "+fld_file+" -r "+myoutfile;

end


%% Run Ginger3D Simulation
myCommand = char(myCommand);
[status,cmdout] = system(myCommand);                          % Run command

% All files for a run are generated in the working directory. Here we move
% files to tmp (output) storage
file_suffix = myoutfile+".hdf5";
movefile("*"+file_suffix, out_path);

% Point to newly-generated Ginger field/out file
field_output_file=out_path+"fld"+myoutfile+".hdf5";
ginger3d_out_file=out_path+myoutfile+".hdf5";

% Create temp_file for Ginger field-file structure (for future passes)
copyfile(field_output_file, out_path+"tmp_fld.hdf5");


%% Func output variables
% Get parameters for beam from simulation
real_params= h5read(ginger3d_out_file, '/base_param/real_params');
ebm_s=real_params.value1(15);          % start time of electron-beam  [sec]
ebm_e=real_params.value2(15);            % end time of electron-beam  [sec]

% Get parameters for radiation field
real_param_buf = h5read(field_output_file, '/header/real_param_buf');
i0=real_param_buf(2);             % intensity->power coefficient % [W/m**2]
delt=real_param_buf(5);                               % slice spacing [sec]
dx=real_param_buf(9);         % transverse grid size, assume same y [meter]
rad_s=real_param_buf(16);             %start time of radiation % head [sec]
rad_e=real_param_buf(17);             % end time of radiation % tail  [sec]

% Extract output field amplitudes
fund_field = h5read(field_output_file, '/radiation/2D_xy_field');
E=fund_field;
E1=E(1:end/2,:,:);                         % get real components of E field
E2=E(end/2+1:end,:,:);                         % get E field imaginary data


end
