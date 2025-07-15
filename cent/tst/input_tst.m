function [E1,E2,delt,i0,dx,rad_s,rad_e,ebm_s,ebm_e]= input_tst(reflect,ar,ai,inp_struc)
%% function for Generating GINGER FIELD files
my_dir='/Users/amirweinberg/ginger_int/tst/';                           %Working Directory 


if reflect==1                                                              % First pass through undulator
    %myinfile = strcat(my_dir, 'ginger/G3din',num2str(reflect));
    %MK3_make_G3d_in_file(myinfile,inp_struc);                              % Generate in file for GINGER3D

    %myoutfile = 'UHM_G3D1';                                                % Name outfile
    %myCommand = "COLUMNS=1000 /Users/amirweinberg/Documents/ginger/3d/xg3D-static-1.1.0-a5  -i "+myinfile+" -r "+myoutfile;              %Command for running GINGER3D
    
    %myCommand = char(myCommand);
    %[status,cmdout] = system(['osascript -e ''tell application "Terminal" to do script "' myCommand '; exec bash"''']);    %Open Terminal debug purposes
    %[status,cmdout] = system([myCommand]);                                 % Run command
    %pause(.1)
    out_file=my_dir+"fldUHM_G3D"+num2str(reflect)+".hdf5";
    hdf_file=my_dir+"UHM_G3D"+num2str(reflect)+".hdf5";
    %h5disp(out_file);
    %h5disp(hdf_file);
    real_params= h5read(hdf_file, '/base_param/real_params');
    %h5disp(out_file);
    
    real_param_buf = h5read(out_file, '/header/real_param_buf');           % Get paramters 
    fund_field = h5read(out_file, '/radiation/2D_xy_field');               % Get Field
    
    E=fund_field;
    E1=E(1:end/2,:,:);                                                     % get  E field to real data 
    E2=E(end/2+1:end,:,:);                                                 % get  E field imaginary data
    % delr=real_param_buf(9);
    i0=real_param_buf(2);                                                  % get power coefficient 
    delt=real_param_buf(5);                                                % get time delta
    dx=real_param_buf(9);                                                  % get delta in x grid
    %y=real_param_buf(8);
    rad_s=real_param_buf(16)
    rad_e=real_param_buf(17)
    ebm_s=real_params.value1(15);
    ebm_e=real_params.value1(15);


else
     %if (mod(reflect,20))
     %    reflect=2;
     %end
    fld_file=my_dir+"fldUHM_G3D11.hdf5";                                    % Input field file
    E3 = cat(1, ar, ai);                                                   % Concatenate complex field
    h5write(fld_file,'/radiation/2D_xy_field',E3);                         % Write back to field file

    %myinfile = strcat(my_dir,'ginger/G3d2in',num2str(reflect));                   % Name the in file 
    %MK3_make_G3d_2in_file(myinfile,inp_struc);                       % Generate in file for 2 pass in oscillator and on
    %myoutfile = strcat('UHM_G3D2',num2str(reflect));                       % Name out file
    
    %myCommand = "COLUMNS=400 /Users/amirweinberg/Documents/ginger/3d/xg3D-static-1.1.0-a5 i="+myinfile+" f="+fld_file+" r="+myoutfile;  %Command for running GINGER3D
    %myCommand = char(myCommand);
       
  % [status,cmdout] = system(['osascript -e ''tell application "Terminal" to do script "' myCommand '; exec bash"''']);    %Open Terminal Debug purposes
    %[status,cmdout] = system([myCommand]);                                 % Run command

    out_file=my_dir+"fldUHM_G3D2"+num2str(reflect)+".hdf5";
    hdf_file=my_dir+"UHM_G3D2"+num2str(reflect)+".hdf5";
    %h5disp(out_file);
    %h5disp(hdf_file);
    real_params= h5read(hdf_file, '/base_param/real_params');
    %onaxis_field = h5read(hdf_file, '/radiation/scalar_data');
   % input_env_vars = h5read(hdf_file, '/input/input_env_vars'); 
    real_param_buf = h5read(out_file, '/header/real_param_buf');
    fund_field = h5read(out_file, '/radiation/2D_xy_field');

    E=fund_field;
    E1=E(1:end/2,:,:);                                                     % get  E field to real  data 
    E2=E(end/2+1:end,:,:);                                                 % get  E field imaginary data
    % delr=real_param_buf(9);
    i0=real_param_buf(2);                                                  % get power coefficient 
    dx=real_param_buf(10);                                                 % get delta in x grid
    delt=real_param_buf(5);                                                % get time delta
    rad_s=real_param_buf(16)
    rad_e=real_param_buf(17)
    ebm_s=real_params.value1(15);
    ebm_e=real_params.value1(14);

    %rads_s=real_param_buf(18);
    %rads_e=real_param_buf(19);
    % figure(1)
    % hold off;
    % plot([E1' E2']')
    % title('input')

end
end
