function [arn,ain]= inser_lens2(ar,ai,lambdaN,Z,LL,fN,dxN,desyn_s,DISP)
%% Function for Propagating the beam in a cavity
% Incorporates natural diffraction and focusing lenses

% Allows a display of the integrated transverse intensity profile, shown at
% undulator exit, before second lens, before first lens, and entering undu

%% Parameters
lambda = lambdaN;                          % radiation Wavelength  [meters]
N1 = length(ar(:,1,1));                     % Number of spatial grid points
padd=100;
N=N1+padd;                           % Number of x or y points with padding
k0 = 2 * pi / lambda;                               % Free-space wavenumber

dx = dxN;                                        % Sampling resolution in x
dy = dxN;                                        % Sampling resolution in y
Lx = dx*N;                                             % Physical size in x
Ly = dy*N;                                             % Physical size in y 

f = fN;                                        % Lens focal length (meters)

% Propagation Distance (cavity)
z1 = (Z-LL-desyn_s)/2;                      % from Undulator to Second lens
z2 = Z-desyn_s/2;                 % from second lens and back to first lens
z11=(Z-LL)/2;                                % from first lens to Undulator


%% Spatial/Frequency Grid
% Frequency coordinates in x,y
kx = (-(N-1)/2:(N-1)/2) * (2*pi / Lx);
ky = (-(N-1)/2:(N-1)/2) * (2*pi / Ly);
% Create meshgrid for frequency coordinates (kx, ky)
[KX, KY] = meshgrid(kx, ky);

% Create meshgrid for spatial coordinates 
x=(-(N-1)/2:(N-1)/2)*dx;
y=(-(N-1)/2:(N-1)/2)*dy;
[X, Y] = meshgrid(x, y);


%% Define Input Field (or TEST Gaussian Beam)
% Input Field (Field exiting the undulator)
E_in1=ar+1i*ai;
E_in=padarray(E_in1,[padd/2,padd/2],0,'both');

% Gaussian beam (test)
% w0=0.1e-3
% E_in = exp(- (X.^2 + Y.^2) / (2 * w0^2));

E_in_p=(abs(E_in)).^2;


%% Fourier Transform of Input Field
E_fft = fftshift( fftshift(  fft2(E_in),  1), 2 );


%% Step 1: Propagate the Field for Distance z1
% Free-space propagation kernel
H1 = exp(1i * z1 * sqrt(k0^2 - KX.^2 - KY.^2));

E_fft_z1 = E_fft .* H1;                                 % Apply propagation

% Transform x,y back to spatial domain
E_z1 = ifft2(  ifftshift( ifftshift(E_fft_z1,1), 2 )  ); 

% Radiation field intensity just before the second lens
E_z1_p=(abs(E_z1)).^2;


%% Step 2: Apply Lens Focusing Effect (2nd mirror)
% Thin lens transfer function
H_lens = exp(-1i * (k0 / f / 2) * (X.^2 + Y.^2));

E_lens = E_z1 .* H_lens;                     % Apply lens in spatial domain

% Transform to Fourier domain
E_fft_lens = fftshift( fftshift(  fft2(E_lens),  1), 2 );


%% Step 3: Propagate the Field for Distance z2 =Lc (backward)
% Free-space propagation kernel z2
H2 = exp(1i * z2 * sqrt(k0^2 - KX.^2 - KY.^2));

E_fft_focus = E_fft_lens .* H2;              % Apply propagation after lens

% Transform back to spatial domain
E_z2 = ifft2(  ifftshift( ifftshift(E_fft_focus,1), 2 )  );

% Radiation field intensity just before the first lens
E_z2_p=(abs(E_z2)).^2;


%% Step 2: Apply Lens Focusing Effect (1st mirror)
% Lens 2 = Lens 1
E_lens2 = E_z2 .* H_lens;                    % Apply lens in spatial domain

% Transform to Fourier domain
E_fft_lens2 = fftshift(  fftshift( fft2(E_lens2), 1 ),  2  );


%% Step 1: Propagate the Field for Distance z11
% Free-space propagation kernel z11
H11 = exp(1i * z11 * sqrt(k0^2 - KX.^2 - KY.^2));

E_fft_z1 = E_fft_lens2.* H11;          % Apply propagation in Fourier space

% Transform back to spatial domain
E_z1b = ifft2( ifftshift(  ifftshift(E_fft_z1,1),  2  ) );

% Radiation field intensity entering the undulator
E_z1b_p=(abs(E_z1b)).^2;


% Function Returns the fields
arn=real(E_z1b);                                                           
ain=imag(E_z1b);

% De-pad field arrays
arn=arn((padd/2+1):N-padd/2,(padd/2+1):N-padd/2,:);
ain=ain((padd/2+1):N-padd/2,(padd/2+1):N-padd/2,:);


%% Visualization
if (DISP)
    ff= figure(81);

subplot(1,4,1, 'Parent', ff); 
imagesc(x,y,sum(E_in_p,3)); title('After Undulator'); axis equal tight;

title2 = "After z1 "+num2str(z1)+" Propagation"; titl2 = char(title2);
subplot(1,4,2, 'Parent', ff); imagesc(x,y,sum(E_z1_p,3)); title(titl2); 

title3="After mirror & Lc"+num2str(z2)+" Propagation"; titl3= char(title3);
subplot(1,4,3, 'Parent', ff); imagesc(x,y,sum(E_z2_p,3)); title(titl3); 

title4="Front mirror Propogation z11 "+num2str(z11)+"(Entering Undulator)"; 
titl4 = char(title4);
subplot(1,4,4, 'Parent', ff); imagesc(x,y,sum(E_z1b_p,3)); title(titl4);

colormap jet; colorbar;

end

end
