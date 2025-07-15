function [arn,ain,wa1,wa2]= inser_lens2(ar,ai,lambdaN,Z,LL,fN,dxN,desyn_s)
%% Function for Propgating the beam in a cavity, Natural diffraction and focusing lenses 

%% Parameters
lambda = lambdaN;                                                          % Wavelength (meters)
N = length(ar(:,1,1));                                                     % Grid size
mid=floor(length(ar(1,1,:))/2);
k0 = 2 * pi / lambda;                                                      % Free-space wavenumber

dx = dxN;                                                                  % Sampling resolution in x
dy = dxN;                                                                  % Sampling resolution in y
Lx = dx*N;                                                                 % Physical size in x (1 mm)
Ly = dy*N;                                                                 % Physical size in y (1 mm)

f = fN;                                                                    % Lens focal length (meters)
z1 = (Z-LL-desyn_s)/2;                                                     % Propagation distance from Undulator to Second lens (at Z) (LL-Undul length, desyn_s-desynchornization shift)
z2 = Z-desyn_s/2;                                                          % Propagation distance from second lens and back to first lens
z11=(Z-LL)/2;                                                              % Propagation from first lens to Undulator


%% Spatial Frequency Grid (kx, ky)
kx = (-N/2:N/2-1) * (2*pi / Lx);                                           % Frequency coordinates in x
ky = (-N/2:N/2-1) * (2*pi / Ly);                                           % Frequency coordinates in y
[KX, KY] = meshgrid(kx, ky);                                               % Create meshgrid for (kx, ky)

x = linspace(-Lx/2, Lx/2, N);
xtk= linspace(-Lx/2, Lx/2, N/10);
y = linspace(-Ly/2, Ly/2, N);
ytk= linspace(-Ly/2, Ly/2, N/10);
[X, Y] = meshgrid(x, y);                                                   % Create meshgrid for spatial coordinates  



%% Define Input Field (can TEST Gaussian Beam)
E_in=ar+1i*ai;                                                             % Input Field
% w0=0.1e-3                     
% E_in = exp(- (X.^2 + Y.^2) / (2 * w0^2));                                % Gaussian beam

myE=squeeze(E_in(:,:,mid));
[waist_x_in, waist_y_in]=find_waistxy(myE, N,x,y);


%% Fourier Transform of Input Field
E_fft = fftshift(fft2(E_in));

%% Step 1: Propagate the Field for Distance z1 (w desynchronization) 
H1 = exp(1i * z1 * sqrt(k0^2 - KX.^2 - KY.^2));                            % Free-space propagation kernel z1
E_fft_z1 = E_fft .* H1;                                                    % Apply propagation
E_z1 = ifft2(ifftshift(E_fft_z1));                                         % Transform back to spatial domain

myE=squeeze(E_z1(:,:,mid));
[waist_x_z1, waist_y_z1]=find_waistxy(myE, N,x,y)


%% Step 2: Apply Lens Focusing Effect (2nd mirror)
H_lens = exp(-1i * (k0 / f / 2) * (X.^2 + Y.^2));                          % Thin lens transfer function
E_fft_lens = E_z1 .* H_lens;                                               % Apply lens in spatial domain
E_lens = fftshift(fft2(E_fft_lens));                                       % Transform to Fourier domain

%% Step 3: Propagate the Field for Distance z2 =Lc (backward)  (w desynchronization)  
H2 = exp(1i * z2 * sqrt(k0^2 - KX.^2 - KY.^2));                            % Free-space propagation kernel z2
E_fft_focus = E_lens .* H2;                                                % Apply propagation after lens
E_z2 = ifft2(ifftshift(E_fft_focus));                                      % Transform back to spatial domain

myE=squeeze(E_z2(:,:,mid));
[waist_x_z2 waist_y_z2]=find_waistxy(myE, N,x,y);

%% Step 2: Apply Lens Focusing Effect (1st mirror)
E_fft_lens2 = E_z2 .* H_lens;                                              % Apply lens in spatial domain
E_lens2 = fftshift(fft2(E_fft_lens2));                                     % Transform to Fourier domain

%% Step 1: Propagate the Field for Distance z1 again (w/o desynchronization) 
H11 = exp(1i * z11 * sqrt(k0^2 - KX.^2 - KY.^2));                          % Free-space propagation kernel z11
E_fft_z1 = E_lens2.* H11;                                                  % Apply propagation in Fourier space
E_z1b = ifft2(ifftshift(E_fft_z1));                                        % Transform back to spatial domain

myE=squeeze(E_z1b(:,:,mid));
[waist_x_z1b waist_y_z1b]=find_waistxy(myE, N,x,y);

arn=real(E_z1b);                                                           %Fucntion Return back the fields
ain=imag(E_z1b);

wa1=waist_x_in
wa2=waist_x_z1

%% Visualization
figure(2);
subplot(1,4,1); imagesc(x,y,squeeze(abs(E_in(:,:,mid)))); title('After Undulator'); axis equal tight; grid on; set(gca, 'XTick',xtk); set(gca, 'YTick',ytk); xline(waist_x_in,'k',num2str(waist_x_in),'LineWidth', 1.5); yline(waist_y_in,'k',num2str(waist_y_in),'LineWidth', 1.5);
title2 = "After z1 "+num2str(z1)+" Propagation"; 
titl2 = char(title2);
subplot(1,4,2); imagesc(x,y,squeeze(abs(E_z1(:,:,mid)))); title(titl2); axis equal tight; grid on; set(gca, 'XTick',xtk); set(gca, 'YTick',ytk); xline(waist_x_z1,'k',num2str(waist_x_z1),'LineWidth', 1.5); yline(waist_y_z1,'k',num2str(waist_y_z1),'LineWidth', 1.5);
title3 = "After second mirror & Lc"+num2str(z2)+" Propagation"; 
titl3 = char(title3);
subplot(1,4,3); imagesc(x,y,squeeze(abs(E_z2(:,:,mid)))); title(titl3); axis equal tight; grid on; set(gca, 'XTick',xtk); set(gca, 'YTick',ytk); xline(waist_x_z2,'k',num2str(waist_x_z2),'LineWidth', 1.5); yline(waist_y_z2,'k',num2str(waist_y_z2),'LineWidth', 1.5);
title4 = "After first mirror & z11 "+num2str(z11)+" Propagation before undulator"; 
titl4 = char(title4);
subplot(1,4,4); imagesc(x,y,squeeze(abs(E_z1b(:,:,mid)))); title(titl4); axis equal tight; grid on; set(gca, 'XTick',xtk); set(gca, 'YTick',ytk); xline(waist_x_z1b,'k',num2str(waist_x_z1b),'LineWidth', 1.5); yline(waist_y_z1b,'k',num2str(waist_y_z1b),'LineWidth', 1.5);

colormap jet; colorbar;
% %% Visualization E gaussain
% figure(2);
% subplot(1,4,1); imagesc(squeeze(abs(E_in))); title('Input Gaussian Beam'); axis equal tight;
% subplot(1,4,2); imagesc(squeeze(abs(E_z1))); title('After z1 Propagation'); axis equal tight;
% subplot(1,4,3); imagesc(squeeze(abs(E_z2))); title('After Lens & Lc Propagation'); axis equal tight;
% subplot(1,4,4); imagesc(squeeze(abs(E_z1b))); title('After 2Lens & z1 Propagation'); axis equal tight;
% colormap jet; colorbar;