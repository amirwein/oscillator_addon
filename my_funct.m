function [power_z,power_s,angl]= my_funct(ar,ai,i0,dx)
%% Function to calculate the Radiation pulse in time and frequency domain
%ar real part electric field
%ai imaginary part electric field
%i0 normalization constant from ginger
%dx delt transverse grid

 power_s = abs(ar+1i*ai).^2;
 power_s = sum(sum(power_s,1),2);                                          % Sum over dx dy 
 power_s = squeeze(power_s(1,1,:))*i0*dx^2;
 
 Ea=ar+1i*ai;
 Ea= sum(sum(Ea,1),2);
 angl=squeeze(Ea(1,1,:));
 angl=unwrap(angle(angl));
 
 E=ar+1i*ai;  
 Ef=abs(fftshift(fft(E, [], 3))).^2;                                       % For Spectrum of Radiation pulse    
 power_z=sum(sum(Ef,1),2);
 power_z=squeeze(power_z(1,1,:));
end
