function [power_z,power_s,angl]= find_pwr(ar,ai,i0,dx)
%% Function to calculate the Radiation pulse in time and frequency domain

 % For Power-profile of Radiation pulse
 power_s = abs(ar+1i*ai).^2;
 power_s = sum(sum(power_s,1),2);                          % Sum over dx dy 
 power_s = squeeze(power_s(1,1,:))*i0*dx^2;
 
 % For Phase of Radiation pulse 
 Ea=ar+1i*ai;
 Ea= sum(sum(Ea,1),2);
 angl=squeeze(Ea(1,1,:));
 angl=unwrap(angle(angl));

 % For Power-spectrum of Radiation pulse
 E=ar+1i*ai;    
 Ef=abs(fftshift(fft(E, [], 3))).^2;                                         
 power_z=sum(sum(Ef,1),2);
 power_z=squeeze(power_z(1,1,:));
end