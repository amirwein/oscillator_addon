function [w0SOL, ZrSOL] = Zray(w11, w22, LL, unduL, desyn_s)
%% This function calculates and returns the Zrayleigh and waist size 

z1=unduL/2;                                            
z2= z1+(LL-unduL-desyn_s)/2;

syms w0 Zr
eq1 = w0*sqrt(1+(z1/Zr)^2)==abs(w11);
eq2 = w0*sqrt(1+(z2/Zr)^2)==abs(w22);

solution = solve([eq1, eq2], [w0, Zr]);

w0_sol=solution.w0;
Zr_sol=solution.Zr;
w0SOL=double(w0_sol)
ZrSOL=double(Zr_sol)
end
