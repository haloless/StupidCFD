function [sdev] = voigt_dev(sigma)
% Deviatoric part
% Both input and output must be in Voigt form

I = [1,1,1,0,0,0]';

sm = sum(sigma(1:3)) / 3.0;
sdev = sigma - sm * I;

return
end


