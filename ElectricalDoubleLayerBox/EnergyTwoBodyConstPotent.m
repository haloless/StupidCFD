function [u] = EnergyTwoBodyConstPotent(a1,a2,H,kappa,cutoff, acoef,bcoef, phi1,phi2)

% currently only equal sphere
if a1~=a2 || phi1~=phi2
	error('Equal sphere only!');
end

if 1
    ka = kappa * a1;
    kh = kappa * H;
    piso = phi1;
    a0 = acoef(1);
    aend = acoef(end);
    
    u = 4*pi * (-(pi/2)*(a0/piso)*csch(ka) + ka + ka*coth(ka));
end


return
end


