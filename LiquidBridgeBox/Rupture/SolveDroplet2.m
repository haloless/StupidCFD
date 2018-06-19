function [phi,ssph,sdrop] = SolveDroplet2(R,V, alpha)

phi0 = pi/6;

phi = fsolve(@(aa) (vol_func(R,V,alpha,aa)./V-1), phi0);

b = R * sin(alpha);
r = b / sin(phi);

hdrop = r * (1-cos(phi));
hsph = R * (1-cos(alpha));

sdrop = 2*pi*r*hdrop;
ssph = 2*pi*R*hsph;

return
end

function [vol] = vol_func(R,V, alpha, phi)
    b = R * sin(alpha);
    r = b / sin(phi);
    
    hdrop = r * (1-cos(phi));
    hsph = R * (1-cos(alpha));
    
    vdrop = SphereCapVolume(b, hdrop);
    vsph = SphereCapVolume(b, hsph);
    vol = vdrop - vsph;
return
end

