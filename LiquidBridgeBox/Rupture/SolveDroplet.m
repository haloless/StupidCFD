function [alpha,ssph,sdrop] = SolveDroplet(R,theta,V)

alpha0 = pi/6;

alpha = fsolve(@(aa) (vol_func(R,theta,V,aa)./V-1), alpha0);

b = R * sin(alpha);
r = b / sin(alpha+theta);

hdrop = r * (1-cos(alpha+theta));
hsph = R * (1-cos(alpha));

sdrop = 2*pi*r*hdrop;
ssph = 2*pi*R*hsph;

return
end

function [vol] = vol_func(R,theta,V, alpha)
    b = R * sin(alpha);
    r = b / sin(alpha+theta);
    
    hdrop = r * (1-cos(alpha+theta));
    hsph = R * (1-cos(alpha));
    
    vdrop = SphereCapVolume(b, hdrop);
    vsph = SphereCapVolume(b, hsph);
    vol = vdrop - vsph;
return
end

