%% 
%% Assume the bridge is a straight cylinder connecting two spheres.
%% Calculate the cylinder radius for given volume.
%%
function [alpha1,alpha2,r1,r2] = SolveFlatBridge(R1,R2,H,theta1,theta2,V)

assert(theta1+theta2<pi, 'Invalid contact angles');

disp(mfilename);

tol = 1.0e-3;
% normalize by input volume for numerical reason
vfunc = @(a1) (1/V).*FlatVolume(R1,R2,H,theta1,theta2,a1) - 1;

if R2 > 0
    % particle-particle, embracing angle need solve 
    
    % initial guess
    a1guess = (pi-theta1-theta2)/2;
    % Newton's method
    alpha1 = SolveFunc(vfrunc, a1guess, tol);
    
    alpha2 = pi - theta1 - theta2 - alpha1;
    
    r1 = R1 * sin(alpha1);
    r2 = R2 * sin(alpha2);
else
    % particle-wall, embracing angle known a priori
    alpha1 = pi - theta1 - theta2;
    alpha2 = 0;
    
    r1 = R1 * sin(alpha1);
    r2 = (R1*(1-cos(alpha1)) + H) / tan(theta2) + r1;
    
    if 1
        % check volume ok
        volerr = vfunc(alpha1);
        if abs(volerr) > 1.0e-2
            warning(['Particle-Wall volume error: ', num2str(volerr)]);
        end
    end
end

return
end

function [vol] = FlatVolume(R1,R2,H,theta1,theta2, alpha1)
    
    r1 = R1 * sin(alpha1);
    h1 = R1 * (1-cos(alpha1));
    vol1 = SphereCapVolume(r1,h1);
    
    if R2 > 0
        alpha2 = pi-theta1-theta2-alpha1;
        r2 = R2 * sin(alhpa2);
        h2 = R2 * (1-cos(alpha2));
        vol2 = SphereCapVolume(r1,h2);
    else
        alpha2 = 0;
        r2 = (h1+H)/tan(theta2) + r1;
        h2 = 0;
        vol2 = 0;
    end
    
    
    
    % frustum cone volume
    vol = ConeVolume(r1,r2,H+h1+h2);
    % minus two sphere cap
    vol = vol - vol1 - vol2;
    
    return
end


