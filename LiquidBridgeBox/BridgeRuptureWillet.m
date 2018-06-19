function [Hrup] = BridgeRuptureWillet(R1,R2,theta,V)

Rm = DerjaguinRadius(R1,R2);

if R2 > 0
    % size ratio must < 1
    if R2 > R1
        rr = R1 / R2;
    else
        rr = R2 / R1;
    end
else
    % sphere-flat
    rr = 0;
end

% nondimensional
Vstar = V / Rm^3;

Vhat = Vstar^(1/3);

Hstar = (1 + 0.25*theta*(rr+1)) * (Vhat + (0.5*rr-0.4)*Vhat^2);

% restore dimension
Hrup = Hstar * Rm;


return
end

