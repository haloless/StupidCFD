function [F1,F1p,F1t,F2,F2p,F2t] = AxisymEvolverDirectForce(bridge, np,rp,xp)
% Obtain force by direct integration

% exec evolver
[rp,xp,pres] = AxisymEvolver(bridge, np,rp,xp);

% eval force
[F1,F1p,F1t,F2,F2p,F2t] = AxisymEvolverEvalForce(bridge, np,rp,xp,pres);

return
end

