%
function [F1,F1p,F1t,F2,F2p,F2t] = AxisymEvolverDirectForce(bridge, np,rp,xp)

% exec evolver
[rp,xp,pres] = AxisymEvolver(bridge, np,rp,xp);

[F1,F1p,F1t,F2,F2p,F2t] = AxisymEvolverEvalForce(bridge, np,rp,xp,pres);
[~,~,~,F2,F2p,F2t] = AxisymEvolverEvalForce(bridge, np,rp,xp,-pres);

return
end

