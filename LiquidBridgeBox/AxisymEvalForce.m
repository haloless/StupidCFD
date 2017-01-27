%% Eval liquid bridge action by pressure and tension.
%% Attractive +, Repulsive -
%%  alpha: embracing angle, wall=0
%%  theta: contact angle
%%  r: contact ring, not sphere!!
%%  pres:
%%  sigma: 
function [f,fp,ft] = AxisymEvalForce(alpha,theta,r,pres,sigma)

% laplace pressure
fp = -pi*r^2*pres;

% tension on contact ring
ft = 2*pi*r*sin(alpha+theta) * sigma;

% total force
f = fp + ft;


return
end
