%% Directly calculate the force from the Laplace part and the tension part.
%%
function [F1,Fp1,Ft1,F2,Fp2,Ft2] = AxisymEvolverEvalForce(bridge, np,rp,xp,pres)

R1 = bridge.R1;
R2 = bridge.R2;
% H = bridge.H;
theta1 = bridge.theta1;
theta2 = bridge.theta2;
sigma = bridge.sigma;

r1 = rp(1);
alpha1 = AxisymEmbraceAngle(R1,r1);
[F1,Fp1,Ft1] = AxisymEvalForce(alpha1,theta1,r1,pres,sigma);

r2 = rp(np);
alpha2 = AxisymEmbraceAngle(R2,r2);
[F2,Fp2,Ft2] = AxisymEvalForce(alpha2,theta2,r2,pres,sigma);

return
end


