function [hrup] = BridgeRuptureLian(theta,V)
hrup = (1+0.5*theta) * V^(1.0/3.0);
return
end

