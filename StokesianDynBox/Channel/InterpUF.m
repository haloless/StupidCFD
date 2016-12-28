function [f1,f3,f5,g1,g3,g5] = InterpUF(theta)

ChannelMobGlobals;

tmp = interp1(MUF_f(:,1),MUF_f(:,2:end), theta);
f1 = tmp(1);
f3 = tmp(2);
f5 = tmp(3);

tmp = interp1(MUF_g(:,1),MUF_g(:,2:end), theta);
g1 = tmp(1);
g3 = tmp(2);
g5 = tmp(3);


return
end



