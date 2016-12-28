function [f2,f4] = InterpUL(theta)

ChannelMobGlobals;

tmp = interp1(MUL_f(:,1),MUL_f(:,2:end), theta);
f2 = tmp(1);
f4 = tmp(2);


return
end



