function [f2,f4,f6,g2,g4,g6] = InterpUS(theta)

ChannelMobGlobals;

tmp = interp1(MUS_f(:,1), MUS_f(:,2:end), theta);
f2 = tmp(1);
f4 = tmp(2);
f6 = tmp(3);

tmp = interp1(MUS_g(:,1), MUS_g(:,2:end), theta);
g2 = tmp(1);
g4 = tmp(2);
g6 = tmp(3);

return
end

