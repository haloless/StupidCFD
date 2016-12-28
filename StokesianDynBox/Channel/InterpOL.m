function [f3,g3] = InterpOL(theta)

ChannelMobGlobals;

tmp = interp1(MOL_fg(:,1), MOL_fg(:,2:end), theta);
f3 = tmp(1);
g3 = tmp(2);

return
end

