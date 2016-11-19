
function [ range ] = LevelNodeRange1D(ilevel)
% Node ID range on given level.

range = 2^(ilevel-1):(2^ilevel-1);

return
end

