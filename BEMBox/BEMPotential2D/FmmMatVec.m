
function [ ax ] = FmmMatVec(u,ax,a,b)

BemMeshGlobals;
FmmTreeGlobals;

a = FmmUpward(u,a);
[ax,b] = FmmDownward(u,ax,a,b);

return
end

