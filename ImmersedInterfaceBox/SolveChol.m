
function [x] = SolveChol(b,R,Rt,perm)

if isempty(Rt)
    Rt = R';
end

if isempty(perm)
    x = R \ (Rt \ b);
else
    x = zeros(size(b));
    x(perm) = R \ (Rt \ b(perm));
end

return
end

