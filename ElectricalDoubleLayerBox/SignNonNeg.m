function [s] = SignNonNeg(x)

s = sign(x);
s(x==0) = 1;

return
end
