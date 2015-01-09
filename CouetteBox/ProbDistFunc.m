
function [ dist ] = ProbDistFunc(x,y)

EBGlobals;

r = sqrt(x.^2 + y.^2);
dist = min(r-R0, R1-r);

return
end




