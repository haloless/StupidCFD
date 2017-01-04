%% The asymptotic resistance equations for wall-sphere.
%% The resistance matrix in F=RU form
%% All normalized by 6*pi*mu
%%
%%
%% | Fxt   0     0     | 0     -Fxr  0     |
%% | 0     Fxt   0     | Fxr   0     0     |
%% | 0     0     Fzt   | 0     0     0     |
%% |-------------------|-------------------|
%% | 0     Tyt   0     | Tyr   0     0     |
%% | -Tyt  0     0     | 0     Tyr   0     |
%% | 0     0     0     | 0     0     Tzr   |
%%
%%


function [fxt,fxr,tyr,fzt,tzr] = WallAsymptotic(h)

s = h-1;
sinv = 1/s;
logs = log(s);

%
% parallel part
%
fxt = -8/15 * logs + 0.9543;
fxr = -2/15 * logs - 0.2573;
tyr = -8/15 * logs + 0.4945;

%
% perpendicular part
%

fzt = sinv - 0.2*logs + 0.9713;

% this is the Zeta-function as zeta(x=3)
zeta3 = 1.202056903159594;
tzr = 4/3*(zeta3 - 0.5*s*logs);


return
end

