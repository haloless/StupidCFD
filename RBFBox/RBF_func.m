% ## Copyright (C) 2013 homu
% ## 
% ## This program is free software; you can redistribute it and/or modify
% ## it under the terms of the GNU General Public License as published by
% ## the Free Software Foundation; either version 3 of the License, or
% ## (at your option) any later version.
% ## 
% ## This program is distributed in the hope that it will be useful,
% ## but WITHOUT ANY WARRANTY; without even the implied warranty of
% ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% ## GNU General Public License for more details.
% ## 
% ## You should have received a copy of the GNU General Public License
% ## along with Octave; see the file COPYING.  If not, see
% ## <http://www.gnu.org/licenses/>.

% ## RBF_func

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-24

function [ R dRdx dRdy ] = RBF_func (x,y,xs,ys,re)

% C^2k smoothness
% k = 1;
k = 2;

switch k
    case {1} % C2 smoothness
        rf = @(q) (1-q).^5 .* (8 + 40*q + 48*q.^2 + 25*q.^3 + 5*q.^4);
        drf = @(q) -5*(1-q).^4 .* (8 + 40*q + 48*q.^2 + 25*q.^3 + 5*q.^4) + ...
        (1-q).^5 .* (40 + 96*q + 75*q.^2 + 20*q.^3);
    case {2} % C4 smoothness
        rf = @(q) (1-q).^6 .* (6 + 36*q + 82*q.^2 + 72*q.^3 + 30*q.^4 + 5*q.^5);
        drf = @(q) -6*(1-q).^5 .* (6 + 36*q + 82*q.^2 + 72*q.^3 + 30*q.^4 + 5*q.^5) + ...
        (1-q).^6 .* (36 + 164*q + 216*q.^2 + 120*q.^3 + 25*q.^4);
    otherwise
        error('Unsupported RBF smoothness C(2k) with k=%d', k);
end

epsilon = re * 1e-4; % avoid zero division

rx = x - xs;
ry = y - ys;
rs = sqrt(rx.^2 + ry.^2);
qs = 1/re * rs;

%
R = (qs<=1) .* rf(qs);
%
dRdq = (qs<=1) .* drf(qs);
dRdx = 1/re * dRdq .* rx ./ (rs+epsilon);
dRdy = 1/re * dRdq .* ry ./ (rs+epsilon);

return
end




