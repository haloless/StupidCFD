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

% ## meshfree_MonoBasis

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-26

function [ p,px,py,m ] = meshfree_MonoBasis (order,x,y)

if iscolumn(x); x=x'; end
if iscolumn(y); y=y'; end

m = (order+1)*(order+2) / 2;

len = length(x);
e0 = zeros(1,len);
e1 = ones(1,len);

switch order
    case {0}
        % [1]
        p = e1;
        px = e0;
        py = e0;
    case {1}
        % [1 x y]
        p = [e1; x; y];
        px = [e0; e1; e0];
        py = [e0; e0; e1];
    case {2}
        % [1 x y x2 xy y2]
        p = [e1; x; y; x.^2; x.*y; y.^2];
        px = [e0; e1; e0; 2*x; y; e0];
        py = [e0; e0; e1; e0; x; 2*y];
    case {3}
        % [1 x y x2 xy y2 x3 x2y xy2 y3]
        p = [e1; x; y; x.^2; x.*y; y.^2; x.^3; (x.^2).*y; x.*(y.^2); y.^3];
        px = [e0; e1; e0; 2*x; y; e0; 3*x.^2; 2*x.*y; y.^2; e0];
        py = [e0; e0; e1; e0; x; 2*y; e0; x.^2; 2*x.*y; 3*y.^2];
    otherwise
        error('Unsupported order k=%d m=%d',order,m);
end


return
end




