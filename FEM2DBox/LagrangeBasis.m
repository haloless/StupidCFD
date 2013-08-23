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

% ## LagrangeBasis

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-30

function [ N dN ] = LagrangeBasis (elem_type,coord)
% Description
% return column vector

xi = coord(1);
eta = coord(2);

switch elem_type
    case {'Q4'}
        N = 1/4 * [(1-xi)*(1-eta); ...
                   (1+xi)*(1-eta); ...
                   (1+xi)*(1+eta); ...
                   (1-xi)*(1+eta)];
        dNdxi = 1/4 * [-(1-eta); ...
                       1-eta; ...
                       1+eta; ...
                       -(1+eta)];
        dNdeta = 1/4 * [-(1-xi); ...
                        -(1+xi); ...
                        1+xi; ...
                        1-xi];
        dN = [dNdxi, dNdeta];
    % case {'Q9'}
        % N = 1/4 * [1/4 * ];
    otherwise
        error('Unknown element type: %s', elem_type);
end

return
end
