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

% ## ../Utils/plot_mesh

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-30

function [  ] = plot_mesh (nodes,conn,elem)
% Description

numElems = size(conn,1);
numNodesPerElem = size(conn,2);


switch elem
    case {'Q4'}
        order = [1:numNodesPerElem, 1];
    otherwise
        error('Unsupported element type: %s', elem);
end

node_style = 'o';
mesh_style = '-';

xs = nodes(:,1);
ys = nodes(:,2);

% draw nodes
plot(xs,ys,node_style);
hold on;

% draw connectivity
for i = 1:numElems
    inodes = conn(i,order);
    plot(xs(inodes),ys(inodes),mesh_style);
end

% finish up
axis equal;
% axis off;
hold off;

return
end
