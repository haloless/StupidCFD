% ## Copyright (C) 2013 admin_2
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

% ## MLS_basis

% ## Author: admin_2 <admin_2@KOSHIZUKA>
% ## Created: 2013-07-19

function [ p,dpdx,dpdy ] = MLS_basis (x,y,neigh)

% Description:
% return column-based representation
% [p(x1), p(x2), p(x3), ...], with p(x) = [1Å x y ...]^T

do_cutoff = true;
if nargin==2 % no neighborhood info., use full vector instead
	do_cutoff = false;
	neigh = ones(size(x));
end

if iscolumn(x); x=x'; end
if iscolumn(y); y=y'; end
if iscolumn(neigh); neigh=neigh'; end


MLS_order = 2;
k = (MLS_order+1)*(MLS_order+2)/2;

sz = size(x);
pad0 = zeros(sz);
pad1 = ones(sz);

switch MLS_order
	case {1}
		% 1st order, [1 x y]
		p = [ones(sz); x; y];
		dpdx = [zeros(sz); ones(sz); zeros(sz)];
		dpdy = [zeros(sz); zeros(sz); ones(sz)];
	case {2}
		% 2nd order, [1 x y x^2 xy y^2]
		p = [ones(sz); x; y; x.^2; x.*y; y.^2];
		dpdx = [zeros(sz); ones(sz); zeros(sz); 2*x; y; zeros(sz)];
		dpdy = [zeros(sz); zeros(sz); ones(sz); zeros(sz); x; 2*y];
	case {3}
		% 3rd order, [1 x y x^2 xy y^2 x^3 x2y xy2 y3]
		p = [ones(sz); x; y; x.^2; x.*y; y.^2; x.^3; (x.^2).*y; x.*(y.^2); y.^3];
		dpdx = [zeros(sz); ones(sz); zeros(sz); 2*x; y; zeros(sz); 3*x.^2; 2*x.*y; y.^2; zeros(sz)];
		dpdy = [zeros(sz); zeros(sz); ones(sz); zeros(sz); x; 2*y; zeros(sz); x.^2; 2*x.*y; 3*y.^2];
	otherwise
		error('Unsupported MLS_order=%d',MLS_order);
end

if do_cutoff
	for i = 1:rows(p)
		p(i,:) = neigh .* p(i,:);
		dpdx(i,:) = neigh .* dpdx(i,:);
		dpdy(i,:) = neigh .* dpdy(i,:);
	end
	p = sparse(p);
	dpdx = sparse(dpdx);
	dpdy = sparse(dpdy);
end

return
end
