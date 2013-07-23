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

% ## MLS_weight

% ## Author: admin_2 <admin_2@KOSHIZUKA>
% ## Created: 2013-07-19

function [ w,dwdx,dwdy ] = MLS_weight (x,y,xs,ys,re,neigh)

if nargin == 5 % no neighborhood info., use a full vector instead
	neigh = ones(size(xs));
end

if isrow(xs); xs=xs'; end
if isrow(ys); ys=ys'; end
if isrow(neigh); neigh=neigh'; end

wf = @(q) neigh .* (q<=1) .* (1-6*q.^2+8*q.^3-3*q.^4);
dwf = @(q) neigh .* (q<=1) .* (-12*q + 24*q.^2 - 12*q.^3);

% avoid zero division
eps = 1e-6;

rx = x - xs;
ry = y - ys;
rs = sqrt(rx.^2 + ry.^2);
qs = 1/re * rs;

w = wf(qs);

dwdq = dwf(qs);
dwdx = 1/re * dwdq .* rx ./ (rs+eps);
dwdy = 1/re * dwdq .* ry ./ (rs+eps);

return
end


