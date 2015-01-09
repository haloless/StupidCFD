% ## Copyright (C) 2014 homu
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

% ## LSNormGrad2D

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2014-04-27

function [ yy,max_radius,d_error,d_count ] = LSNormGrad2D (dd,d4,lo,hi,ng,dh, error_width)
% Description
% DD is the old LS, D4 is the new LS.

LSGlobals2D;

yy = zeros(size(dd));
d_error = 0.0;
d_count = 0;
max_radius = 1.0;

% d/dx, d2/dx2, d/dy, d2/dy2
[b00,b01,b10,b11] = LSDiffOp(d4,lo,hi,ng,dh);

for j = lo(2):hi(2)
for i = lo(1):hi(1)
    current_radius = 0;
    main_level = dd(i,j);
    
    for dim = 0:1 % d/dx & d/dy
        for iup = 0:1
            k1 = iup - 1;
            k1hold = -(2*k1 + 1);
            % 1st order
            if (dim == 0)
                dtable = b00(i+iup,j);
            else
                dtable = b10(i,j+iup);
            end
            % 2nd order
            if (dim == 0) 
                acmp = b01(i+iup-1,j);
                bcmp = b01(i+iup,j);
            else
                acmp = b11(i,j+iup-1);
                bcmp = b11(i,j+iup);
            end
            if (abs(acmp) < abs(bcmp))
                ddtable = acmp;
            else
                ddtable = bcmp;
            end
            
            dupwind = dtable + ddtable*k1hold;
            if (iup == 0)
                dminus = dupwind;
            else
                dplus = dupwind;
            end
        end % IUP
        
        dright = sign(main_level) * dplus;
        dleft = sign(main_level) * dminus;
        
        if (dleft>0 && dright<0) % compression
            current_radius = current_radius + max(dleft^2,dright^2);
        elseif (dleft>0) % upwind is left
            current_radius = current_radius + dleft^2;
        elseif (dright<0) % upwind is right
            current_radius = current_radius + dright^2;
        else % expansion
            current_radius = current_radius + ((dleft+dright)*0.5)^2;
        end
        
        if (dim == 0)
            dxplus = dplus; dxminus = dminus;
        else
            dyplus = dplus; dyminus = dminus;
        end
    end % DIM
    
    % norm of gradient LS
    current_radius = sqrt(current_radius);
    yy(i,j) = current_radius;
    
    current_level = d4(i,j);
    if (abs(current_level) <= error_width)
        % current cell is close to zero LS
        max_radius = max(max_radius, current_radius);
        d_count = d_count + 1;
        
        % with normalized error
        if (current_radius < eps_small)
            temp_error = 1.0 / eps_small;
        elseif (current_radius < 1.0)
            temp_error = 1.0 / current_radius;
        else
            temp_error = current_radius;
        end
        d_error = d_error + abs(1.0-temp_error);
    end
end
end

if (d_count > 0)
    d_error = d_error / d_count;
end

return
end % LSNormGrad2D

function [ b00,b01,b10,b11 ] = LSDiffOp(d4,lo,hi,ng,dh)
% dphi/dx
b00 = zeros(size(d4) + [1, 0]); 
% d2phi/dx2
b01 = zeros(size(d4));
% dphi/dy
b10 = zeros(size(d4) + [0, 1]);
% d2phi/dy2
b11 = zeros(size(d4));

for j = lo(2)-ng:hi(2)+ng
    for i = lo(1)-ng+1:hi(1)+ng
        b00(i,j) = (d4(i,j) - d4(i-1,j)) / dh;
    end
    for i = lo(1)-ng+1:hi(1)+ng-1
        b01(i,j) = 0.5 * (b00(i+1,j) - b00(i,j));
    end
end
for i = lo(1)-ng:hi(1)+ng
    for j = lo(2)-ng+1:hi(2)+ng
        b10(i,j) = (d4(i,j) - d4(i,j-1)) / dh;
    end
    for j = lo(2)-ng+1:hi(2)+ng-1
        b11(i,j) = 0.5 * (b10(i,j+1) - b10(i,j));
    end
end

return 
end % LSDiffOp



