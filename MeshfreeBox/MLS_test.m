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

% ## MLS_test

% ## Author: admin_2 <admin_2@KOSHIZUKA>
% ## Created: 2013-07-19

clc;
clear all;



% random test
N = 64;
re = 0.25;
nodes = rand(N,2);
xs = nodes(:,1);
ys = nodes(:,2);

% xc = nodes(N/2,1);
% yc = nodes(N/2,2);
xc = 0.5;
yc = 0.5;

% figure; plot(xs,ys,'o',xc,yc,'x');

neigh = MLS_neigh(xc,yc,xs,ys,re);
[phi,phi_x,phi_y] = MLS_shape(xc,yc,xs,ys,re,neigh);

reproduce = [phi; phi_x; phi_y] * [ones(N,1),xs,ys]
theoretic = [1, xc, yc; 0, 1, 0; 0, 0, 1]


if (1) % test grid data
	clear all;
	
	n = 33;
	h0 = 4 / (n-1);
	dilation = 2.5;
	re = h0 * dilation;
	[X,Y] = ndgrid(linspace(-2,2,n),linspace(-2,2,n));
	Z = X .* exp(-X.^2 - Y.^2);
	% Zx = 
	% plot this
	figure; surf(X',Y',Z'); title('analytic');

	% convert to vector data
	N = n^2;
	xs = reshape(X,[],1);
	ys = reshape(Y,[],1);
	zs = reshape(Z,[],1);

	C = sparse(N,N); % connectivity matrix
	S = sparse(N,N); % shape matrix
	Dx = sparse(N,N); % gradient matrix
	Dy = sparse(N,N);
	for i = 1:N
		xc = xs(i);
		yc = ys(i);
		
		% connectivity (or neighborhood)
		neigh = MLS_neigh(xc,yc,xs,ys,re);
		C(i,:) = neigh;
		
		% shape function
		[phi,phi_x,phi_y] = MLS_shape(xc,yc,xs,ys,re,neigh);
		S(i,:) = phi;
		
		Dx(i,:) = phi_x;
		Dy(i,:) = phi_y;
	end
	
	if (1)
		disp('put in true nodal values space');
		MT = S';
		S = speye(N,N);
		Dx = (MT \ Dx')';
		Dy = (MT \ Dy')';
	end
	
	zh = S * zs;
	Zh = reshape(zh,n,n);
	figure; surf(X',Y',Zh'); title('MLS: Z');

	% dzh = D * zs;
	% dzdxh = reshape(dzh,2,[])(1,:);
	% dzdyh = reshape(dzh,2,[])(2,:);
	dzdxh = Dx * zs;
	dzdyh = Dy * zs;
	dZdXh = reshape(dzdxh,n,n);
	dZdYh = reshape(dzdyh,n,n);
	figure;
	subplot(2,2,1); surf(X',Y',dZdXh'); title('MLS: dZ/dX');
	subplot(2,2,2); surf(X',Y',dZdYh'); title('MLS: dZ/dY');
    
	dZdX = (1-2*X.^2) .* exp(-X.^2-Y.^2);
	subplot(2,2,3); surf(X',Y',dZdX'); title('Analytical: dZ/dX');
	dZdY = -2*X.*Y .* exp(-X.^2-Y.^2);
	subplot(2,2,4); surf(X',Y',dZdY'); title('Analytical: dZ/dY');
end

if (1) % test singular
    clear all;
    
    nx = 33;
    ny = 2;
    h0 = 4 / (nx-1);
    dilation = 2.5;
    re = h0 * dilation;
    
    [X,Y] = ndgrid(linspace(-2,2,nx),linspace(-0.5,0.5,ny));
    Z = X .* exp(-X.^2 - Y.^2);
    dZdX = (1-2*X.^2) .* exp(-X.^2-Y.^2);
    dZdY = -2*X.*Y .* exp(-X.^2-Y.^2);
    
    N = nx * ny;
    xs = X(:);
    ys = Y(:);
    zs = Z(:);
    
    S = sparse(N,N); % shape matrix
    Dx = sparse(N,N); % gradient matrix
    Dy = sparse(N,N);
    for i = 1:N
        xc = xs(i);
        yc = ys(i);
        
        % connectivity (or neighborhood)
        neigh = MLS_neigh(xc,yc,xs,ys,re);
        
        % shape function
        [phi,phi_x,phi_y] = MLS_shape(xc,yc,xs,ys,re,neigh);
        S(i,:) = phi;
        
        Dx(i,:) = phi_x;
        Dy(i,:) = phi_y;
    end
    
    % MLS approximation to Z
    zh = S * zs;
    Zh = reshape(zh,nx,ny);
    
    if (1)
        figure;
        subplot(2,2,1); surf(X',Y',Z'); title('z');
        subplot(2,2,2); surf(X',Y',Zh'); title('z-MLS');
        subplot(2,2,3); surf(X',Y',(Zh-Z)'); title('z-err');
    end
    
    %
    dzdxh = Dx * zs;
    dzdyh = Dy * zs;
    dZdXh = reshape(dzdxh,nx,ny);
    dZdYh = reshape(dzdyh,nx,ny);
    
    if(1)
        figure;
        subplot(2,2,1); surf(X',Y',dZdX'); title('dz/dx');
        subplot(2,2,2); surf(X',Y',dZdXh'); title('dz/dx-MLS');
        subplot(2,2,3); surf(X',Y',(dZdXh-dZdX)'); title('dz/dx-err');
        figure;
        subplot(2,2,1); surf(X',Y',dZdY'); title('dz/dy');
        subplot(2,2,2); surf(X',Y',dZdYh'); title('dz/dy-MLS');
        subplot(2,2,3); surf(X',Y',(dZdYh-dZdY)'); title('dz/dy-err');
    end
end




