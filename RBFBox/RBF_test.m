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

% ## RBF_test

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-24

clc;
clear all;

% random test
N = 64;
re = 0.25;
nodes = rand(N,2);
xs = nodes(:,1);
ys = nodes(:,2);
xc = 0.5;
yc = 0.5;

neigh = meshfree_neigh(xc,yc,xs,ys,re);

if (0)
    figure; 
    plot(xs,ys,'o', xs(neigh),ys(neigh),'+', xc,yc,'x'); 
    axis equal;
    legend('xs','neigh','x');
    hold on;
    plot(xc+re*cos(linspace(0,2*pi,36)),yc+re*sin(linspace(0,2*pi,36)),'-');
end


if (0) % grid test
    clear all;
    
    refine = 2;
    n = 16*refine+1; N = n^2;
    h0 = 4 / (n-1);
    dilation = 2.5;
    re = h0 * dilation;
    [X,Y] = ndgrid(linspace(-2,2,n),linspace(-2,2,n));
    Z = X .* exp(-X.^2 - Y.^2);
    dZdX = (1-2*X.^2) .* exp(-X.^2-Y.^2);
    dZdY = -2*X.*Y .* exp(-X.^2-Y.^2);
    
    xs = reshape(X,[],1);
    ys = reshape(Y,[],1);
    zs = reshape(Z,[],1);
    
    % build RBF
    C0 = sparse(N,N);
    R0 = sparse(N,N);
    Ns = sparse(N,N); % shape func.
    Nx = sparse(N,N);
    Ny = sparse(N,N);
    for i = 1:N
        xc = xs(i);
        yc = ys(i);
        
        neigh = meshfree_neigh(xc,yc,xs,ys,re);
        % connectivity
        C0(i,neigh) = 1;
        
        [R,Rx,Ry] = RBF_func(xc,yc,xs(neigh),ys(neigh),re);
        R0(i,neigh) = R';
        
        % remember to transpose when use
        [Ns(:,i),Nx(:,i),Ny(:,i)] = RBF_spfunc(xc,yc,xs,ys,re,neigh);
    end
    
    Ns = R0' \ Ns;
    Nx = R0' \ Nx;
    Ny = R0' \ Ny;
    
    zh = Ns' * zs;
    Zh = reshape(zh,n,n);
    if (1)
        figure;
        subplot(2,2,1); surf(X',Y',Z'); title('z');
        subplot(2,2,2); surf(X',Y',Zh'); title('z-RBF');
        subplot(2,2,3); surf(X',Y',(Zh-Z)'); title('z-err (must be 0)');
    end
    
    dzdxh = Nx' * zs;
    dzdyh = Ny' * zs;
    dZdXh = reshape(dzdxh,n,n);
    dZdYh = reshape(dzdyh,n,n);
    if (1)
        figure;
        subplot(2,3,1); surf(X',Y',dZdX'); title('dz/dx');
        subplot(2,3,2); surf(X',Y',dZdXh'); title('dz/dx-RBF');
        subplot(2,3,3); surf(X',Y',(dZdXh-dZdX)'); title('dz/dx-err');
        subplot(2,3,4); surf(X',Y',dZdY'); title('dz/dy');
        subplot(2,3,5); surf(X',Y',dZdYh'); title('dz/dy-RBF');
        subplot(2,3,6); surf(X',Y',(dZdYh-dZdY)'); title('dz/dy-err');
    end

end % end of grid test

if (1) % augmented test
    clear all;
    
    refine = 1;
    n = 16*refine+1; 
    N = n^2;
    h0 = 4 / (n-1);
    dilation = 2.5;
    re = h0 * dilation;
    [X,Y] = ndgrid(linspace(-2,2,n),linspace(-2,2,n));
    
    Z = X .* exp(-X.^2 - Y.^2);
    dZdX = (1-2*X.^2) .* exp(-X.^2-Y.^2);
    dZdY = -2*X.*Y .* exp(-X.^2-Y.^2);
    
    Z = X.^2 + Y.^2;
    dZdX = 2*X;
    dZdY = 2*Y;
    
    xs = reshape(X,[],1);
    ys = reshape(Y,[],1);
    zs = reshape(Z,[],1);
    
    % build RBF-part
    C0 = sparse(N,N);
    R0 = sparse(N,N);
    Ns = sparse(N,N); % shape func.
    Nx = sparse(N,N);
    Ny = sparse(N,N);
    for i = 1:N
        xc = xs(i);
        yc = ys(i);
        
        neigh = meshfree_neigh(xc,yc,xs,ys,re);
        % connectivity
        C0(i,neigh) = 1;
        
        [R,Rx,Ry] = RBF_func(xc,yc,xs(neigh),ys(neigh),re);
        R0(i,neigh) = R';
        
        % remember to transpose when use
        [Ns(:,i),Nx(:,i),Ny(:,i)] = RBF_spfunc(xc,yc,xs,ys,re,neigh);
    end
    
    % build poly-part
    p_order = 2;
    [P0,P0x,P0y,M] = meshfree_MonoBasis(p_order,xs,ys);
    
    % augmented system
    A0 = [R0 P0'; P0 zeros(M,M)];
    
    Ns = A0' \ [Ns; P0];
    Nx = A0' \ [Nx; P0x];
    Ny = A0' \ [Ny; P0y];
    % shrink to square
    Ns = Ns(1:N,1:N);
    Nx = Nx(1:N,1:N);
    Ny = Ny(1:N,1:N);
    
    % Ns = R0' \ Ns;
    % Nx = R0' \ Nx;
    % Ny = R0' \ Ny;
    
    zh = Ns' * zs;
    % zh = Ns' * [zs; zeros(M,1)]; zh = zh(1:N);
    Zh = reshape(zh,n,n);
    if (0)
        figure;
        subplot(2,2,1); surf(X',Y',Z'); title('z');
        subplot(2,2,2); surf(X',Y',Zh'); title('z-RBF');
        subplot(2,2,3); surf(X',Y',(Zh-Z)'); title('z-err (must be 0)');
    end
    
    dzdxh = Nx' * zs;
    dzdyh = Ny' * zs;
    % dzdxh = Nx' * [zs; zeros(M,1)]; dzdxh = dzdxh(1:N);
    % dzdyh = Ny' * [zs; zeros(M,1)]; dzdyh = dzdyh(1:N);
    dZdXh = reshape(dzdxh,n,n);
    dZdYh = reshape(dzdyh,n,n);
    if (1)
        figure;
        subplot(2,3,1); surf(X',Y',dZdX'); title('dz/dx');
        subplot(2,3,2); surf(X',Y',dZdXh'); title('dz/dx-RBF');
        subplot(2,3,3); surf(X',Y',(dZdXh-dZdX)'); title('dz/dx-err');
        subplot(2,3,4); surf(X',Y',dZdY'); title('dz/dy');
        subplot(2,3,5); surf(X',Y',dZdYh'); title('dz/dy-RBF');
        subplot(2,3,6); surf(X',Y',(dZdYh-dZdY)'); title('dz/dy-err');
    end
    
    % %
    % coefs = A0 \ [zs; zeros(M,1)];
    % lambda = coefs(1:N);
    % beta = coefs(N+1:N+M);
    % %
    % zh = Ns'*lambda + P0'*beta;
    % Zh = reshape(zh,n,n);
    % if (0)
        % figure;
        % subplot(2,2,1); surf(X',Y',Z'); title('z');
        % subplot(2,2,2); surf(X',Y',Zh'); title('z-RBF');
        % subplot(2,2,3); surf(X',Y',(Zh-Z)'); title('z-err (must be 0)');
    % end
    % %
    % dzdxh = Nx'*lambda + P0x'*beta;
    % dzdyh = Ny'*lambda + P0y'*beta;
    % dZdXh = reshape(dzdxh,n,n);
    % dZdYh = reshape(dzdyh,n,n);
    % if (1)
        % figure;
        % subplot(2,3,1); surf(X',Y',dZdX'); title('dz/dx');
        % subplot(2,3,2); surf(X',Y',dZdXh'); title('dz/dx-RBF');
        % subplot(2,3,3); surf(X',Y',(dZdXh-dZdX)'); title('dz/dx-err');
        % subplot(2,3,4); surf(X',Y',dZdY'); title('dz/dy');
        % subplot(2,3,5); surf(X',Y',dZdYh'); title('dz/dy-RBF');
        % subplot(2,3,6); surf(X',Y',(dZdYh-dZdY)'); title('dz/dy-err');
    % end
    
end % end of augmented test



