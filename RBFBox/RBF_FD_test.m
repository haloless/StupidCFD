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

% ## RBF_FD_test

% ## Author: homu <homu@HOMU-PC>
% ## Created: 2013-07-26

clc;
clear all;

% Gaussian
RBF_GA = @(r,eps) exp(-(eps.*r).^2);
% inverse multiquadric
RBF_IMQ = @(r,eps) 1 ./ sqrt(1 + (eps.*r).^2);
% inverse quadratic
RBF_IQ = @(r,eps) 1 ./ (1+(eps.*r).^2);
% multiquadric, not positive definite
RBF_MQ = @(r,eps) sqrt(1 + (eps.*r).^2);

if (1) % 1D Hermite FD
    % xa = -2; xb = 2;
    % xa = 0; xb = sqrt(2*pi);
    xa = 0; xb = sqrt(4*pi);
    n = 128 + 1;
    h0 = (xb-xa) / (n-1);
    re = h0 * 14.5;
    % shape_param = h0 * 12.5;
    % shape_param = 2.5;
    shape_param = 1 / re;
    X = linspace(xa,xb,n)';
    % X = sort(rand(n,1)) * (xb-xa) + xa;
    
    %
    Z = X .* exp(-X.^2);
    dZdX = (1-2*X.^2) .* exp(-X.^2);
    %
    Z = sin(X.^2);
    dZdX = 2*X .* cos(X.^2);
    
    % Gaussian
    rbf = @(rs,e2) exp(-e2 .* rs.^2);
    Lrbf = @(rs,dxs,e2) -2*e2 .* exp(-e2.*rs.^2) .* dxs;
    L2rbf = @(rs,dxs,e2) -2*e2 .* exp(-e2*rs.^2) .* (1-2*e2*dxs.^2);
    % % MQ
    % rbf = @(rs,e2) sqrt(rs.^2 + e2);
    % Lrbf = @(rs,dxs,e2) dxs ./ sqrt(rs.^2+e2);
    % L2rbf = @(rs,dxs,e2) e2 ./ (rs.^2+e2).^(3/2);
    % % Wendland C2
    % rbf = @(rs,re) (rs<re) .* (1-rs/re).^4 .* (4*rs/re + 1);
    % Lrbf = @(rs,dxs,re) (rs<re) .* (-20 .* (1-rs/re).^3) .* 1/re^2 .* dxs;
    % L2rbf = @(rs,dxs,re) (rs<re) .* 20/re^2 .* (rs/re-1).^2 .* (4*rs/re-1);
    % Wendland C6
    WendlandC6 = @(q) (1-q).^8 .* (32*q.^3 + 25*q.^2 + 8*q + 1);
    WendlandC6dq = @(q) 22*q .* (q-1).^7 .* (16*q.^2 + 7*q + 1);
    WendlandC6dq2 = @(q) 22*(q-1).^6 .* (160*q.^3 + 15*q.^2 - 6*q - 1);
    rbf = @(rs,re) (rs<re) .* WendlandC6(rs/re);
    Lrbf = @(rs,dxs,re) (rs<re) .* WendlandC6dq(rs/re) .* (1/re) .* (dxs./(rs+1e-10));
    L2rbf = @(rs,dxs,re) (rs<re) .* WendlandC6dq2(rs/re) .* (1/re^2);
    
    R0 = zeros(n,n);
    LN = zeros(n,n);
    L2N = zeros(n,n);
    for i = 1:n
        xc = X(i);
        dxs = xc - X;
        rs = abs(dxs);
        
        % % GA
        % eps2 = shape_param^2;
        % phi = rbf(rs,eps2);
        % Lphi = Lrbf(rs,dxs,eps2);
        % L2phi = L2rbf(rs,dxs,eps2);
        
        % % MQ
        % eps2 = re^2;
        % phi = rbf(rs,eps2);
        % Lphi = Lrbf(rs,dxs,eps2);
        % L2phi = L2rbf(rs,dxs,eps2);
                
        % Wendland
        phi = rbf(rs,re);
        Lphi = Lrbf(rs,dxs,re);
        L2phi = L2rbf(rs,dxs,re);
        
        % phi = (rs<re) .* phi;
        % Lphi = (rs<re) .* Lphi;
        % L2phi = (rs<re) .* L2phi;
        
        R0(:,i) = phi;
        LN(:,i) = Lphi;
        L2N(:,i) = L2phi;
        
        % neigh = find(rs<=shape_param*40)';
        % R0(i,neigh) = phi(neigh)';
        % LN(neigh,i) = Lphi(neigh)';
    end
    
    p_order = 0;
    [P,Px] = MonomialBasis1D(p_order,X);
    m = MonomialBasisSize(1,p_order);
    
    % A = [R0', P'; P, zeros(m,m)];
    % B = [LN; zeros(m,n)];
    
    % A = [R0, LN', P';
        % -LN', L2N', zeros(n,m);
        % P, zeros(m,n), zeros(m,m)];
    A = [R0', LN', P'; 
    LN, L2N, zeros(n,m);
    P, zeros(m,n), zeros(m,m)];
    B = [LN; L2N; zeros(m,n)];
    
    sol = A \ B;
    %
    alpha = sol(1:n,1:n);
    %
    beta = sol(n+1:n+n,1:n);
    %
    mu = sol(end-m+1:end,1:n);
    cs = diag(mu' * Px);
    % dzdxh = alpha'*Z + cs;
    % dzdxh = alpha'*Z + beta'*dZdX + cs;
    dzdxh = (eye(n,n)-beta') \ (alpha'*Z + cs);
    
    if (1)
        figure;
        subplot(2,1,1);
        plot(X,Z, X,dZdX, X,dzdxh,'x'); legend('Z','dZ/dX','dZ/dX-RBF');
        subplot(2,1,2);
        plot(X, (dzdxh-dZdX)); title('dZ/dX-err');
    end
    
    abs(dzdxh((n-1)/2+1) - dZdX((n-1)/2+1))
    
end % end of 1D HFD

if (0) % augmented test
    
    refine = 2;
    n = 16*refine+1; 
    N = n^2;
    h0 = 4 / (n-1);
    dilation = 10.5;
    re = h0 * dilation;
    [X,Y] = ndgrid(linspace(-2,2,n),linspace(-2,2,n));
    
    % % exp
    % Z = X .* exp(-X.^2 - Y.^2);
    % dZdX = (1-2*X.^2) .* exp(-X.^2-Y.^2);
    % dZdY = -2*X.*Y .* exp(-X.^2-Y.^2);
    
    % % parabolic
    % Z = X.^2 + Y.^2;
    % dZdX = 2*X;
    % dZdY = 2*Y;
    
    % sin, cos
    kx = pi/4; ky = pi/4;
    Z = sin(kx*X) .* cos(ky*Y);
    dZdX = kx * cos(kx*X) .* cos(ky*Y);
    dZdY = -ky * sin(kx*X) .* sin(ky*Y);
    
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
    p_order = 0;
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
    if (1)
        figure;
        subplot(2,2,1); surf(X',Y',Z'); title('z');
        subplot(2,2,2); surf(X',Y',Zh'); title('z-RBF');
        subplot(2,2,3); surf(X',Y',(Zh-Z)'); title('z-err (must be 0)');
    end
    
    ngrow = ceil(dilation)+1;
    Is = ngrow+1:n-ngrow;
    Js = Is;
    
    dzdxh = Nx' * zs;
    dzdyh = Ny' * zs;
    % dzdxh = Nx' * [zs; zeros(M,1)]; dzdxh = dzdxh(1:N);
    % dzdyh = Ny' * [zs; zeros(M,1)]; dzdyh = dzdyh(1:N);
    % dZdXh = dZdX; dZdXh(Is,Js) = reshape(dzdxh,n,n)(Is,Js);
    % dZdYh = dZdY; dZdYh(Is,Js) = reshape(dzdyh,n,n)(Is,Js);
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














