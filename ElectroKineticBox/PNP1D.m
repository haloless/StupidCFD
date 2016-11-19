%
% Poisson-Nernst-Planck system
%


clear all;

% coord=0, cart; coord=1, polar
coord = 0;
% coord = 1;

xlo = 0.5;
xhi = 10.0;
Lx = xhi-xlo;
dh = 0.05;
nx = round(Lx / dh);
dx = Lx / nx;
dx2 = dx^2;


xcell = linspace(xlo+dx/2,xhi-dx/2,nx)';
xnode = linspace(xlo,xhi,nx+1)';

if coord == 0
    rcell = ones(nx,1);
    rnode = ones(nx+1,1);
elseif coord == 1
    rcell = xcell;
    rnode = xnode;
end

lambda = 2;
lambda2 = lambda^2;

nspec = 2;
specz(1) = 1.0;
specz(2) = -1.0;
specc = zeros(nx,nspec);
specc(:,1) = 1.0;
specc(:,2) = 1.0;
for ispec = 1:nspec
    specq(ispec) = sum(specc(:,ispec).*rcell);
end

psi = zeros(nx,1);

%% zeta-potential BC
bctype = [ 1, 1 ];
bcval = [ -1, 1 ];
%% surface charge BC
% bctype = [ 0, 1 ];
% bcval = [ -lambda2/2*(-0.5), 0 ];

disp(['Build Psi']);
Apsi = sparse(nx,nx);
rpsi = zeros(nx,1);
for i = 1:nx
    aa = 0;
    if i > 1
        cc = 1.0/dx2 * rnode(i);
        Apsi(i,i-1) = cc;
        aa = aa - cc;
    else
        if bctype(1) == 0
            rpsi(i) = rpsi(i) - bcval(1)*rnode(i)/dx;
        elseif bctype(1) == 1
            cc = 1.0/dx2 * rnode(i);
            aa = aa - cc*2;
            rpsi(i) = rpsi(i) + bcval(1)*cc*2;
        end
    end
    if i < nx
        cc = 1.0/dx2 * rnode(i+1);
        Apsi(i,i+1) = cc;
        aa = aa - cc;
    else
        if bctype(2) == 0
            rpsi(i) = rpsi(i) + bcval(2)*rnode(i+1)/dx;
        elseif bctype(2) == 1
            cc = 1.0/dx2 * rnode(i+1);
            aa = aa - cc*2;
            rpsi(i) = rpsi(i) + bcval(2)*cc*2;
        end
    end
    Apsi(i,i) = aa;
end

figdbg = figure;

% wspec = 0.1;
% wpsi = 0.1;
wspec = 0.2;
wpsi = 0.2;
% wpsi = 1.0 - wspec;

for iter = 1:40
    % charge density
    rho = zeros(nx,1);
    for ispec = 1:nspec
        rho = rho + specz(ispec).*specc(:,ispec);
    end
    
    bpsi = -lambda2/2 .* rho;
    bpsi = bpsi .* rcell - rpsi;
    sol = Apsi \ bpsi;
    psi = wpsi.*sol + (1-wpsi).*psi;
    
    %
    dpsi = zeros(nx+1,1);
    dpsi(2:nx) = 1.0/dx .* (psi(2:nx)-psi(1:nx-1));
    if bctype(1) == 0
        dpsi(1) = bcval(1);
    elseif bctype(1) == 1
        dpsi(1) = (psi(2)-bcval(1)) / (dx/2);
    end
    if bctype(2) == 0
        dpsi(nx+1) = bcval(2);
    elseif bctype(2) == 1
        dpsi(nx+1) = (bcval(2)-psi(nx)) / (dx/2);
    end
    
    for ispec = 1:nspec
        zz = specz(ispec);
        
        Aspec = sparse(nx,nx);
        bspec = zeros(nx,1);
        for i = 1:nx
            aa = 0;
            if i > 1
                Aspec(i,i-1) = rnode(i)*(1.0/dx2 - zz*dpsi(i)/(dx*2));
                aa = aa + rnode(i)*(-1.0/dx2 - zz*dpsi(i)/(dx*2));
            else
                % do nothing, nonflux BC
            end
            if i < nx
                Aspec(i,i+1) = rnode(i+1)*(1.0/dx2 + zz*dpsi(i+1)/(dx*2));
                aa = aa + rnode(i+1)*(-1.0/dx2 + zz*dpsi(i+1)/(dx*2));
            else
                % do nothing, nonflux BC
            end
            Aspec(i,i) = aa;
        end
        
        % NOTE Aspec is singular
        % use density conservation 
        qq = specq(ispec);
        Aspec(1,:) = rcell';
        bspec(1) = qq;
        
        cc = Aspec \ bspec;
        specc(:,ispec) = wspec.*cc + (1-wspec).*specc(:,ispec);
    end
    
    % figure;
    % plot(xcell,specc(:,1),'o-',xcell,specc(:,2),'x-');
    % legend('C1','C2');
    
    figure(figdbg)
    subplot(2,1,1);
    plot(xcell,specc(:,1),'o-',xcell,specc(:,2),'x-');
    legend('C1','C2');
    subplot(2,1,2);
    plot(xcell,psi,'x-');
    
    disp(['iter=',int2str(iter)]);
    
    pause;
end


figure;
plot(xcell,specc(:,1),'o-',xcell,specc(:,2),'x-');
% axis([0 2.5 0 3])
xlim([0 2.5])
legend('C1','C2');
figure;
plot(xcell,psi,'x-');
% axis([0 2.5 -1.5 0.0])
xlim([0 2.5])






















