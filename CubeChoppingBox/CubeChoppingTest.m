
clc;
clear all;

% global
CubeChoppingGlobals;
ChopEps = 1.0e-8;

% nvec = [ 0, 1, 1 ];
nvec = [ 0, 0.0, 1 ];
% nvec = nvec + ChopEps;
nvec = nvec / norm(nvec,2);
% d = 0.25*sqrt(2); % 7/8
% d = 0; % 1/2
% d = sqrt(2); % 1

% nvec = [ 1, 1, 1 ];
% nvec = nvec / norm(nvec,2);
% d = -sqrt(3) / 3; % 47/48

cen = [ 0.5, 0.5, 0.5 ];
len = [ 1.0, 1.0, 1.0 ];

%
alpha0 = dot(nvec,cen);
% normalization denominator
denom0 = dot(nvec,len);

% find the max distance allowed
dmax = 0.5;
for kk = 0:1
for jj = 0:1
for ii = 0:1
    pt = [ ii, jj, kk ] .* len;
    dd = dot(pt-cen,nvec);
    dmax = max(dmax,abs(dd));
end
end
end

nsample = 20;

for nx = 0.0:0.05:1.0
for ny = nx:0.05:1.0
for nz = ny:0.05:1.0

nz = max(nz,ChopEps);
nvec = [ nx, ny, nz ];
nvec = nvec / norm(nvec,2);
nvec

%
alpha0 = dot(nvec,cen);
% normalization denominator
denom0 = dot(nvec,len);

% find the max distance allowed
dmax = 0.5;
for kk = 0:1
for jj = 0:1
for ii = 0:1
    pt = [ ii, jj, kk ] .* len;
    dd = dot(pt-cen,nvec);
    dmax = max(dmax,abs(dd));
end
end
end

% intercept -> volume -> intercept
% disp('intercept -> volume -> intercept');
ds = linspace(-dmax,dmax,nsample+1);
for d = ds
    m1 = nvec(1);
    m2 = nvec(2);
    m3 = nvec(3);
    alpha = d + alpha0;
    % normalize to STD case
    m1 = m1 / denom0;
    m2 = m2 / denom0;
    m3 = m3 / denom0;
    alpha = alpha / denom0;
    
    V = StdChopGetVolume(m1,m2,m3,alpha);
    a = StdChopGetIntercept(m1,m2,m3,V);
    
    % [ alpha, V, a ]
    
    % disp(['input=',num2str(alpha),', v=',num2str(V), ', a=',num2str(a)]);
    
    if (abs(a-alpha)>ChopEps) 
        Error(['input=',num2str(alpha),', v=',num2str(V), ', a=',num2str(a)]);
    end
end

% volume -> intercept -> volume
% disp('volume -> intercept -> volume');
vs = linspace(0.0, 1.0, nsample+1);
for vol = vs
    m1 = nvec(1);
    m2 = nvec(2);
    m3 = nvec(3);
    % normalize to STD case
    m1 = m1 / denom0;
    m2 = m2 / denom0;
    m3 = m3 / denom0;
    
    
    a = StdChopGetIntercept(m1,m2,m3,vol);
    V = StdChopGetVolume(m1,m2,m3,a);
    
    % [ vol, a, V ]
    
    % disp(['input=',num2str(alpha),', v=',num2str(V), ', a=',num2str(a)]);
    
    if (abs(V-vol)>ChopEps)
        Error(['input=',num2str(vol),', a=',num2str(a), ', v=',num2str(V) ])
    end
end

end
end
end

% m1 = nvec(1);
% m2 = nvec(2);
% m3 = nvec(3);
% alpha = d + dot(nvec,cen);
% % normalize to STD case
% denom = dot(nvec,len);
% m1 = m1 / denom;
% m2 = m2 / denom;
% m3 = m3 / denom;
% alpha = alpha / denom;

% V = StdChopGetVolume(m1,m2,m3,alpha)
% a = StdChopGetIntercept(m1,m2,m3,V)
