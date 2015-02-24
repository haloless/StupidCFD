
clear all;

% global space
FTRegridGlobals;
%

% this contains reference cell length
load unit_sphere;
%
cen0 = [0.5,0.75,0.5]; 
rad0 = 0.15;

%
dh = reflen * rad0;
am = dh;
amax = 1.3 * dh;
amin = 0.4 * dh;
aspmax = 1.54;

%
maxpt = 4000;
maxel = 8000;
%
ptcon = zeros(maxpt,1);
bptcon = zeros(maxpt,1);
pt = zeros(maxpt,3);
elcon = zeros(maxel,1);
belcon = zeros(maxel,1);
icp = zeros(maxpt,3);
ine = zeros(maxel,3);

% load initial data
% transform unit sphere to problem setup
pt0 = pt0.*rad0 + repmat(cen0,size(pt0,1),1);
FTInit(pt0,tri0,nbr0);
clear reflen pt0 tri0 nbr0;






if (1)
    figure;
    trimesh(icp(1:nelem,:),pt(1:np,1),pt(1:np,2),pt(1:np,3));
    axis equal;
    axis([0 1 0 1 0 1]);
end

%
ProblemGlobals;
U = 1.0;
T = 4.0;

Tmax = T/2;
% Tmax = T;
max_step = 500;
dt = Tmax / max_step;
time = 0.0;
for step = 1:max_step
    
    t1 = time;
    [u1,v1,w1] = VortVel(pt(:,1),pt(:,2),pt(:,3), t1);
    pt1 = pt + 0.5*dt*[u1,v1,w1];
    
    t2 = time + dt*0.5;
    [u2,v2,w2] = VortVel(pt1(:,1),pt1(:,2),pt1(:,3), t2);
    
    pt = pt + dt*[u2,v2,w2];
    
    time = time + dt;
    
    if (1 && mod(step,5)==0)
        FTRegrid();
    end
    
    if mod(step,10)==0 || step==max_step
        prompt = ['step=',int2str(step),';time=',num2str(time),...
        ';np=',int2str(np),';nelem=',int2str(nelem)];
        disp(prompt);
        
        icp1 = icp;
        ke = ffe;
        for kk = 1:nelem
            icp1(kk,:) = icp(ke,:);
            ke = elcon(ke);
        end
        
        trimesh(icp1(1:nelem,:),pt(:,1),pt(:,2),pt(:,3));
        title(prompt);
        axis equal;
        axis([0 1 0 1 0 1]);
        
        view(0,90);
        drawnow
    end
end





