%%
%% 
%%

clear;

sigma = 1.0;

thetaa = deg2rad(0);
thetab = deg2rad(60);

dist = linspace(0, 0.1, 11);

% 1. particle-particle force
R = 1;
V = 0.005;

datapp = [];
datapw = [];
for H = dist
    fppa = BridgeForceHR2(R,R,H,thetaa,thetaa,V,sigma);
    fppb = BridgeForceHR2(R,R,H,thetab,thetab,V,sigma);
    % fppb = BridgeForceHR2(R,R,H,thetaa,thetaa,V/2,sigma);
    
    fpwa = BridgeForceHR2(R,-1,H,thetaa,thetaa,V,sigma);
    fpwb = BridgeForceHR2(R,-1,H,thetab,thetab,V,sigma);
    
    datapp(end+1,:) = [fppa,fppb];
    
    datapw(end+1,:) = [fpwa,fpwb];
end


data0 = [];
% for H = dist
    % np = 51;
    
    % bridge = MakeBridge(R,R,H,thetaa,thetaa,V,sigma);
    % fppa = AxisymEvolverDirectForce(bridge,np,[],[]);
    
    % bridge = MakeBridge(R,R,H,thetab,thetab,V,sigma);
    % fppb = AxisymEvolverDirectForce(bridge,np,[],[]);
    
    % bridge = MakeBridge(R,-1,H,thetaa,thetaa,V,sigma);
    % fpwa = AxisymEvolverDirectForce(bridge,np,[],[]);
    
    % bridge = MakeBridge(R,-1,H,thetab,thetab,V,sigma);
    % fpwb = AxisymEvolverDirectForce(bridge,np,[],[]);
    
    
    % data0(end+1,:) = [fppa,fppb,fpwa,fpwb];
% end

data = [datapp,datapw,data0];





