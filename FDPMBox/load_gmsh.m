 %----------- Description ---------
% Transfer Mesh from GMSH to Matlab
%---------------------------------
% Written by 
% Shah, 18 November 2011

clc
clear all
close all

load Elements.txt
load Nodes.txt

nrelm=size(Elements,1);
nnod=size(Nodes,1);

Edof=zeros(nrelm,4);
Ex=zeros(nrelm,3); Ey=Ex; Ez=Ex;
Edof(:,1)=sort(1:nrelm);
Elm(:,1)=sort(1:nrelm);

for i=1:nrelm
    Edof(i,2:7)=[2*Elements(i,6)-1 Elements(i,6)*2 2*Elements(i,7)-1 ...
        Elements(i,7)*2 2*Elements(i,8)-1 2*Elements(i,8)];  
    Elm(i,2:4)=Elements(i,6:8);  
    Ex(i,:)=[Nodes(Elm(i,2),2) Nodes(Elm(i,3),2) Nodes(Elm(i,4),2)];
    Ey(i,:)=[Nodes(Elm(i,2),3) Nodes(Elm(i,3),3) Nodes(Elm(i,4),3)];
    Ez(i,:)=[Nodes(Elm(i,2),4) Nodes(Elm(i,3),4) Nodes(Elm(i,4),4)];    
end

nrdof=max(max(Edof));


figure(1);
hold on;
for i=1:length(Edof(:,1));
    plot(Ex(i,[1:end 1]),Ey(i,[1:end 1]));
end
axis equal


