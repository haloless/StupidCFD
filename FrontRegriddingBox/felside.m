
function [len,ind] = felside(ke,icp,pt)

% FTRegridGlobals;

kp1 = icp(ke,1);
kp2 = icp(ke,2);
kp3 = icp(ke,3);

% s1 = norm(pt(kp2,:)-pt(kp1,:));
% s2 = norm(pt(kp3,:)-pt(kp2,:));
% s3 = norm(pt(kp1,:)-pt(kp3,:));

len = [ norm(pt(kp2,:)-pt(kp1,:)),norm(pt(kp3,:)-pt(kp2,:)),norm(pt(kp1,:)-pt(kp3,:))];

[len,ind] = sort(len);


return
end


