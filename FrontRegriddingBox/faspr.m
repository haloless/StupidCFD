
function [ aspr ] = faspr(ke,icp,pt)

kp1 = icp(ke,1);
kp2 = icp(ke,2);
kp3 = icp(ke,3);

p1 = pt(kp1,:);
p2 = pt(kp2,:);
p3 = pt(kp3,:);

s1 = norm(p2-p1);
s2 = norm(p3-p2);
s3 = norm(p1-p3);
s = (s1+s2+s3)/3.0;

at = norm(cross(p2-p1,p3-p1));

aspr = 0.5*s*s*sqrt(3)/at;

return
end


