
function [n1,n2,n3,nc1,nc2,nc3,kp] = flocal(ke,n1,icp,ine)

% NOTE the relationship between n1 and ine(n1)

kec = ine(ke,n1);

% [i,i+1,i+2]
n2 = mod(n1,3)+1;
n3 = mod(n2,3)+1;

kp = zeros(1,4);
%
kp(1) = icp(ke,n1);
kp(2) = icp(ke,n2);
kp(3) = icp(ke,n3);
%
nc1 = 0;
for i = 1:3
    if icp(kec,i) == kp(1)
        nc1 = i;
    end
end
if ~nc1
    fprintf('KE=%d,N1=%d,KP1=%d,KEC=%d\n', ke,n1,kp(1),kec);
    error('Node not match');
end
nc2 = mod(nc1,3)+1;
nc3 = mod(nc2,3)+1;
%
kp(4) = icp(kec,nc2);

if ine(kec,nc3)~=ke
    fprintf('KE=%d,N1=%d,KP1=%d,KEC=%d\n', ke,n1,kp(1),kec);
    error('Element not match');
end


return
end

