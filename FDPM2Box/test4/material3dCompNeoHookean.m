function [] = material3dCompNeoHookean(F, lambda,mu)

Cstr = F.' * F;
Estr = 0.5*(Cstr - eye(3));

% J = det(F);

[W,S,D] = matconst(Cstr, lambda,mu);

dh = 1.0e-4;

Sdiff = zeros(3);
for i = 1:3
for j = 1:3
    cc = Cstr;
    cc(i,j) = cc(i,j) + dh;
    if i~=j; cc(j,i) = cc(j,i) + dh; end
    [wa] = matconst(cc, lambda,mu);
    
    cc = Cstr;
    cc(i,j) = cc(i,j) - dh;
    if i~=j; cc(j,i) = cc(j,i) - dh; end
    [wb] = matconst(cc, lambda,mu);
    
    Sdiff(i,j) = 2 * (wa-wb)/(dh*2);
    if i~=j; Sdiff(i,j) = 0.5 * Sdiff(i,j); end
end
end

[S, Sdiff]

Ddiff = zeros(9);
for i = 1:3
for j = 1:3
    cc = Cstr;
    cc(i,j) = cc(i,j) + dh;
    if i~=j; cc(j,i) = cc(j,i) + dh; end
    [~,sa] = matconst(cc, lambda,mu);
    
    cc = Cstr;
    cc(i,j) = cc(i,j) - dh;
    if i~=j; cc(j,i) = cc(j,i) - dh; end
    [~,sb] = matconst(cc, lambda,mu);
    
    Ddiff(:,flatindex(i,j)) = 2 * (sa(:)-sb(:))./(dh*2);
    if i~=j
        Ddiff(:,flatindex(i,j)) = 0.5 * Ddiff(:,flatindex(i,j));
    end
end
end

D
Ddiff

return
end

function [W,S,D] = matconst(Cstr, lambda,mu)

I = eye(3);

c1 = sum(diag(Cstr));
c2 = sum(Cstr(:).^2);
c3 = det(Cstr);

J = sqrt(c3);
logJ = log(J);

Cinv = inv(Cstr);

W = mu/2*(c1-3) - mu*logJ + lambda/2*logJ^2;

S = mu*(I-Cinv) + lambda*logJ*Cinv;

L = zeros(9);
for i = 1:3
for j = 1:3
    for k = 1:3
    for l = 1:3
        L(flatindex(i,j),flatindex(k,l)) = 0.5 * (Cinv(i,k)*Cinv(j,l) + Cinv(i,l)*Cinv(j,k));
    end
    end
end
end

D = lambda*Cinv(:)*Cinv(:)' + 2*(mu-lambda*logJ)*L;


return
end

function [ind] = flatindex(i,j)
ind = i + (j-1)*3;
return
end




