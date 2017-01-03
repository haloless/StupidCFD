
clear all;

a = 1;
% h = 3;
% h = 10.0677;
% h = 1.005;
% h = 1.0018;
h = 1.0002;

gap = h/a - 1;

alpha = acosh(h/a);

nmax = 500;
nmax2 = nmax + 2;

ks = zeros(nmax2,1);
for n = 0:nmax2-1
    ks(n+1) = (n+0.5) * coth((n+0.5)*alpha) - coth(alpha);
end

mat = sparse(nmax,nmax);
rhs = zeros(nmax,1);
for i = 1:nmax
    ca = (2*i-1)*ks(i) - (2*i-3)*ks(i+1);
    cb = (2*i+5)*ks(i+1) - (2*i+3)*ks(i+2);
    
    c0 = ca * (i-1)/(2*i-1);
    c1 = -ca*i/(2*i+1) - cb*(i+1)/(2*i+1);
    c2 = cb * (i+2)/(2*i+3);
    
    if i-1 >= 1
        mat(i,i-1) = c0;
    end
    mat(i,i) = c1;
    if i+1 <= nmax  
        mat(i,i+1) = c2;
    end
    
    % mat(i,i) = c0;
    % if i+1 <= nmax
        % mat(i,i+1) = c1;
    % end
    % if i+2 <= nmax
        % mat(i,i+2) = c2;
    % end
    
    a1 = (i+0.5) * alpha;
    a2 = (i-0.5) * alpha;
    a3 = (i+1.5) * alpha;
    rhs(i) = sqrt(2) * (2*coth(a1) - coth(a2) - coth(a3));
end

As = mat \ rhs;

%
fstar = 0;
for n = 0:nmax-1
    En = 2*sqrt(2)*exp(-(n+0.5)*alpha) / sinh((n+0.5)*alpha);
    if n >= 2
        En = En + ks(n+1) * (n-1)*n/(2*n-1) * As(n-1);
    end
    En = En - ks(n+1) * (n+1)*(n+2)/(2*n+3) * As(n+1);
    
    Cn = 0;
    if n >= 1
        if n >= 2
            Cn = Cn + (n-1)/(2*n-1)*As(n-1);
        end
        Cn = Cn - As(n) + (n+2)/(2*n+3)*As(n+1);
        Cn = -2*ks(n+1) * Cn;
    end
    
    fstar = fstar + En + n*(n+1)*Cn;
end
fstar = 1/6*sqrt(2) * sinh(alpha) * fstar;
fstar

fasymp = -8/15*log(gap) + 0.9543
