
function [ C ] = C2b(a,h,H)

hsmall = 0.1;

% low wall
xi = (h-a) / a;
b = sqrt(h^2/a^2-1);
mu = log(h/a+b);
if xi > hsmall
    c11 = b * Km(0,mu);
    c12 = b^2 * Km(1,mu);
    c22 = b^3 * Km(2,mu);
    C = 2 .* [c11,c12;c12,c22];
else
    C = Csmall(xi);
end

p11 = 1 - a/(2*h);
p12 = (a/(2*h))^2;
p22 = 1 - 2*(a/(2*h))^3;
P = [p11,p12;p12,p22];

Clo = C;
Plo = P;

% high wall
xi = (H-h-a) / a;
b = sqrt((H-h)^2/a^2-1);
mu = log((H-h)/a + b);
% mu = -log(h/a + b);
if xi > hsmall 
    c11 = b * Km(0,mu);
    c12 = -b^2 * Km(1,mu);
    c22 = b^3 * Km(2,mu);
    C = 2 .* [c11,c12;c12,c22];
else
    C = Csmall(xi);
    C(1,2) = -C(1,2);
    C(2,1) = C(1,2);
end

p11 = 1 - a/(2*(H-h));
p12 = -(a/(2*(H-h)))^2;
p22 = 1 - 2*(a/(2*(H-h)))^3;
P = [p11,p12;p12,p22];

Chi = C;
Phi = P;
% Phi = zeros(2,2);

% total C2b
C = Clo - inv(Plo);
C = C + Chi - inv(Phi);
% C = C + Chi;

return
end

function [ K ] = Km(m,mu)
    nlarge = 200;
    K = 0;
    for n = 0:nlarge
        kn = (2*n+1-coth(mu))^m / (exp((2*n+1)*mu)-1);
        K = K + kn;
    end
    return
end


function [ C ] = Csmall(xi)

euler_const = 0.5772156649;
zeta3 = 1.202056903159594;

f = 2*euler_const + log(2) - log(xi);

c11 = f + 1/3*(f+2/3)*xi;
c12 = 1/3*pi^2 - f + 1/9*(2*pi^2-5-12*f)*xi;
c22 = f - 2/3*pi^2 + 4*zeta3 + 1/9*(21*f+14-10*pi^2+36*zeta3)*xi;

C = 0.5 .* [c11,c12;c12,c22];
return
end

