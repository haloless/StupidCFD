
clear all;

a = 1;
H = a * 4;

is = 0:50;
dd = 5.5;
ss = 0.5 + 0.5*tanh(dd.*(is./max(is)-0.5)) ./ tanh(dd/2);

hsmall = a * 1.0e-3;
% remember to exclude the sphere radius
hmin = a + hsmall;
hmax = H - a - hsmall;
hs = hmin + (hmax-hmin).*ss;

c11 = [];
c12 = [];
c22 = [];
cs = [];

% hs = [1.01];
for h = hs
    
    P = InductionMat(a,h,H);
    
    Cinf = inv(P);
    
    C2 = C2b(a,h,H);
    
    C = Cinf + C2;
    
    cs(end+1,:) = reshape(C,1,[]);
end

figure;
plot(hs,cs(:,1),'x-', hs,cs(:,2),'x-', hs,cs(:,4),'x-');
legend('c11','c12','c22');





