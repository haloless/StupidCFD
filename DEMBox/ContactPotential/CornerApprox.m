function [phi] = CornerApprox(xin,yin,dw,gw,hw,xbase)

x = xin - xbase(1);
y = yin - xbase(2);

g1 = gw(1);
g2 = gw(2);
h11 = hw(1,1);
h22 = hw(2,2);
h12 = hw(1,2);

phi = 0.5*h11*x.^2 + 0.5*h22*y.^2 + h12*x.*y + g1*x + g2*y + dw;

return
end

