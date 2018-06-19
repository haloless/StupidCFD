function [theta] = DerjaguinAngle(theta1,theta2)
% Average different contact angles

theta = acos(0.5*(cos(theta1)+cos(theta2)));

return
end


