function [eta,etabar,xi] = materialDPSetAngle(friction_angle, dilation_angle, approx_type)

friction_angle = friction_angle / 180 * pi;
dilation_angle = dilation_angle / 180 * pi;


switch approx_type
case 'outer-cone'
    eta    = 6/sqrt(3) * sin(friction_angle) / (3-sin(friction_angle));
    etabar = 6/sqrt(3) * sin(dilation_angle) / (3-sin(dilation_angle));
    xi     = 6/sqrt(3) * cos(friction_angle) / (3-sin(friction_angle));
case 'inner-cone'
    eta    = 6/sqrt(3) * sin(friction_angle) / (3+sin(friction_angle));
    etabar = 6/sqrt(3) * sin(dilation_angle) / (3+sin(dilation_angle));
    xi     = 6/sqrt(3) * cos(friction_angle) / (3+sin(friction_angle));
case 'plain-strain'
    eta    = 3 * tan(friction_angle) / sqrt(9 + 12*tan(friction_angle)^2);
    etabar = 3 * tan(dilation_angle) / sqrt(9 + 12*tan(dilation_angle)^2);
    xi     = 3 / sqrt(9 + 12*tan(friction_angle)^2);
otherwise
    error('Unknown DP approx = %s', approx_type);
end


% outer approx.
% plain-strain approx.

return
end


