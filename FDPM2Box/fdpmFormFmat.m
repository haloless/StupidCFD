function [F3] = fdpmFormFmat(F2, prob_type, xref, xpos)
%fdpmFormFmat: Fill the (3,3) out-of-plane component of deform grad F.

error('bad');

F3 = F2;

if prob_type == 1 
    % plane strain
    % no deformation in Z-dir
    F3(3,3) = 1;
elseif prob_type == 2 
    % axisymmetric
    % the hoop deform = x/X
    if xref == 0 
        % for point just on axis
        F3(3,3) = 1;
    else
        F3(3,3) = xpos / xref;
    end
end


return
end
