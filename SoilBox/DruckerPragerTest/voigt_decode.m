function [afull] = voigt_decode(avoigt)

a1 = avoigt([1,6,5])';
a2 = avoigt([6,2,4])';
a3 = avoigt([5,4,3])';

afull = [ a1; a2; a3 ];

return
end

