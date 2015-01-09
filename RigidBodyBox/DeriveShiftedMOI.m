

function [ Jall ] = DeriveShiftedMOI(xcen,ycen,zcen,npart,xpart,ypart,zpart,mpart,Jpart)

xref = xpart - xcen;
yref = ypart - ycen;
zref = zpart - zcen;

J11 = (yref.^2 + zref.^2) .* mpart;
J12 = -xref .* yref .* mpart;
J13 = -xref .* zref .* mpart;
J21 = -yref .* xref .* mpart;
J22 = (zref.^2 + xref.^2) .* mpart;
J23 = -yref .* zref .* mpart;
J31 = -zref .* xref .* mpart;
J32 = -zref .* yref .* mpart;
J33 = (xref.^2 + yref.^2) .* mpart;

Jsft = [sum(J11),sum(J12),sum(J13); sum(J21),sum(J22),sum(J23); sum(J31),sum(J32),sum(J33)];
Jall = sum(Jpart,3) + Jsft;



return
end



