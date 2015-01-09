
function [xpos,ypos,zpos] = UpdatePartPosition(xcen,ycen,zcen,omega,dt,npart,xpart,ypart,zpart)

xpos = xpart;
ypos = ypart;
zpos = zpart;

xref = xpart - xcen;
yref = ypart - ycen;
zref = zpart - zcen;

theta = norm(omega*dt);
ctheta = cos(theta);
stheta = sin(theta);

axis = omega / norm(omega);


for i = 1:npart
    rold = [xref(i); yref(i); zref(i)];
    proj = dot(rold, axis);
    orth = cross(axis, rold);
    
    rnew = ctheta*rold + proj*(1-ctheta)*axis + stheta*orth;
    xpos(i) = rnew(1) + xcen;
    ypos(i) = rnew(2) + ycen;
    zpos(i) = rnew(3) + zcen;
end

return
end

