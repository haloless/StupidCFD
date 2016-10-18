
fid = fopen('tmp.dat','w');

ncirc = 2;
ndiv = [ 24,    48 ];
xcirc = [ 0.0,  0.0 ];
ycirc = [ 0.0,  0.0 ];
rcirc = [ 40,  80 ];
scirc = [ -1,   1 ];
ucirc = [ 1.5,  1.0 ];


npoint = sum(ndiv(1:ncirc));
nelem = npoint;

xcoord = zeros(npoint,1);
ycoord = zeros(npoint,1);
elem = zeros(2,nelem);
bc = zeros(2,nelem);

ioff = 0;
for icirc = 1:ncirc
    nn = ndiv(icirc);
    ind = 1:nn;
    I = ind + ioff;
    
    ang = scirc(icirc) * 2*pi/nn .* (ind-1);
    xx = xcirc(icirc);
    yy = ycirc(icirc);
    rr = rcirc(icirc);
    xcoord(I) = xx + rr.*cos(ang);
    ycoord(I) = yy + rr.*sin(ang);
    
    elem(1,I) = I;
    elem(2,I) = I+1; elem(2,I(end)) = I(1);
    
    bc(1,I) = 1;
    bc(2,I) = ucirc(icirc);
    
    ioff = ioff + ndiv(icirc);
end


nfield = 41;
xfield = zeros(nfield,1);
yfield = zeros(nfield,1);
xfield(:) = linspace(rcirc(1),rcirc(2),nfield);
yfield(:) = 0.0;
xfield(1) = rcirc(1) + 1.0e-3;
xfield(nfield) = rcirc(2) - 1.0e-3;


fprintf(fid,'circle\n');
fprintf(fid,'%d %d\n', npoint,nfield);

fprintf(fid,'# Nodes:\n');
fprintf(fid,'%d %f %f\n',[1:npoint; xcoord'; ycoord']);

fprintf(fid,'# Elements and Boundary Conditions:\n');
fprintf(fid,'%d %d %d %d %f\n', [1:nelem; elem; bc]);

fprintf(fid,'# Field Points Inside Domain:\n');
fprintf(fid,'%d %f %f\n',[1:nfield; xfield'; yfield']);

fclose(fid);


