
fid = fopen('tmp.dat','w');

xlo = 0.0;
xhi = 1.0;
ylo = 0.0;
yhi = 1.0;
nx = 20;
ny = 5;

npoint = (nx+ny)*2;
nelem = npoint;

xs = linspace(xlo,xhi,nx+1)';
ys = linspace(ylo,yhi,ny+1)';

xcoord = zeros(npoint,1);
ycoord = zeros(npoint,1);
I = 1:nx+1;
xcoord(I) = xs;
ycoord(I) = ylo;
I = nx+1:nx+ny+1;
xcoord(I) = xhi;
ycoord(I) = ys;
I = nx+ny+1:nx*2+ny+1;
xcoord(I) = flipud(xs);
ycoord(I) = yhi;
I = [nx*2+ny+1:nx*2+ny*2,1];
xcoord(I) = xlo;
ycoord(I) = flipud(ys);

elem = zeros(2,nelem);
bc = zeros(2,nelem);
I = 1:nx;
elem(1,I) = I; elem(2,I) = I+1;
bc(1,I) = 2; bc(2,I) = 0;
I = nx+1:nx+ny;
elem(1,I) = I; elem(2,I) = I+1;
bc(1,I) = 1; bc(2,I) = 1;
I = nx+ny+1:nx*2+ny;
elem(1,I) = I; elem(2,I) = I+1;
bc(1,I) = 2; bc(2,I) = 0;
I = nx*2+ny+1:nx*2+ny*2;
elem(1,I) = I; elem(2,I) = I+1; elem(2,nelem)=1;
bc(1,I) = 1; bc(2,I) = 1;

nfield = 41;
xfield = zeros(nfield,1);
yfield = zeros(nfield,1);
xfield(:) = linspace(xlo,xhi,nfield);
yfield(:) = 0.5*(ylo+yhi);
xfield(1) = xlo + 1.0e-3;
xfield(nfield) = xhi - 1.0e-3;


fprintf(fid,'square\n');
fprintf(fid,'%d %d\n', npoint,nfield);

fprintf(fid,'# Nodes:\n');
fprintf(fid,'%d %f %f\n',[1:npoint; xcoord'; ycoord']);

fprintf(fid,'# Elements and Boundary Conditions:\n');
fprintf(fid,'%d %d %d %d %f\n', [1:nelem; elem; bc]);

fprintf(fid,'# Field Points Inside Domain:\n');
fprintf(fid,'%d %f %f\n',[1:nfield; xfield'; yfield']);

fclose(fid);


