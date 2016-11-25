
function [us,vs,ps] = FastSolve2d(kx,ky,fs,gs, nx,ny,dx,dy,nu)

ImUnit = 1.0i;

fh = fft2(fs);
gh = fft2(gs);

if 0
%
uh = zeros(nx,ny);
vh = zeros(nx,ny);
ph = zeros(nx,ny);

for j = 1:ny
for i = 1:nx
    ki = kx(i,j);
    kj = ky(i,j);
    % kn = k2(i,j);
	kn = ki^2 + kj^2;
	
    fhat = fh(i,j);
    ghat = gh(i,j);
    
    if kn == 0
        uh(i,j) = 0;
        vh(i,j) = 0;
        ph(i,j) = 0;
    else
        uh(i,j) = 1.0/kn / (kn*nu) * (kj*kj*fhat - ki*kj*ghat);
        vh(i,j) = 1.0/kn / (kn*nu) * (-ki*kj*fhat + ki*ki*ghat);
        ph(i,j) = 1.0/kn * (-ImUnit*ki*fhat - ImUnit*kj*ghat);
    end
end
end

else
	kxx = kx.^2;
	kyy = ky.^2;
	kxy = kx.*ky;
	
	k2 = kxx + kyy;
	k2(1,1) = 1;
	k4 = k2.*k2;
	
	uh = 1.0/nu ./ k4 .* (kyy.*fh - kxy.*gh);
	vh = 1.0/nu ./ k4 .* (kxx.*gh - kxy.*fh);
	ph = 1.0 ./ k2 .* (-ImUnit.*kx.*fh - ImUnit.*ky.*gh);
	
	uh(1,1) = 0;
	vh(1,1) = 0;
	ph(1,1) = 0;
end


%
us = real(ifft2(uh));
vs = real(ifft2(vh));
ps = real(ifft2(ph));


return
end

