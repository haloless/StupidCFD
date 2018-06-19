function [eigvec,eigdir,eigval] = voigt_eig(vm)

% convert to full matrix
m = voigt_decode(vm);

% eigen vectors and values
[v,d] = eig(m);


eigvec = v;
eigval = d([1,5,9])';

eigdir = zeros(6,3);
for i = 1:3
	vi = v(:,i);
	eigdir(:,i) = voigt_encode(vi * vi');
end


return
end


