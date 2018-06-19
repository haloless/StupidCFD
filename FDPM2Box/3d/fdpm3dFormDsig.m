function [ a ] = fdpm3dFormDsig(Be,sig,D,F,detF)
%fdpm3dFormDsig
% Form the spatial tangent operator [a]
% Note the [a]ijkl is encoded in 9x9.
% Input
% Be
%

% tensor func derivative
% L = d(log(b))/db
L = mathDerLogm(Be);


%% the following part depends on our tensor ordering
%% they can be confirmed by test/TestFormDsig

% B = dil*bjl + djk*bil
% note the order of B [11,22,33,12,21,23,32,31,13]
% Be is 3x3 form
B = zeros(9);
B(1,[1,4,9]) = 2 * Be([1,2,3]);
B(2,[2,5,6]) = 2 * Be([5,4,6]);
B(3,[3,7,8]) = 2 * Be([9,8,7]);
B(4,[1,2,4,5,6,9]) = Be([4,2,5,1,3,6]);
B(5,:) = B(4,:);
B(6,[2,3,5,6,7,8]) = Be([8,6,7,9,5,4]);
B(7,:) = B(6,:);
B(8,[1,3,4,7,8,9]) = Be([7,3,8,2,1,9]);
B(9,:) = B(8,:);

% Sijkl = sigma_il * djk
% note sigma is in voigt form, i.e. [11,22,33,12,23,31]
S = zeros(9);
S(1,[1,4,9]) = sig([1,4,6]);
S(2,[2,5,6]) = sig([2,4,5]);
S(3,[3,7,8]) = sig([3,5,6]);
S(4,[2,5,6]) = sig([4,1,6]);
S(5,[1,4,9]) = sig([4,2,5]);
S(6,[3,7,8]) = sig([5,2,4]);
S(7,[2,5,6]) = sig([5,6,3]);
S(8,[1,4,9]) = sig([6,5,3]);
S(9,[3,7,8]) = sig([6,4,1]);

D9 = [D(1:3,1:3), D(1:3,[4,4,5,5,6,6]); D([4,4,5,5,6,6],1:3), D([4,4,5,5,6,6],[4,4,5,5,6,6])];

a = 1/(2*detF) * D9 * L * B - S;







return
end


