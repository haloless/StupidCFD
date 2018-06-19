function [a] = fdpmFormAmat(mBe, vsig, D, detF, prob_type)

if prob_type == 1
    % plane-strain
    ncomp = 4;
    out_of_plane = false;
elseif prob_type == 2
    % axisymmetric
    ncomp = 5;
    out_of_plane = true;
end

L = fdpmCalcLmat(mBe, out_of_plane);

Lfull = zeros(ncomp);
Lfull(1,1:4) = L(1,[1,3,3,2]);
Lfull(2,1:4) = L(3,[1,3,3,2]);
Lfull(3,1:4) = L(3,[1,3,3,2]);
Lfull(4,1:4) = L(2,[1,3,3,2]);
if ncomp == 5
    Lfull(5,5) = L(4,4);
end

Dfull = zeros(ncomp);
if ncomp == 4
    Dfull = D([1,3,3,2],[1,3,3,2]);
elseif ncomp == 5
    Dfull = D([1,3,3,2,4],[1,3,3,2,4]);
end

B = zeros(ncomp);
B(1:4,1:4) = [...
2*mBe(1,1), 0, 2*mBe(1,2), 0;
mBe(2,1), mBe(1,1), mBe(2,2), mBe(1,2);
mBe(2,1), mBe(1,1), mBe(2,2), mBe(1,2);
0, 2*mBe(2,1), 0, 2*mBe(2,2)];
if ncomp == 5
    B(5,5) = 2*mBe(3,3);
end

S = zeros(ncomp);
S(1:4,1:4) = [...
vsig(1), 0, vsig(3), 0;
vsig(3), 0, vsig(2), 0;
0, vsig(1), 0, vsig(3);
0, vsig(3), 0, vsig(2)];
if ncomp == 5
    S(5,5) = vsig(4);
end

% 
a = 1/(2*detF) * Dfull * Lfull * B  - S;

return
end


