% Output index to composite spatial modulus matrix
% [a] = 1/2J * [D][L][B] - [S]
%
%

clear;

mapsym = [ 1,3,0; 3,2,0; 0,0,4 ];
rmapsym = [ 1,1; 2,2; 1,2; 3,3 ];

mapfull = [ 1,3,0; 2,4,0; 0,0,5 ];
rmapfull = [ 1,1; 2,1; 1,2; 2,2; 3,3 ];

% B = dik*bjl + djk*bil
for p = 1:5
for q = 1:5
    i = rmapfull(p,1);
    j = rmapfull(p,2);
    k = rmapfull(q,1);
    l = rmapfull(q,2);
    fprintf('B%d%d=',p,q);
    
    if i == k
        fprintf('b%d + ', mapsym(j,l));
    end
    if j == k
        fprintf('b%d + ', mapsym(i,l));
    end
    fprintf('\n');
end
end

% S = sil*djk
for p = 1:5
for q = 1:5
    i = rmapfull(p,1);
    j = rmapfull(p,2);
    k = rmapfull(q,1);
    l = rmapfull(q,2);
    fprintf('S%d%d=',p,q);
    
    if j == k
        fprintf('s%d + ', mapsym(i,l));
    end
    fprintf('\n');    
end
end

