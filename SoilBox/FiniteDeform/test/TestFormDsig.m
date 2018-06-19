
clear;


%
enc3 = [1,4,7;2,5,8;3,6,9];
enc3 = enc3';


%
enc9 = [1,4,9;5,2,6;8,7,3];
dec9 = zeros(9,2);
for i = 1:3
for j = 1:3
	dec9(enc9(i,j),:) = [i,j];
end
end

% 
enc6 = [1,4,6;4,2,5;6,5,3];


if 0
	% Bijkl = dik*bjl + djk*bil
	for i = 1:3
	for j = 1:3
		for k = 1:3
		for l = 1:3
			fprintf('B%d%d%d%d=',i,j,k,l);
			
			if i==k
				fprintf('b%d%d + ',j,l)
			end
			if j==k
				fprintf('b%d%d + ',i,l)
			end
			fprintf('\n')
		end
		end
	end
	end
end
if 1
	% Bijkl = dik*bjl + djk*bil
	% encode to 9x9
	for m = 1:9
		i = dec9(m,1);
		j = dec9(m,2);
		for n = 1:9
			k = dec9(n,1);
			l = dec9(n,2);
			ss = [];
			
			if i==k
				ss = [ss, sprintf('b%d + ',enc3(j,l))];
			end
			if j==k
				ss = [ss, sprintf('b%d + ',enc3(i,l))];
			end
			
			if ~isempty(ss)
				fprintf('B%d%d=%s\n',enc9(i,j),enc9(k,l),ss)
			end
		end
	end
end

if 1
	% Sijkl = sigma(il) * djk
	% sigma is in voigt form 6x1
	% encode to 9x9
	for m = 1:9
		i = dec9(m,1);
		j = dec9(m,2);
		for n =1:9
			k = dec9(n,1);
			l = dec9(n,2);
			ss = [];
			
			if j==k
				ss = [ss, sprintf('s%d',enc6(i,l))];
			end
			
			if ~isempty(ss)
				fprintf('S%d%d=%s\n',m,n,ss)
			end
		end
		% fprintf('\n')
	end
end
