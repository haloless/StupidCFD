function [vec] = utilZipVec(varargin)
%utilZipVec
% Input
% [a1,a2,...,an], [b1,b2,...,bn], [c1,c2,...,cn]
% Output
% a column vector 
% [a1,b1,c1,a2,b2,c3,...,an,bn,cn]
%

if nargin <= 0
	error('Nothing to zip');
end

% number of arrays to zip
nvarargin = length(varargin);

% find the max length of inputs
nlen = 0;
for i = 1:nvarargin
	vin = varargin{i};
	if ~isempty(vin)
		nlen = max(nlen, numel(vin));
	end
end

if nlen <= 0
	error('Cannot zip zero length');
end

% allocate buffer
vec = zeros(nvarargin, nlen);

% stack input arrays in buffer
% empty [] will be pad by zero
% short array will be filled by zero
for i = 1:nvarargin
	vin = varargin{i};
	if ~isempty(vin)
		vec(i,1:numel(vin)) = vin(:);
	end
end

% ensure column vector
vec = reshape(vec, [],1);






return
end

