function x = rmtrailzeros(x)
%RMTRAILZEROS Remove trailing zeros from numerical array
%
%   X = RMTRAILZEROS(X) removes trailing zeros from numerical array X.

index = find(x~=0,1,'last');
if isempty(index); x = [];
elseif index < length(x); x(index+1:end) = []; end
