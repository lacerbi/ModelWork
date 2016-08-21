function [C,Q] = compatibility(varargin)
%COMPATIBILITY Compatibility metric between distributions.
%   C = COMPATIBILITY(X_1,...X_N) returns the minimum quantile interval, symmetric
%   around the median, at which two empirical distributions X and Y begin 
%   to overlap.

N = numel(varargin);
nx = 2^12;

ns = zeros(1,N);
lb = zeros(1,N);
ub = zeros(1,N);

for i = 1:N
    X{i} = sort(varargin{i});
    ns(i) = numel(X{i});
    lb(i) = X{i}(1);
    ub(i) = X{i}(end);
end

LB = min(lb);
UB = max(ub);

% Compute KDE
C = zeros(N,nx);
b = zeros(1,N);
for i = 1:N
    [b(i),pdf{i},xmesh,cdf{i}] = kde(X{i},nx,LB,UB);    
    if 0
        cdf{i} = max(0,min(1,cdf{i}));
        C(i,:) = 1 - 2*abs(cdf{i} - 0.5);
    else
        % [~,ord,rnk] = unique(pdf{i}); % Get pdf ranking (HPD regions)
        [sortpdf,ord] = sort(pdf{i},'descend');
        
        tmp = qcumtrapz(sortpdf);
        Q{i}(ord) = tmp / tmp(end);
        
        %cc = 0;
        %for j = 1:nx
        %    cc = cc + pdf{i}(ord(j));
        %    Q{i}(ord(j)) = cc;
        %end
        % Q{i} = Q{i}/cc;
        C(i,:) = 1 - Q{i};
    end
end

C = max(prod(C,1)).^(1/N);

end

