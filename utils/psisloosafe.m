function loosafe = psisloosafe(loos,ks,ktol)
%PSISLOOSAFE Safety check for Pareto smoothed importance sampling

if nargin < 3 || isempty(ktol); ktol = 0.7; end

n=numel(loos);
loo=sum(loos);
idx = ks < ktol;
loosafe=sum(loos(idx))/sum(idx)*n;

loosafe - loo

end