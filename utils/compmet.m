function [C,P] = compmet(varargin)
%COMPMET Compatibility metric between distributions.
%   C = COMPMET(X_1,...X_N) returns the minimum quantile interval, symmetric
%   around the median, at which two empirical distributions X and Y begin 
%   to overlap.

method = 'bayes';
TolFun = 1e-8;

if ischar(varargin{nargin}) && strcmpi(varargin{nargin},'pdf')
    pdfflag = 1;
    N = numel(varargin)-1;
else
    pdfflag = 0;
    N = numel(varargin);
end

ns = zeros(1,N);
lb = zeros(1,N);
ub = zeros(1,N);

if ~pdfflag
    % Distributions provided as samples; compute KDE
    for i = 1:N
        X{i} = sort(varargin{i});
        ns(i) = numel(X{i});
        lb(i) = X{i}(1);
        ub(i) = X{i}(end);
    end
    LB = min(lb);
    UB = max(ub);
    nx = 2^12;    
    
    b = zeros(1,N);
    for i = 1:N
        [b(i),pdf{i},xmesh,cdf{i}] = kde(X{i},nx,LB,UB);        
        pdf{i} = pdf{i}(:)';
    end
else
    % Pdfs already provided    
    nx = numel(varargin{1});
    for i = 1:N
        pdf{i} = varargin{i}(:)';
        if numel(pdf{i}) ~= nx
            error('All input PDFs should be defined on the same grid.');
        end
    end
    for i = 1:N
        cdf{i} = qcumtrapz(pdf{i});
        cdf{i} = cdf{i}/cdf{i}(end);
    end
end
    
% Compute overlap metric
if strcmpi(method,'joint')
    C = 0;
else
    C = zeros(N,nx);
end

% Compute joint posterior
Nx = numel(pdf{1});
joint_pdf = ones(1,Nx);
for i = 1:N; joint_pdf = joint_pdf .* pdf{i}; end

% Select 95% HPD of the joint posterior
[sortjoint,ordjoint] = sort(joint_pdf,'descend');
tmp = qcumtrapz(sortjoint);
tmp = tmp/tmp(end);
hpd = zeros(1,Nx);
hpd_idx = find(tmp > 0.95,1);
hpd(ordjoint(1:hpd_idx-1)) = 1;

H0 = ones(1,Nx);
H1 = 1;

Np = zeros(1,Nx);
pdfall = zeros(1,Nx);
for i = 1:N; pdfall = pdfall + pdf{i}/sum(pdf{i}); end
pdfall = cumsum(pdfall);
pdfall = pdfall / pdfall(end);
% for i = 1:N; Np = Np | (pdf{i} > TolFun); end
%Pleft = find(Np,1,'first');
%Pright = find(Np,1,'last');
Pleft = find(pdfall >= TolFun,1,'first');
Pright = find(pdfall <= 1-TolFun,1,'last');

Np = Pright - Pleft + 1;
% Np = Nx;
% Np / Nx

for i = 1:N
    
    switch lower(method)
        case 'joint'
            C(i) = sum(pdf{i}.*hpd)/sum(pdf{i});            
        case 'hpd'
            % [~,ord,rnk] = unique(pdf{i}); % Get pdf ranking (HPD regions)
            [sortpdf,ord] = sort(pdf{i},'descend');
            tmp = qcumtrapz(sortpdf);
            Q{i}(ord) = tmp / tmp(end);
            C(i,:) = 1 - Q{i};
        case 'quantile'
            cdf{i} = max(0,min(1,cdf{i}));
            C(i,:) = 1 - 2*abs(cdf{i} - 0.5);
        case 'bayes'
            pdf{i} = pdf{i}/sum(pdf{i});            
            H0 = H0 .* pdf{i};            
            H1 = H1 / Np;
    end
end

switch lower(method)
    case 'joint'
        C = prod(C).^(1/N);
    case 'bayes'
        H2 = 0;
        if N > 2 && 0
            for i = 1:N
                temp = ones(1,Nx);
                for j = 1:N
                    if i == j; continue; end
                    temp = temp .* pdf{j};                    
                end
                H2(i) = (1 / Np) .* sum(temp) / Np;
            end
        end
        
        H0 = sum(H0) / Np;
        Nf = H0 + H1 + sum(H2);
        H0 = H0 / Nf;
        H1 = H1 / Nf;
        H2 = H2 / Nf;
        C = H0;
        
        if all(H2 == 0); H2 = []; end
        P = [H0, H1, H2];
        
    otherwise
        C = max(prod(C,1)).^(1/N);
end

end

