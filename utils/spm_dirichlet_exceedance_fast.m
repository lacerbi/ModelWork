function [xp,bor,pxp] = spm_dirichlet_exceedance_fast(alpha, Nsamp, lme, g, alpha0)
% Compute exceedance probabilities and related quantities for a Dirichlet distribution
% FORMAT xp = spm_dirichlet_exceedance_fast(alpha,Nsamp)
% 
% Input:
% alpha     - Dirichlet parameters
% Nsamp     - number of samples used to compute xp [default = 1e6]
% lme      - array of log model evidences 
%              rows: subjects
%              columns: models (1..Nk)
% g       - matrix of individual posterior probabilities
% alpha0   - [1 x Nk] vector of prior model counts
% 
% Output:
% xp        - exceedance probability
% bor       - Bayes Omnibus Risk (probability that model frequencies 
%             are equal)
% pxp       - protected exceedance probabilities
%__________________________________________________________________________
%
% This function computes exceedance probabilities, i.e. for any given model
% k1, the probability that it is more likely than any other model k2.  
% More formally, for k1=1..Nk and for all k2~=k1, it returns p(x_k1>x_k2) 
% given that p(x)=dirichlet(alpha).
% If requested (second and third outputs), it also computes the Bayes
% Omnibus Risk (probability that model frequencies are equal), and protected
% exceedance probabilities (probability that any one model is more 
% frequent than the others, above and beyond chance).
% 
% Refs:
% Stephan KE, Penny WD, Daunizeau J, Moran RJ, Friston KJ
% Bayesian Model Selection for Group Studies. NeuroImage (2008)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Based on code by Will Penny & Klaas Enno Stephan
% Modified by Luigi Acerbi 2016
% $Id: spm_dirichlet_exceedance.m 3118 2009-05-12 17:37:32Z guillaume $

if nargin < 2
    Nsamp = 1e6;
end

Nk = length(alpha);

% Compute exceedance probabilities p(r_i>r_j)
%--------------------------------------------------------------------------
if Nk == 2
    % comparison of 2 models
    xp(1) = spm_Bcdf_private(0.5,alpha(2),alpha(1));
    xp(2) = spm_Bcdf_private(0.5,alpha(1),alpha(2));
else
    % comparison of >2 models: use sampling approach

    % Perform sampling in blocks
    %--------------------------------------------------------------------------
    blk = ceil(Nsamp*Nk*8 / 2^28);
    blk = floor(Nsamp/blk * ones(1,blk));
    blk(end) = Nsamp - sum(blk(1:end-1));

    xp = zeros(1,Nk);
    for i=1:length(blk)

        % Sample from univariate gamma densities then normalise
        % (see Dirichlet entry in Wikipedia or Ferguson (1973) Ann. Stat. 1,
        % 209-230)
        %----------------------------------------------------------------------

        r = gamrnd(repmat(alpha,[blk(i),1]),1,blk(i),Nk);    
        r = bsxfun(@rdivide, r, sum(r,2));

        % Exceedance probabilities:
        % For any given model k1, compute the probability that it is more
        % likely than any other model k2~=k1
        %----------------------------------------------------------------------
        [y, idx] = max(r,[],2);
        xp = xp + histc(idx, 1:Nk)';

    end
    xp = xp / Nsamp;
end

% Compute Bayes Omnibus Risk - use functions from VBA toolbox
if nargout > 1
    posterior.a=alpha;
    posterior.r=g';
    priors.a=alpha0;

    F1 = FE(lme',posterior,priors); % Evidence of alternative
    
    options.families=[];
    F0 = FE_null(lme',options); % Evidence of null (equal model freqs)
    
    % Implied by Eq 5 (see also p39) in Rigoux et al.
    % See also, last equation in Appendix 2
    bor=1/(1+exp(F1-F0));

    % Compute protected exceedance probs - Eq 7 in Rigoux et al.
    if nargout > 2
        pxp=(1-bor)*xp+bor/Nk;
    end
end


function F = spm_Bcdf_private(x,v,w)
% Inverse Cumulative Distribution Function (CDF) of Beta distribution
% FORMAT F = spm_Bcdf(x,v,w)
%
% x   - Beta variates (Beta has range [0,1])
% v   - Shape parameter (v>0)
% w   - Shape parameter (w>0)
% F   - CDF of Beta distribution with shape parameters [v,w] at points x
%__________________________________________________________________________
%
% spm_Bcdf implements the Cumulative Distribution Function for Beta
% distributions.
%
% Definition:
%--------------------------------------------------------------------------
% The Beta distribution has two shape parameters, v and w, and is
% defined for v>0 & w>0 and for x in [0,1] (See Evans et al., Ch5).
% The Cumulative Distribution Function (CDF) F(x) is the probability
% that a realisation of a Beta random variable X has value less than
% x. F(x)=Pr{X<x}: This function is usually known as the incomplete Beta
% function. See Abramowitz & Stegun, 26.5; Press et al., Sec6.4 for
% definitions of the incomplete beta function.
%
% Variate relationships:
%--------------------------------------------------------------------------
% Many: See Evans et al., Ch5
%
% Algorithm:
%--------------------------------------------------------------------------
% Using MATLAB's implementation of the incomplete beta finction (betainc).
%
% References:
%--------------------------------------------------------------------------
% Evans M, Hastings N, Peacock B (1993)
%       "Statistical Distributions"
%        2nd Ed. Wiley, New York
%
% Abramowitz M, Stegun IA, (1964)
%       "Handbook of Mathematical Functions"
%        US Government Printing Office
%
% Press WH, Teukolsky SA, Vetterling AT, Flannery BP (1992)
%       "Numerical Recipes in C"
%        Cambridge
%__________________________________________________________________________
% Copyright (C) 1999-2011 Wellcome Trust Centre for Neuroimaging

% Andrew Holmes
% $Id: spm_Bcdf.m 4182 2011-02-01 12:29:09Z guillaume $


%-Format arguments, note & check sizes
%--------------------------------------------------------------------------
if nargin<3, error('Insufficient arguments'), end

ad = [ndims(x);ndims(v);ndims(w)];
rd = max(ad);
as = [[size(x),ones(1,rd-ad(1))];...
      [size(v),ones(1,rd-ad(2))];...
      [size(w),ones(1,rd-ad(3))]];
rs = max(as);
xa = prod(as,2)>1;
if sum(xa)>1 && any(any(diff(as(xa,:)),1))
    error('non-scalar args must match in size');
end

%-Computation
%--------------------------------------------------------------------------
%-Initialise result to zeros
F = zeros(rs);

%-Only defined for x in [0,1] & strictly positive v & w.
% Return NaN if undefined.
md = ( x>=0  &  x<=1  &  v>0  &  w>0 );
if any(~md(:))
    F(~md) = NaN;
    warning('Returning NaN for out of range arguments');
end

%-Special cases: F=1 when x=1
F(md & x==1) = 1;

%-Non-zero where defined & x>0, avoid special cases
Q  = find( md  &  x>0  &  x<1 );
if isempty(Q), return, end
if xa(1), Qx=Q; else Qx=1; end
if xa(2), Qv=Q; else Qv=1; end
if xa(3), Qw=Q; else Qw=1; end

%-Compute
F(Q) = betainc(x(Qx),v(Qv),w(Qw));


function [F,ELJ,Sqf,Sqm] = FE(L,posterior,priors)
% derives the free energy for the current approximate posterior
% This routine has been copied from the VBA_groupBMC function
% of the VBA toolbox http://code.google.com/p/mbb-vb-toolbox/ 
% and was written by Lionel Rigoux and J. Daunizeau
%
% See equation A.20 in Rigoux et al. (should be F1 on LHS)

[K,n] = size(L);
a0 = sum(posterior.a);
Elogr = psi(posterior.a) - psi(sum(posterior.a));
Sqf = sum(gammaln(posterior.a)) - gammaln(a0) - sum((posterior.a-1).*Elogr);
Sqm = 0;
for i=1:n
    Sqm = Sqm - sum(posterior.r(:,i).*log(posterior.r(:,i)+eps));
end
ELJ = gammaln(sum(priors.a)) - sum(gammaln(priors.a)) + sum((priors.a-1).*Elogr);
for i=1:n
    for k=1:K
        ELJ = ELJ + posterior.r(k,i).*(Elogr(k)+L(k,i));
    end
end
F = ELJ + Sqf + Sqm;



function [F0m,F0f] = FE_null(L,options)
% derives the free energy of the 'null' (H0: equal model frequencies)
% This routine has been copied from the VBA_groupBMC function
% of the VBA toolbox http://code.google.com/p/mbb-vb-toolbox/ 
% and was written by Lionel Rigoux and J. Daunizeau
%
% See Equation A.17 in Rigoux et al.

[K,n] = size(L);
if ~isempty(options.families)
    f0 = options.C*sum(options.C,1)'.^-1/size(options.C,2);
    F0f = 0;
else
    F0f = [];
end
F0m = 0;
for i=1:n
    tmp = L(:,i) - max(L(:,i));
    g = exp(tmp)./sum(exp(tmp));
    for k=1:K
        F0m = F0m + g(k).*(L(k,i)-log(K)-log(g(k)+eps));
        if ~isempty(options.families)
            F0f = F0f + g(k).*(L(k,i)-log(g(k))+log(f0(k)));
        end
    end
end
