function [xc,nc] = concatchains(x1,n1,x2,n2,storedmax)
%CONCATCHAINS Concatenate observables from two MCMC chains.
%
%  [XA,NA] = CONCATCHAINS(X1,N1,X2,N2) concatenate observables X1 and X2 
%  from two MCMC chains of same or different lengths. N1 and N2 are the 
%  real sample lengths of each chain. XC is the matrix of concatenated 
%  observables and NA the real sample length of the concatenated chain. 
%  CONCATCHAINS respects the relative thinning factor of each chain.
%
%  [XA,NA] = CONCATCHAINS(X1,N1,X2,N2,STOREDMAX) limits the stored length 
%  of the concatenated chain at STOREMAX (default is no maximum length).

if nargin < 5 || isempty(storedmax); storedmax = Inf; end

len1 = size(x1,1);
len2 = size(x2,1);

if n1 < len1    
    error(['The real chain length (N1 = ' num2str(n1) ') cannot be lesser than the number of rows of X1 (' num2str(len1) ').']);
end
if n2 < len2    
    error(['The real chain length (N2 = ' num2str(n1) ') cannot be lesser than the number of rows of X2 (' num2str(len2) ').']);
end

if isempty(x1) && isempty(x2)
    xc = [];
    nc = 0;
    return;
end

% Compute thinning factor of each chain
thin1 = n1/len1;
thin2 = n2/len2;

% Put both chains at the higher thinning level
if thin1 >= thin2
    len2_thin = round(n2/thin1);
    idx = round(linspace(1,len2,len2_thin));
    x2 = x2(idx,:);
else
    len1_thin = round(n1/thin2);
    idx = round(linspace(1,len1,len1_thin));
    x1 = x1(idx,:);
end

% Concatenate chains
xc = [x1; x2];
nc = n1 + n2;

% Thin concatenated chain if needed
lena = size(xc,1);
if lena > storedmax
    idx = round(linspace(lena/storedmax,lena,storedmax));
    xc = xc(idx,:);
end

end