function f = num2flags(x,type,n)
%NUM2FLAGS Convert a non-negative integer into a flag array.
%
%   F = NUM2FLAGS(X) converts a numeric array X into a logical flags array. 
%
%   F = NUM2FLAGS(X,TYPE) converts X to the specified representation for
%   the flags array, depending on TYPE: 
%      'c'  char
%      'l'  logical (default).
%
%   F = NUM2FLAGS(X,TYPE,N) returns a flags array of length at least N
%   (padding with appropriate zeros).
%
%   See also NUM2FLAGS.

if nargin<2 || isempty(type); type = 'l'; end
if nargin<3; n = []; end

if ~isnumeric(x); error('X needs to be a numeric array.'); end
if (x > 2^32 - 1); error('Integer is too large (maximum length of flags array is 32 bits).'); end

f = fliplr(dec2bin(x));
if ~isempty(n); f = [f, ones(1,round(n(1))-length(f))*'0']; end

switch lower(type(1))
    case {'c','s'};
    case 'l'; f = logical(double(f-'0'));
    otherwise; error('TYPE can be ''char'' or ''logical''.');
end