function g = setflag(f,p,v,typeout)
%SETFLAG Set the value of a flag of a flags variable.
%
%   G = SETFLAG(F,P) sets the value of the P-th flag of flag variable F to 
%   true and returns the modified flags in the same format.
%   F can be a flags array in char or logical format, or a flags in
%   numerical format. 
%   If F is in numerical format, SETFLAG supports vectorized input/output.
%
%   G = SETFLAG(F,P,V) sets the value of the P-th flag of flag variable F 
%   to V. V needs to convert to a logical scalar (either true or false).
%
%   G = SETFLAG(F,P,V,TYPE) converts G to the specified TYPE: 
%      'c'  single char array
%      'l'  single logical array
%      'n'  numeric array
%   By default, G is the same type as F.
%
%   See also FLAGS2NUM, NUM2FLAGS.

if nargin<3 || isempty(v); v = true; end
if nargin<4 || isempty(typeout); typeout = []; end

% Find flag format
typein = [];
if isnumeric(f); typein = 'n';
elseif ischar(f); typein = 'c';
elseif islogical(f); typein = 'l';
end
if isempty(typein); error('F is not a flag variable.'); end
if isempty(typeout); typeout = typein; end

if typein == 'n' && typeout == 'n' && numel(f) > 1
    g = zeros(1,numel(f));
    for k = 1:numel(f); g(k) = setflag(f(k),p,v,'n'); end
    g = reshape(g,size(f));
    return;
elseif typein == 'n' && numel(f) > 1
    error('Single char/logical array output does not support vectorized input for F.');
end

% Convert to logical format
switch typein
    case 'n'; f = num2flags(f,'l');
    case {'c','s'}; f = logical(double(f-'0'));
    case 'l';
end

% Set flag
f(p) = logical(v);

% Convert back to original format, unless otherwise specified
switch typeout
    case 'n'; g = flags2num(f);
    case {'c','s'}; g = char(f+'0');
    case 'l';
end




