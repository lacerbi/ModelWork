function c = packuint(x)
%PACKUINT Convert integer array to packed string representation.
%
%   C = PACKUINT(X) converts a numeric array of non-negative integers to
%   a 'packed' character array representation. The elements of X need to
%   be non-negative integers (from 0 to 127).
%
%   See also UNPACKUINT.

% 21/08/2015 First version

symbol = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ',char(256:290)];

if ~isnumeric(x); error('X needs to be a numeric array.'); end
    
x = uint8(rmtrailzeros(x));    % Remove trailing zeros
c = zeros(1,length(x));    
for i = 1:length(x); c(i) = symbol(x(i)+1); end
c = char(c);
