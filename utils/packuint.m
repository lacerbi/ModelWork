function c = packuint(x)
%PACKUINT Convert integer array to packed string representation.
%
%   C = PACKUINT(X) converts a numeric array of non-negative integers to
%   a 'packed' character array representation. The elements of X need to
%   be non-negative integers (from 0 to 1024).
%
%   See also UNPACKUINT.

% 21/08/2015 First version

if ~isnumeric(x); error('X needs to be a numeric array.'); end
    
if 0
    symbol = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ',char(256:290)];
    x = uint8(rmtrailzeros(x));    % Remove trailing zeros
    c = zeros(1,length(x));    
    for i = 1:length(x); c(i) = symbol(x(i)+1); end
    c = char(c);
else
    units = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    tens = ['abcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ'];
    base = numel(units);
    c = [];
    x = rmtrailzeros(x);    % Remove trailing zeros
    for i = 1:numel(x)
        if x(i) >= base
            c2 = mod(uint16(x(i)), base) + 1;
            c1 = floor(x(i) / base);
            c = [c,tens(c1),units(c2)];
        else
            c = [c,units(x(i)+1)];
        end
    end
    c = char(c);
end