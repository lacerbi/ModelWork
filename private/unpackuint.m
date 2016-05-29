function x = unpackuint(c)
%UNPACKUINT Convert packed string representation into integer array.
%
%   X = UNPACKUINT(C) converts a 'packed' character array representation 
%   into a numeric array of non-negative integers (from 0 to 127). 
%
%   See also PACKUINT.

%
%   MODEL = MODELWORK_GETMODELSTRING(MSTRING) returns the model vector 
%   MODEL for model string MSTRING.
%
%   A model vector element can take integer values from 0 to 127.

% 21/08/2015 First version

symbol = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ',char(256:290)];

if ~ischar(c); error('C needs to be a character array.'); end

x = zeros(1,length(c));
for i = 1:length(c)
    index = find(c(i) == symbol,1);
    if isempty(index); error('Unknown character is character array representation C.'); end
    x(i) = index - 1;        
end
x = uint8(rmtrailzeros(x));    % Remove trailing zeros