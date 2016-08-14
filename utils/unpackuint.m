function x = unpackuint(c)
%UNPACKUINT Convert packed string representation into integer array.
%
%   X = UNPACKUINT(C) converts a 'packed' character array representation 
%   into a numeric array of non-negative integers (from 0 to 1024). 
%
%   See also PACKUINT.

%
%   MODEL = MODELWORK_GETMODELSTRING(MSTRING) returns the model vector 
%   MODEL for model string MSTRING.
%
%   A model vector element can take integer values from 0 to 127.

% 21/08/2015 First version

if ~ischar(c); error('C needs to be a character array.'); end

if 0
    symbol = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ',char(256:290)];
    x = zeros(1,length(c));
    for i = 1:length(c)
        index = find(c(i) == symbol,1);
        if isempty(index); error('Unknown character in character array representation C.'); end
        x(i) = index - 1;        
    end
    x = uint8(rmtrailzeros(x));    % Remove trailing zeros
else
    units = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'];
    tens = ['abcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ'];
    
    base = numel(units);
    x = [];
    tmp = 0;
    for i = 1:numel(c)
        index = find(c(i) == units,1);
        if isempty(index)
            index = find(c(i) == tens,1);
            if isempty(index) || tmp > 0
                error('Unknown character in character array representation C.');
            end
            tmp = base * index;
            i = i+1;
        else
            x = [x, index-1 + tmp];
            tmp = 0;
        end        
    end
    x = rmtrailzeros(double(x));
end