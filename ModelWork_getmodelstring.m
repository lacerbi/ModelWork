function y = ModelWork_getmodelstring(x)
% MODELWORK_GETMODELSTRING Convert from model to string and vice versa.
%
%   MSTRING = MODELWORK_GETMODELSTRING(MODEL) returns the model string
%   MSTRING for model vector MODEL.
%
%   MODEL = MODELWORK_GETMODELSTRING(MSTRING) returns the model vector 
%   MODEL for model string MSTRING.
%
%   A model vector element can take integer values from 0 to 127.

% 28/10/2013 First version
% 20/08/2015 Complete revision

symbol = ['0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyzàáâãäåæçèéêëìíîïğñòóôõöøùúûüışÿ',char(256:290)];

if ischar(x)
    y = zeros(1,length(x));
    % Convert model string back to base number    
    for i = 1:length(x)
        index = find(x(i) == symbol,1);
        if isempty(index); error('Model number out of bounds.'); end
        y(i) = index - 1;        
    end
    y = rmtrailzeros(y);    % Remove trailing zeros
else
    x = rmtrailzeros(x);    % Remove trailing zeros
    y = zeros(1,length(x));    
    for i = 1:length(x); y(i) = symbol(x(i)+1); end
    y = char(y);
end