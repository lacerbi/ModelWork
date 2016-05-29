function x = flags2num(f)
%FLAGS2NUM Convert a flag array into an integer.
%
%   X = FLAGS2NUM(F) converts a flag array F into a non-negative integer X. 
%   F can be a character array of '0' and '1' or a logical array. 
%   If F is a character array, characters different than '0' and '1' are 
%   ignored.
%
%   See also NUM2FLAGS.

if ischar(f)
    f(f ~= '0' & f ~= '1') = [];
    f(f == '0') = 0; f(f == '1') = 1;
    f = logical(double(f));
end
        
if islogical(f)
    x = sum((f(:)' ~= 0).*2.^(0:numel(f)-1));
    if (x > 2^32 - 1); error('Flag array is too long (maximum length 32 bits).'); end
    x = uint32(x);
else
    error('F needs to be either a string or a numerical or logical array.');
end