function x = modelhash(dataid,model,cnd,tablesize)
%MODELHASH Return hash number for given model.
%
%   X = MODELHASH(DATAID,MODEL,CND,TABLESIZE) returns hash number of given 
%   model with fields DATAID, MODEL and CND, for table size TABLESIZE.

key = uint64([dataid,model,cnd]);
x = uint64(5381);
for i = 1:size(key, 2); x = mod(x*33 + key(i), 2^32-1); end
x = mod(x, tablesize) + uint64(1);

