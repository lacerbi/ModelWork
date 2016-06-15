function idx = findmodelhash(mbag,dataid,model,cnd)
% FINDMODELHASH return index of given model.
%
%   I = FINDMODELHASH(MBAG,DATAID,MODEL,CND) returns index of given model
%   in MBAG through hashing the fields DATAID, MODEL and CND.

h = modelhash(dataid,model,cnd,mbag.tablesize);
nnlist = mbag.table{h};    
for i = 1:size(nnlist, 2)
    idx = nnlist(i);
    if isfield(mbag.bag{idx},'dataid')
        dataid = mbag.bag{idx}.dataid;
    else    % Keep for retrocompatibility
        dataid = mbag.bag{idx}.nid;
    end    
    if all(dataid == dataid) && ...
            all(model == mbag.bag{idx}.model) && ...
            all(cnd == mbag.bag{idx}.cnd)
        return;
    end
end
idx = [];    % No match found