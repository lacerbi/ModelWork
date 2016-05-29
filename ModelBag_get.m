function [mfit,idx] = ModelBag_get(mbag,dataid,model,cnd)
%MODELBAG_GET get a model from a model bag.
%
%   MFIT = MODELBAG_GET(MBAG, DATAID, MODEL, CND) returns a model fit 
%   struct MFIT with data DATAID, model MODEL and condition CND taken from 
%   model bag MBAG. DATAID, MODEL and CND may specify multiple datasets, 
%   models or conditions, in which case MFIT is an array of structs.
%
%   [MFIT,I] = MODELBAG_GET(...) also returns the index I of the model 
%   struct in the model bag.

if nargin < 1
    help ModelBag_get;
    return;
end

if ~iscell(dataid); dataid = {dataid}; end
if ~iscell(cnd); cnd = {cnd}; end

if length(dataid) == 1 && size(model, 1) == 1 && length(cnd) == 1
    idx = findmodelhash(mbag,removetrailzeros(dataid{1}),removetrailzeros(model),removetrailzeros(cnd{1}));            
    if isempty(idx); mfit = [];
    else mfit = mbag.bag{idx}; end                
else
    mfit = [];                        
    for i = 1:length(dataid)
        for j = 1:size(model,1)                
            for k = 1:length(cnd)
                idx = findmodelhash(mbag, ...
                    removetrailzeros(dataid{i}), ...
                    removetrailzeros(model(j, :)), ...
                    removetrailzeros(cnd{k}) ...
                    );
                if isempty(idx); mfit{length(mfit)+1} = [];
                else mfit{length(mfit)+1} = mbag.bag{idx}; end
            end
        end
    end
end