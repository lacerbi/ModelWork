function [mfit,idx] = ModelBag_get(mbag,dataid,model,cnd)
%MODELBAG_GET get a model from a model bag.
%
%   MFIT = MODELBAG_GET(MBAG, DATAID, MODEL, CND) returns a model fit 
%   struct MFIT with data DATAID, model MODEL and condition CND taken from 
%   model bag MBAG. DATAID, MODEL and CND may specify multiple datasets, 
%   models or conditions, in which case MFIT is an array of structs.
%
%   MFIT = MODELBAG_GET(MBAG, MODEL) returns model fit MODEL where model
%   can be a model vector or model name, assuming default subjects and
%   conditions. If MODEL is a scalar, it picks the MODEL-th model in the
%   generated model summary.
%
%   [MFIT,I] = MODELBAG_GET(...) also returns the index I of the model 
%   struct in the model bag.

if nargin < 1
    help ModelBag_get;
    return;
end

if nargin < 3
    modelsummary = ModelWork_summary(mbag);
    model = dataid;
    dataid = modelsummary.dataid;
    cnd = modelsummary.cnd;
    
    if ischar(model)
        idx = find(strcmp(model,modelsummary.modelnames));
        if isempty(idx) || ~isscalar(idx)
            error(['Cannot find unique match for model ' model ' in model summary.']);
        end
        model = modelsummary.models(idx,:);
    elseif isscalar(model)
        model = modelsummary.models(model,:);        
    end
end

if ~iscell(dataid)  % Cellify
    temp = dataid;
    dataid = [];
    for i = 1:size(temp,1); dataid{i} = temp(i,:); end
end
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