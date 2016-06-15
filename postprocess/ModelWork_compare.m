function [ranktable scoretable] = ModelWork_compare(modelsummary,cnd,metric,models)
% MODELWORK_COMPARE analyze results from model comparison.
%
% RANKTABLE = MODELWORK_COMPARE(MODELSUMMARY,CND,METRIC) compare models in 
% summarized model table MODELSUMMARY for condition CND under metric METRIC.
% Each element of RANKTABLE corresponds to the ranking of a certain model
% (column) for a given subject (row). Models are ordered according to
% MODELSUMMARY.MODELS.
%
% RANKTABLE = MODELWORK_COMPARE(MODELSUMMARY,CND,METRIC,MODELS) only use
% models in MODELS list.
%
% [RANKTABLE,SCORETABLE] = MODELWORK_COMPARE(MODELSUMMARY,CND,METRIC,MODELS) 
% also return a table of performance scores SCORETABLE (in which models not 
% in the list MODELS are set to NaN).
%
% Common usage:
% [ranktable,scoretable] = ModelWork_compare(modelsummary,1,'aicc',models)

if nargin < 1
    help ModelWork_compare;
    return;
end

if nargin < 4; models = []; end

table = modelsummary.(metric);

[score, ord] = sort(squeeze(table(:, :, cnd)), 2, 'ascend'); 

nid = modelsummary.nid;
if isempty(models); models = modelsummary.models; end

% Rank models

ranktable = NaN(length(nid), size(modelsummary.models, 1));

for i = 1:length(nid)
    for j = 1:size(models, 1)
        mnumber = findmodelnumber(modelsummary.models, models(j, :));
        if ~isempty(mnumber)
            p = find(ord(nid(i), :) == mnumber, 1);
            ranktable(i, mnumber) = p;
        end
    end    
end


% Adjust rankings if there are missing models
if any(isnan(ranktable(:)))
    for i = 1:length(nid)
        [~, rank] = sort(ranktable(i, :));
        notnans = sum(~isnan(ranktable(i, :)));
        ranktable(i, rank(1:notnans)) = 1:notnans; 
    end
end

% Compile scoretable
if nargout > 1
    scoretable = NaN(length(nid), size(modelsummary.models, 1));    
    for i = 1:length(nid)
        for j = 1:size(models, 1)
            mnumber = findmodelnumber(modelsummary.models, models(j, :));
            if ~isempty(mnumber)
                scoretable(i, mnumber) = table(i, mnumber);
            end
        end    
    end
end

end


% FINDMODEL return index of given model (private function)
function nn = findmodelnumber(models, model)

    nn = 1:size(models, 1);
    for i = 1:length(model)
        if isempty(nn); break; end
        subnn = find(models(nn, i) == model(i));
        nn = nn(subnn);        
    end

    % Multiple matching models?
    if length(nn) > 1
        error('Duplicate model in model list?');
    end
end
