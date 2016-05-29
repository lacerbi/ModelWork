function ModelWork_printParamsTable(mbag,modelsummary,modellist)
% MODELWORK_PRINTPARAMSTABLE LaTeX table of fitted parameters.

nid = modelsummary.nid;
cnd = modelsummary.cnd;

% Get model strings
for i = 1:size(modellist,1)
    % Find model in full model list
    model = modellist(i, :);
    pos = find(all(bsxfun(@eq, modelsummary.models, model),2), 1);
    
    if isfield(modelsummary, 'modelnames') && ~isempty(modelsummary.modelnames)
        modelstr{i} = modelsummary.modelnames{pos};
    else
        modelstr{i} = ['M' num2str(pos)];
    end 
    
    mfit = ModelBag_get(mbag, nid, model, {cnd});
    params = mfit{1}.mp.params;
    
    for j = 1:length(mfit)
        for k = 1:length(params)
            theta(j, k) = mfit{j}.mp.fulltheta{1}.(params{k});
        end
    end
    
    fprintf(' ');
    for k = 1:length(params); fprintf(' & $%s$ ', params{k}); end
    fprintf('\\\\\n');
    fprintf('%s ', modelstr{i});
    for k = 1:length(params)
        fprintf('& $%.2f \\pm %.2f$ ', nanmean(theta(:,k),1), stderr(theta(:,k),[],1));
    end
    fprintf('\\\\\n');
        
end

end