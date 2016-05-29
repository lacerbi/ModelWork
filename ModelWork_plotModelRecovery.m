function [ctab,rmse] = ModelWork_plotModelRecovery(fakedata,mbag,modelsummary)
% MODELWORK_PLOTMODELRECOVERY Plot results of model recovery test.

% List of model comparison metrics
mcmlist = {'aicc','bic','marginallike','dic','waic1'};

% Reconstruct true model list
truemodelindex = zeros(1, length(fakedata));
truemodels = fakedata{1}.truemodel;
truemodelindex(1) = 1;
for iSubj = 2:length(fakedata)
    index = find(all(bsxfun(@eq, truemodels, fakedata{iSubj}.truemodel), 2), 1);
    if isempty(index); 
        truemodels = [truemodels; fakedata{iSubj}.truemodel];
        truemodelindex(iSubj) = size(truemodels, 1);
    else
        truemodelindex(iSubj) = index;
    end
end

% Build model recovery matrix for each model comparison metric
for imcm = 1:length(mcmlist)
    if ~isfield(modelsummary, mcmlist{imcm}) || isempty(modelsummary.(mcmlist{imcm}))
        ctab.(mcmlist{imcm}) = [];
        continue;
    end
    
    mcm = modelsummary.(mcmlist{imcm});
    % Convert to deviance
    if strcmpi(mcmlist{imcm}, 'marginallike') || strcmpi(mcmlist{imcm}, 'mlike')
        mcm = -2*mcm;
    end
    % Find best model per subject
    [~,best] = min(mcm,[],2);

    for iTruemodel = 1:size(truemodels, 1)
        for imodel = 1:size(modelsummary.models, 1)
            ctab.(mcmlist{imcm})(iTruemodel, imodel) = sum(best(truemodelindex == iTruemodel) == imodel);
        end
    end
end

% Build parameter recovery fractional RMSE
for iModel = 1:length(truemodels)
    model = truemodels(iModel,:);
    subjs = find(truemodelindex == iModel);
    truetheta = []; maptheta = [];
    for i = 1:length(subjs)
        mfit = ModelBag_get(mbag,fakedata{subjs(i)}.id,model,{modelsummary.cnd});
        if ~isempty(mfit)
            truetheta(i,:) = fakedata{subjs(i)}.truetheta;
            maptheta(i,:) = mfit.maptheta;
        end
    end
    for k = 1:size(maptheta,2)
        rmse(iModel,k) = sqrt(mean( ((truetheta(:,k)-maptheta(:,k))./maptheta(:,k)).^2)); 
    end
end

end