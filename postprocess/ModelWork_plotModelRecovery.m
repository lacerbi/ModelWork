function [ctab,bias,se,mtab,stab] = ModelWork_plotModelRecovery(fakedata,mbag,modelsummary)
% MODELWORK_PLOTMODELRECOVERY Plot results of model recovery test.

% List of model comparison metrics
mcmlist = {'aicc','bic','marginallike','dic','waic1'};

Ndatasets = numel(fakedata);

% Reconstruct true model list
truemodelindex = zeros(1, Ndatasets);
truemodels = fakedata{1}.truemodel;
truemodelindex(1) = 1;
for iSubj = 2:Ndatasets
    index = find(all(bsxfun(@eq, truemodels, fakedata{iSubj}.truemodel), 2), 1);
    if isempty(index)
        truemodels = [truemodels; fakedata{iSubj}.truemodel];
        truemodelindex(iSubj) = size(truemodels, 1);
    else
        truemodelindex(iSubj) = index;
    end
end
Ntruemodels = size(truemodels,1);

% Rearrange model ordering
newmodels = [];
newmodelindex = NaN(1, Ndatasets);
neworder = NaN(Ntruemodels,1);
for i = 1:size(modelsummary.models,1)
    idx = find(all(bsxfun(@eq, modelsummary.models(i,:), truemodels),2),1);
    if isempty(idx); continue; end
    newmodels = [newmodels; modelsummary.models(i,:)];
    neworder(idx) = size(newmodels,1);
    newmodelindex(truemodelindex == idx) = neworder(idx);
end
    
neworder

if any(isnan(neworder))
    error('Some models were not assigned.');
end

truemodels = newmodels;
truemodelindex = newmodelindex;

% Build model recovery matrix for each model comparison metric
for imcm = 1:length(mcmlist)
    if ~isfield(modelsummary, mcmlist{imcm}) || isempty(modelsummary.(mcmlist{imcm}))
        ctab.(mcmlist{imcm}) = [];
        mtab.(mcmlist{imcm}) = [];
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
            mtab.(mcmlist{imcm})(iTruemodel, imodel) = mean(mcm(truemodelindex == iTruemodel,imodel));
            stab.(mcmlist{imcm})(iTruemodel, imodel) = stderr(mcm(truemodelindex == iTruemodel,imodel) - mcm(truemodelindex == iTruemodel,iTruemodel));
        end
    end
    
    % Remove diagonal
    mtab.(mcmlist{imcm}) = bsxfun(@minus, mtab.(mcmlist{imcm}), diag(mtab.(mcmlist{imcm})));
end


% Build parameter recovery fractional RMSE
for iModel = 1:size(truemodels,1)
    model = truemodels(iModel,:);
    subjs = find(truemodelindex == iModel);
    truetheta = []; maptheta = [];
    for i = 1:length(subjs)
        mfit = ModelBag_get(mbag,fakedata{subjs(i)}.dataid,model,{modelsummary.cnd});
        if ~isempty(mfit)
            if isfield(fakedata{subjs(i)},'truetheta')
                truetheta(i,:) = fakedata{subjs(i)}.truetheta;
            else
                truetheta(i,:) = fakedata{subjs(i)}.trueparams;                
            end
            maptheta(i,:) = mfit.maptheta;
        end
    end
    for k = 1:size(maptheta,2)
        % rmse(iModel,k) = sqrt(nanmean( ((truetheta(:,k)-maptheta(:,k))./maptheta(:,k)).^2)); 
        bias(iModel,k) = nanmean( truetheta(:,k)-maptheta(:,k) ); 
        se(iModel,k) = nanstd( truetheta(:,k)-maptheta(:,k) ) / sqrt(size(truetheta,1));
        % rmse(iModel,k) = sqrt(nanmean( (truetheta(:,k)-maptheta(:,k)).^2)); 
    end
end
bias(bias == 0) = NaN;
se(se == 0) = NaN;

end