function ModelWork_plotParameterRecovery(fakedata,mbag,modelsummary)
% MODELWORK_PLOTPARAMETERRECOVERY Plot results of parameter recovery test.

if nargin < 3; modelsummary = []; end
if isempty(modelsummary); modelsummary = ModelWork_summary(mbag); end

models = modelsummary.models;

fakedata = fakedata(modelsummary.dataid(:,1)');
Ndatasets = numel(fakedata);

% Reconstruct subjects model list
modelindex = zeros(1, Ndatasets);
for iSubj = 1:Ndatasets
    idx = find(all(bsxfun(@eq, models, fakedata{iSubj}.truemodel), 2), 1);
    if isempty(idx)
        modelindex(iSubj) = NaN;
    else
        modelindex(iSubj) = idx;
    end
end

% Loop over models
for iModel = 1:size(models,1)    
    idx = find(modelindex == iModel);
    if isempty(idx); continue; end

    mfit = ModelBag_get(mbag,modelsummary.dataid((idx) - idx(1) + 1,:),models(iModel,:),modelsummary.cnd);
    Nfits = numel(mfit);
    if isempty(mfit{1}.mp)
        mfit{1} = ModelWork_loadFields(mbag.project,mfit{1});
    end
    Nparams = size(mfit{1}.maptheta,2);
    
    trueparams = NaN(Nfits, Nparams);
    fitparams = NaN(Nfits, Nparams);
    for i = 1:numel(mfit)
        trueparams(i,:) = fakedata{idx(i)}.trueparams;
        fitparams(i,:) = mfit{i}.maptheta;
    end
    
    nrows = 2;
    ncols = ceil(Nparams/2);
    
    figure;
    params = mfit{1}.mp.params;
    
    for i = 1:Nparams
        subplot(nrows,ncols,i);
        scatter(trueparams(:,i), fitparams(:,i)); hold on;
        axis square;
        box on;
        set(gca,'TickDir','out');
        p = params{i};
        p(p == '_') = '-';
        xlabel(['True ' p]);
        ylabel(['Fitted ' p]);
        plot(xlim,ylim,'k:','LineWidth',1);
        
        if i == floor(ncols/2)
            title(['Parameter recovery for model ' modelsummary.modelnames{iModel}]);
        end
    end
    
    set(gcf,'Color','w');
end


end
