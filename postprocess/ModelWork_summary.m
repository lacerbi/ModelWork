function modelsummary = ModelWork_summary(mbag)
% MODELWORK_SUMMARY computes summary of model fits.
%
% MODELSUMMARY = MODELWORK_SUMMARY(MBAG) merges batches with identifiers in 
% list ID returning the merged model bag MBAG.
%

if nargin < 1
    help ModelWork_summary;
    return;
end
    
modellist = []; dataidlist = [];    cndlist = [];
modellen = 0;   dataidlen = 0;      cndlen = 0;
thetalen = 0;

for i = 1:length(mbag.bag)
    dataid = removetrailzeros(mbag.bag{i}.dataid);
    model = removetrailzeros(mbag.bag{i}.model);
    cnd = removetrailzeros(mbag.bag{i}.cnd);

    newmodel = 1;
    for j = 1:numel(modellist)
        if any(cellfun(@(x_) isequal(x_,model),modellist)); newmodel = 0; end
    end
    newdataid = 1;
    for j = 1:numel(dataidlist)
        if any(cellfun(@(x_) isequal(x_,dataid),dataidlist)); newdataid = 0; end
    end
    newcnd = 1;
    for j = 1:numel(cndlist)
        if any(cellfun(@(x_) isequal(x_,cnd),cndlist)); newcnd = 0; end
    end
    
    if newmodel
        modellist{end+1} = model;
        modellen = max(modellen,numel(model));
    end
    if newdataid
        dataidlist{end+1} = dataid; 
        dataidlen = max(dataidlen,numel(dataid));
    end
    if newcnd
        cndlist{end+1} = cnd;
        cndlen = max(cndlen,numel(cnd));
    end

    %thetalen = max(length(mbag.bag{i}.maptheta), thetalen);
end

dataidtab = zeros(numel(dataidlist),dataidlen);
for j = 1:numel(dataidlist)
    dataidtab(j,1:numel(dataidlist{j})) = dataidlist{j};
end
modeltab = zeros(numel(modellist),modellen);
for j = 1:numel(modellist)
    modeltab(j,1:numel(modellist{j})) = modellist{j};
end
cndtab = zeros(numel(cndlist),cndlen);
for j = 1:numel(cndlist)
    cndtab(j,1:numel(cndlist{j})) = cndlist{j};
end

% Sort conditions
dataidtab = sortrows(dataidtab);
cndtab = sortrows(cndtab);

modelsummary = ModelSummary(mbag, modeltab, dataidtab, cndtab, thetalen);



end

% Perform basic post-processing of model fit.
function modelfit = PostprocessModel(modelfit)
    
% Run some diagnostics on the sampled chains
nparams = length(modelfit.maptheta);
nstoredsamples = size(modelfit.smpl, 1);
nsamples = modelfit.nsamplesperchain;

if nstoredsamples > 0
    storedchainssize = floor(nstoredsamples/modelfit.nchains);    
    modelfit.kruskalwallisp = zeros(1, nparams);
    for i = 1:nparams
        smpl = zeros(storedchainssize, modelfit.nchains);
        for g = 1:modelfit.nchains
            smpl(:, g) = modelfit.smpl((g-1)*storedchainssize + (1:storedchainssize), i);         
            modelfit.kruskalwallisp(i) = kruskalwallis(smpl, [], 'off');    
        end
    end
end

end

function modelsummary = ModelSummary(mbag, models, dataid, cnd, thetalen)

% Exemplar model
m = mbag.bag{1};

% Convert full models to model numbers
modelsummary.project = mbag.project;
modelsummary.models = models;
modelsummary.dataid = dataid;
modelsummary.cnd = cnd;

nmodels = size(models,1);
ndata = size(dataid,1);
ncnd = size(cnd,1);

% Get model names
getModelNameFun = [modelsummary.project '_getModelName'];
if exist(getModelNameFun, 'file') && ~isempty(modelsummary.project)
    getModelNameFun = str2func(getModelNameFun);
else
    getModelNameFun = [];
end
for iModel = 1:nmodels
    if ~isempty(getModelNameFun)
        modelsummary.modelnames{iModel} = ...
            getModelNameFun(removetrailzeros(modelsummary.models(iModel,:)));
    else
        modelsummary.modelnames{iModel} = ['M' num2str(iModel)];
    end
end

% Create empty model metrics tables
modelsummary.nparams = NaN(1,nmodels);
metrics = {'aic','aicc','bic','dic','waic1','waic2','loocv','maploglike','marginallike','marginallike_rlr','marginallike_whmg','marginallike_whmu'};
for mcm = metrics
    modelsummary.(mcm{:}) = NaN(ndata,nmodels,ncnd);
end

%if isfield(m, 'robusttheta')
%    if ncnd == 1
%        modelsummary.robusttheta = NaN(max(dataid), size(models, 1), thetalen);        
%    else
%        modelsummary.robusttheta = NaN(max(dataid), size(models, 1), ncnd, thetalen);
%    end
%end

for j = 1:nmodels
    display(['Summarizing model ' num2str(j) ' out of ' num2str(nmodels) '.']);
    model = models(j, :);    
    for i = 1:ndata
        for k = 1:ncnd
            mfit = ModelBag_get(mbag, dataid(i, :), model, {cnd(k, :)});
            if ~isempty(mfit)
                m = mfit.metrics;
                
                if isnan(modelsummary.nparams(j))
                    modelsummary.nparams(j) = numel(mfit.maptheta);
                end
                
                for mcm = metrics
                    if isfield(m,mcm{:}) && ~isempty(m.(mcm{:}))
                        modelsummary.(mcm{:})(i, j, k) = m.(mcm{:});
                    end
                end
                
                s = mfit.sampling;
                if ~isempty(s) && ...
                        isfield(s,'sumstats') && ~isempty(s.sumstats)
                        if ~isfield(modelsummary,'sampling')
                            modelsummary.sampling.maxR = NaN(ndata,nmodels,ncnd);
                            modelsummary.sampling.minneff = NaN(ndata,nmodels,ncnd);
                            modelsummary.sampling.multiESS = NaN(ndata,nmodels,ncnd);
                            modelsummary.sampling.nchains = NaN(ndata,nmodels,ncnd);
                            modelsummary.sampling.psis_maxK = NaN(ndata,nmodels,ncnd);
                            modelsummary.sampling.psis_badKfreq = NaN(ndata,nmodels,ncnd);
                            modelsummary.sampling.psis_mehKfreq = NaN(ndata,nmodels,ncnd);
                        end
                        
                        if ~isempty(s.sumstats.R)
                            modelsummary.sampling.maxR(i, j, k) = max(s.sumstats.R);
                        end
                        if ~isempty(s.sumstats.neff)
                            modelsummary.sampling.minneff(i, j, k) = min(s.sumstats.neff);
                        end
                        if isfield(s.sumstats,'mESS_chain') && ~isempty(s.sumstats.mESS_chain)
                            modelsummary.sampling.minneff(i, j, k) = sum(s.sumstats.mESS_chain);
                        end
                        modelsummary.sampling.nchains(i, j, k) = s.nchains;
                        if isfield(s.sumstats,'ks') && ~isempty(s.sumstats.ks)
                            modelsummary.sampling.psis_maxK(i, j, k) = max(s.sumstats.ks);
                            modelsummary.sampling.psis_badKfreq(i, j, k) = sum(s.sumstats.ks > 1)/numel(s.sumstats.ks);                            
                            modelsummary.sampling.psis_mehKfreq(i, j, k) = sum(s.sumstats.ks > 0.5)/numel(s.sumstats.ks);                            
                        end
                end
                    
                %if isfield(m, 'robusttheta')
                %    if ncnd == 1
                %        modelsummary.robusttheta(dataid(i), j, 1:length(m.robusttheta)) = m.robusttheta;                        
                %    else
                %        modelsummary.robusttheta(dataid(i), j, keff, 1:length(m.robusttheta)) = m.robusttheta;
                %    end
                %end
            end
        end
    end
end

% Remove unused fields
for mcm = metrics
    if all(isnan(modelsummary.(mcm{:})(:)))
        modelsummary = rmfield(modelsummary, mcm{:});
    end
end

end