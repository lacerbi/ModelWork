function [mbag,modelsummary] = ModelWork_batchmerge(id,mbag,flag,flatten,prefix)
% MODELWORK_BATCHMERGE merges a batch of results of model analyses and
% also performs some postprocessing.
%
% MBAG = MODELWORK_BATCHMERGE(ID) merges batches with identifiers in 
% list ID returning the merged model bag MBAG.
%
% MBAG = MODELWORK_BATCHMERGE(ID,MBAG) merges batches with already 
% existing model bag MBAG.
%
% MBAG = MODELWORK_BATCHMERGE(ID,MBAG,FLAG) specifies whether to overwrite
% or not in case of overlap (by default FLAG is 1).
%
% MBAG = MODELWORK_BATCHMERGE(ID,MBAG,FLAG,FLATTEN) flattens
% the MODELSUMMARY table ignoring the experimental condition (all subjects
% count as having condition one). By default, FLATTEN is 0.
%
% MBAG = MODELWORK_BATCHMERGE(ID,MBAG,FLAG,FLATTEN,PREFIX) specifies the
% prefix (i.e., program name) to add to the model bag. If empty, use the
% one in the bag or leave it empty.
%
% [MBAG,MODELSUMMARY] = MODELWORK_BATCHMERGE(...) also returns a
% MODELSUMMARY table.
%

if nargin < 1
    help ModelWork_batchmerge;
    return;
end

if nargin < 2; mbag = []; end
if nargin < 3 || isempty(flag); flag = 1; end
if nargin < 4 || isempty(flatten); flatten = 0; end
if nargin < 5; prefix = []; end

bigbag = mbag;

for nid = id
    filename = ['outdata-' num2str(nid) '.mat'];
    display(['Processing file ''' filename '''.'])
    if ~exist(filename, 'file')
        error(['Cannot find file ''' filename '''.']);
    end
    
    % Load file and add its bag of models
    load(filename);    
    if isfield(mbag, 'prefix') && ~isempty(mbag.prefix); prefix = mbag.prefix; end
    
    for ii = 1:length(mbag.bag)
        mbag.bag{ii} = PostprocessModel(mbag.bag{ii});        
    end
        
    bigbag = ModelBag_add(bigbag, mbag, flag);
end

bigbag.prefix = prefix;
mbag = bigbag;
clear bigbag;

% Get model summary
if nargout > 1
    display('Computing model summary...');
    
    mlist = [];
    nidlist = [];
    cndlist = [];
    thetalen = 0;
    
    for i = 1:length(mbag.bag)
        model = mbag.bag{i}.model;
        nid = mbag.bag{i}.nid;
        cnd = mbag.bag{i}.cnd;
        
        newmodel = 1;
        for j = 1:size(mlist, 1)
            if all(model == mlist(j, :)); newmodel = 0; end
        end
        newnid = 1;
        for j = 1:length(nidlist)
            if nid == nidlist(j); newnid = 0; end
        end
        newcnd = 1;
        for j = 1:size(cndlist, 1)
            if all(cnd == cndlist(j, :)); newcnd = 0; end
        end
                
        if newmodel; mlist = [mlist; model]; end        
        if newnid; nidlist = [nidlist, nid]; end
        if newcnd; cndlist = [cndlist; cnd]; end
        
        thetalen = max(length(mbag.bag{i}.maptheta), thetalen);
    end
    
    % Sort conditions
    cndlist = sortrows(cndlist);
    
    modelsummary = ModelSummary(mbag, mlist, nidlist, cndlist, thetalen, flatten);
    
    % Convert model parameters on a session-by-session basis
    % (at the moment it works only for the best model)
    
    %bestmodel = find(mlist == 20561, 1);
    %if ~isempty(bestmodel)
    %    modelsummary = CoSMo13_convertParams(modelsummary, bestmodel);
    %end
end

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

function modelsummary = ModelSummary(mbag, models, nid, cnd, thetalen, flatten)

% Exemplar model
% m = ModelBag_get(mbag, nid(1), models(1, :), cnd(1, :));
m = mbag.bag{1};

% Convert full models to model numbers
modelsummary.prefix = mbag.prefix;
modelsummary.models = models;
modelsummary.nid = nid;
if flatten
    modelsummary.cnd = NaN;
    cndsize = 1;
else
    modelsummary.cnd = cnd;
    cndsize = size(cnd, 1);
end

% Get model names
getModelNameFun = [modelsummary.prefix '_getModelName'];
if exist(getModelNameFun, 'file') && ~isempty(modelsummary.prefix)
    getModelNameFun = str2func(getModelNameFun);
else
    getModelNameFun = [];
end
for i = 1:size(modelsummary.models, 1)
    if ~isempty(getModelNameFun)
        modelsummary.modelnames{i} = getModelNameFun(modelsummary.models(i,:));
    else
        modelsummary.modelnames{i} = ['M' num2str(i)];
    end
end

modelsummary.aic = NaN(max(nid), size(models, 1), cndsize);
if isfield(m, 'aicc')
    modelsummary.aicc = NaN(max(nid), size(models, 1), cndsize);
end
modelsummary.bic = NaN(max(nid), size(models, 1), cndsize);

if isfield(m, 'dic')
    modelsummary.dic = NaN(max(nid), size(models, 1), cndsize);
end
modelsummary.marginallike = NaN(max(nid), size(models, 1), cndsize);

if isfield(m, 'robusttheta')
    if cndsize == 1
        modelsummary.robusttheta = NaN(max(nid), size(models, 1), thetalen);        
    else
        modelsummary.robusttheta = NaN(max(nid), size(models, 1), cndsize, thetalen);
    end
end

for j = 1:size(models, 1)
    display(['Summarizing model ' num2str(j) ' out of ' num2str(size(models, 1)) '.']);
    model = models(j, :);
    for i = 1:length(nid)
        for k = 1:size(cnd, 1)
            if size(cnd, 2) > 1
                m = ModelBag_get(mbag, nid(i), model, {cnd(k, :)});
            else
                m = ModelBag_get(mbag, nid(i), model, cnd(k));                
            end
            if ~isempty(m)
                if flatten; keff = 1; else keff = k; end
                if ~isempty(m.aic); modelsummary.aic(nid(i), j, keff) = m.aic; end
                if isfield(m, 'aicc')
                    if ~isempty(m.aicc); modelsummary.aicc(nid(i), j, keff) = m.aicc; end
                end
                if ~isempty(m.bic); modelsummary.bic(nid(i), j, keff) = m.bic; end
                if isfield(m, 'dic') && ~isempty(m.dic)
                    modelsummary.dic(nid(i), j, keff) = m.dic;
                end
                if isfield(m, 'waic1') && ~isempty(m.waic1)
                    modelsummary.waic1(nid(i), j, keff) = m.waic1;
                end
                if isfield(m, 'waic2') && ~isempty(m.waic2)
                    modelsummary.waic2(nid(i), j, keff) = m.waic2;
                end
                if ~isempty(m.marginallike); modelsummary.marginallike(nid(i), j, keff) = m.marginallike; end
                if isfield(m, 'robusttheta')
                    if cndsize == 1
                        modelsummary.robusttheta(nid(i), j, 1:length(m.robusttheta)) = m.robusttheta;                        
                    else
                        modelsummary.robusttheta(nid(i), j, keff, 1:length(m.robusttheta)) = m.robusttheta;
                    end
                end
            end
        end
    end
end

end