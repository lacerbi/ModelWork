%MODELWORK_GENMODELRECOVERY
function gendata = ModelWork_genModelRecovery(mbag,N,modelnames,seed)
%MODELWORK_GENMODELRECOVERY Generate fake datasets for model recovery.
%   GENDATA = MODELWORK_GENMODELRECOVERY(MBAG) generate fake datasets from
%   all models in model bag MBAG, ten datasets per subject.
%
%   GENDATA = MODELWORK_GENMODELRECOVERY(MBAG,N) generates N datasets per
%   subject.
%
%   GENDATA = MODELWORK_GENMODELRECOVERY(MBAG,N,MODELNAMES) only generates
%   datasets from models in MODELNAMES. MODELWORK_GENMODELRECOVERY throws
%   an error if a model in MODELNAMES is not in the provided MBAG, or if
%   there are multiple matches.
%
%   GENDATA = MODELWORK_GENMODELRECOVERY(MBAG,N,MODELNAMES,SEED) uses random
%   seed SEED for each dataset. Otherwise, the random seed is fixed by
%   a simple hash that involves model and data id.

% Number of fake datasets per subject
if nargin < 2 || isempty(N); N = 10; end

% Models (use all by default)
if nargin < 3; modelnames = []; end

% Random seed
if nargin < 4 || isempty(seed); seed = []; end

project = mbag.project;

% Fake data generation function and analytics function
gendataFun = str2func([project '_gendata']);
analyticsFun = str2func([project '_analytics']);    
gendata = []; % Cell array of fake data

modelsummary = ModelWork_summary(mbag);
if isempty(modelnames); modelnames = modelsummary.modelnames; end

for i = 1:numel(modelnames)
    fprintf('Generating fake datasets for model %s.\n', modelnames{i});
    
    idx = find(strcmp(modelsummary.modelnames,modelnames{i}));
    if isempty(idx)
        error(['Cannot find model ''' modelnames{i} ''' in provided model bag.']);
    elseif numel(idx) > 1
        error(['Multiple models matching the model name ''' modelnames{i} ''' in provided model bag.']);        
    end
    
    model = modelsummary.models(idx,:);
    mfit = ModelBag_get(mbag,modelsummary.dataid,model,modelsummary.cnd);
    
    fprintf('Generating fake datasets for dataid ');
    for j = 1:numel(mfit)
        
        % Fix random seed
        if isempty(seed)
            rnseed = prod(max(mfit{j}.dataid,1)) + prod(max(model,1));
        else
            rnseed = seed;
        end
        try RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',rnseed));
        catch; RandStream.setDefaultStream(RandStream.create('mt19937ar','seed',rnseed)); end
        
        fprintf('#%d', mfit{j}.dataid(1));
        if numel(mfit{j}.dataid > 1); fprintf('(%d)', mfit{j}.dataid(2)); end
        fprintf('...');
        
        % Create a number of raw datasets (data matrices)
        [genmats,trueparams] = gendataFun(N,mfit{j});
        D = analyticsFun(genmats);  % Analyze fake data
        for k = 1:N
            D{k}.dataid = mfit{j}.dataid;
            D{k}.fakeid = k;
            D{k}.truemodel = model;
            D{k}.trueparams = trueparams(k,:);
        end
        gendata = [gendata, D];
    end
    fprintf('\n');
end

end