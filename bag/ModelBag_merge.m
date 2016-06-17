function mfit = ModelBag_merge(project,m1,m2,data)
% MODELBAG_MERGE merges two model fit structures.
%
%   MFIT = MODELBAG_MERGE(PROJECT,M1,M2) merge the model fits M1 and M2 for 
%   project PROJECT and returns the merged MFIT.
%
%   MFIT = MODELBAG_MERGE(PROJECT,M1,M2,DATA) uses dataset DATA to update
%   model bag, if necessary.

if nargin < 4; data = []; end

% The following elements need to match
equalfields = {'type','dataid','model','cnd','nData'};
for i = 1:length(equalfields)
    f1 = m1.(equalfields{i});
    f2 = m2.(equalfields{i});
    if ~all(f1(:) == f2(:)); error(['Cannot merge: Model fits differ in field ''' equalfields{i} '''.']); end    
end

mfit = m1;  % Use M1 as base
hessianflag = 0;

% Merge optimization data
optimization = [];
optimization.chains = sort(unique([m1.optimization.chains,m2.optimization.chains]));
chainorigin = zeros(1,length(optimization.chains));
for iChain = 1:length(optimization.chains)
    c = optimization.chains(iChain);
    index1 = find(m1.optimization.chains == c,1);
    index2 = find(m2.optimization.chains == c,1);
    
    % Check for overlapping chains (weird but possible), take the best one
    if ~isempty(index1) && ~isempty(index2)
        best1 = m1.optimization.output{index1}.fvalmin;
        best2 = m2.optimization.output{index2}.fvalmin;
        if best1 <= best2; index2 = []; else index1 = []; end
    end        
    if ~isempty(index1) && isempty(index2)
        optimization.output{iChain} = m1.optimization.output{index1};
        chainorigin(iChain) = 1;
    else
        optimization.output{iChain} = m2.optimization.output{index2};
        chainorigin(iChain) = 2;
    end    
end
mfit.optimization = optimization;

% Find the best optimization chain
fval = zeros(1,length(optimization.chains));
index = zeros(1,length(optimization.chains));
for iChain = 1:length(optimization.chains)
    fval(iChain) = optimization.output{iChain}.fvalmin;
end
[~,bestChain] = min(fval);

% Update MAP and related model metrics
mfit.maptheta = optimization.output{bestChain}.xmin;
metricFields = {'H','Herr','aic','aicc','bic','maploglike','maploglikeSD',...
    'mapmarginallike','marginallike'};
if chainorigin(bestChain) == 1; metrics = m1.metrics;
else metrics = m2.metrics; end
for iField = 1:length(metricFields)
    mfit.metrics.(metricFields{iField}) = metrics.(metricFields{iField});
end

if isempty(mfit.sampling) && ~isempty(m2.sampling)
    mfit.sampling = m2.sampling;
    mfit = ModelWork_loadFields(project,mfit,data);    
    mfit = ModelWork_modelStats(project,mfit,hessianflag);    

elseif ~isempty(mfit.sampling) && ~isempty(m2.sampling) && ...
        (size(mfit.sampling.samples,1) > 0 || size(m2.sampling.samples,1) > 0)
    sampling = mfit.sampling;
    
    addfields = {'burnin','funccount','nchains','nsamplestot'};    
    for f = 1:numel(addfields)
        sampling.(addfields{f}) = sampling.(addfields{f}) + m2.sampling.(addfields{f});
    end    
    
    % LOGLIKES left for retrocompatibility
    concatfields = {'logliks','loglikes','logpriors','samples'};
    for f = 1:numel(concatfields)
        if isfield(m2.sampling, concatfields{f})
            sampling.(concatfields{f}) = [sampling.(concatfields{f}); m2.sampling.(concatfields{f})];
        end
    end
    
    % Summary sampling statistics
    if isfield(mfit, 'sumstats') && ~isempty(mfit.sumstats)
        n1 = sampling.sumstats.n;
        n2 = m2.sampling.sumstats.n;
        ntot = n1 + n2;        
        sumstats.smplmean1 = (mfit.sumstats.smplmean1*n1 + m2.sumstats.smplmean1*n2)/ntot;
        sumstats.loglikesmean1 = (mfit.sumstats.loglikesmean1*n1 + m2.sumstats.loglikesmean1*n2)/ntot;
        sumstats.loglikesmean1exp = (mfit.sumstats.loglikesmean1exp*n1 + m2.sumstats.loglikesmean1exp*n2)/ntot;
        sumstats.loglikessqsum1 = mfit.sumstats.loglikessqsum1 + m2.sumstats.loglikessqsum1;
        sumstats.n = mfit.sumstats.n + m2.sumstats.n;
        mfit.sampling.sumstats = sumstats;
    else
        mfit.sampling.sumstats = [];
    end
    
    mfit.sampling = sampling;
    
    % Update model statistics
    mfit = ModelWork_loadFields(project,mfit,data);    
    mfit = ModelWork_modelStats(project,mfit,hessianflag);
end

mfit.uptodate = 0;  % Might need a call to MODELWORK_MODELSTATS

end

