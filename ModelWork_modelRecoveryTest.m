% MODELWORK_MODELRECOVERYTEST Analyze results of model recovery test.
%
% MODELWORK_BATCHRUN(PREFIX,ID,TYPE,IDTOT) runs batch with identified ID out 
% of IDTOT. TYPE specifies the TYPE of the run (see below). ID can be a
% vector of process identifiers, in which case the runs are sequential.

% By Luigi Acerbi <luigi.acerbi@gmail.com>
% Last update: Jan/29/2015

function [recomatrix,err] = ModelWork_modelRecoveryTest(fakedata, mbag)

% Check fake data and make list of true generative models
truemodels = fakedata{1}.truemodel;
mdata{1}{1} = fakedata{1};
subjtruemodel(1) = 1;
for i = 2:length(fakedata)
    index = find(all(bsxfun(@eq, fakedata{i}.truemodel, truemodels), 2));
    if isempty(index)
        truemodels = [truemodels; fakedata{i}.truemodel]; 
        index = size(truemodels, 1); 
        mdata{index}{1} = fakedata{i};
        
    else
        mdata{index}{end+1} = fakedata{i};
    end
    subjtruemodel(fakedata{i}.id) = index;
end

% Check model recovery data and make list of used models
models = mbag.bag{1}.model;
mmfit{1}{mbag.bag{1}.nid} = mbag.bag{1};
for i = 2:length(mbag.bag)
    index = find(all(bsxfun(@eq, mbag.bag{i}.model, models), 2));
    if isempty(index)
        models = [models; mbag.bag{i}.model]; 
        index = size(models, 1); 
        mmfit{index}{mbag.bag{i}.nid} = mbag.bag{i};
    else
        mmfit{index}{mbag.bag{i}.nid} = mbag.bag{i};
    end    
end

% Get all subjects id and conditions
for i = 1:length(mbag.bag)
    fakenid(i) = mbag.bag{i}.nid;
    fakecnd{fakenid(i)} = mbag.bag{i}.cnd; 
end
fakenid = sort(unique(fakenid));

% Should maybe order models as truemodels array, but shall do later
% ...

% Build model fit metric for each individual fake subject
for i = 1:length(fakenid)
   for j = 1:size(models, 1)
       nid = fakenid(i);
       mm = ModelBag_get(mbag, nid, models(j, :), {fakecnd{nid}});
       score(j, i) = mm.aic;
   end
end

% Winning model per subject (minimum score)
[~,bestmodel] = min(score);

% Create model recovery matrix
recomatrix = zeros(size(truemodels, 1), size(models, 1));

for i = 1:length(bestmodel)
    nid = fakenid(i);
    recomatrix(subjtruemodel(nid), bestmodel(i)) = recomatrix(subjtruemodel(nid), bestmodel(i)) + 1;    
end
recomatrix = bsxfun(@rdivide, recomatrix, sum(recomatrix, 2));

% Parameter recovery error (for matching models only)
for i = 1:size(truemodels, 1); err{i} = []; end
for i = 1:length(fakedata)
    D = fakedata{i};
    model = find(all(bsxfun(@eq, models, D.truemodel), 2));
    mfit = mmfit{model}{D.id};
    err{model}(end+1, :) = mfit.maptheta - D.truetheta;
end

end
