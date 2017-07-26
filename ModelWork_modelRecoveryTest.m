% MODELWORK_MODELRECOVERYTEST Analyze results of model recovery test.
%
% MODELWORK_BATCHRUN(PREFIX,ID,TYPE,IDTOT) runs batch with identified ID out 
% of IDTOT. TYPE specifies the TYPE of the run (see below). ID can be a
% vector of process identifiers, in which case the runs are sequential.

% By Luigi Acerbi <luigi.acerbi@gmail.com>
% Last update: Mar/06/2016

function [recomatrix,err] = ModelWork_modelRecoveryTest(fakedata,mbag,metric,modelnames,fontsize,axesfontsize)

if nargin < 3 || isempty(metric); metric = 'aicc'; end
if nargin < 4 || isempty(modelnames); modelnames = []; end
if nargin < 5 || isempty(fontsize); fontsize = 18; end
if nargin < 6 || isempty(axesfontsize); axesfontsize = 14; end

% Comparison metric sign
if any(strncmpi(metric,{'aic','bic','dic','wai'},3))
    signum = -1;
else
    signum = 1;
end

project = mbag.project;
getModelName = str2func([project '_getModelName']);

% Number of fake datasets per subject
for i = 1:numel(fakedata)
    dataid(i) = fakedata{i}.dataid;
    fakeid(i) = fakedata{i}.fakeid;
end

Nfaketot = numel(fakedata);
Nsubjs = max(dataid);
Nfakeiter = max(fakeid);

% Check fake data and make list of true generative models
truemodels = fakedata{1}.truemodel;
mdata{1}{1} = fakedata{1};
subjtruemodel(1) = 1;
for i = 2:numel(fakedata)
    % subjid = (fakedata{i}.dataid-1)*Nfakeiter + fakedata{i}.fakeid;
    
    index = find(all(bsxfun(@eq, fakedata{i}.truemodel, truemodels), 2));
    if isempty(index)
        truemodels = [truemodels; fakedata{i}.truemodel]; 
        index = size(truemodels, 1); 
        mdata{index}{1} = fakedata{i};        
    else
        mdata{index}{end+1} = fakedata{i};
    end
    subjtruemodel(i) = index;
end

Ntruemodels = size(truemodels,1);
for i = 1:Ntruemodels
    truemodelnames{i} = getModelName(truemodels(i,:));
end

truemodelnames
if isempty(modelnames); modelnames = truemodelnames; end

% Check model recovery data and make list of fitted models
models = mbag.bag{1}.model;
subjid = mbag.bag{1}.dataid;
mmfit{1}{subjid} = mbag.bag{1};
for i = 2:numel(mbag.bag)
    subjid = mbag.bag{i}.dataid;    
    index = find(all(bsxfun(@eq, mbag.bag{i}.model, models), 2));
    if isempty(index)
        models = [models; mbag.bag{i}.model]; 
        index = size(models, 1); 
        mmfit{index}{subjid} = mbag.bag{i};
    else
        mmfit{index}{subjid} = mbag.bag{i};
    end    
end

Nfitmodels = size(models,1);

% Get all subjects id and conditions
%for i = 1:numel(mbag.bag)
%    fakenid(i) = mbag.bag{i}.nid;
%    fakecnd{fakenid(i)} = mbag.bag{i}.cnd; 
%end
%fakenid = sort(unique(fakenid));

% Order models as truemodels array
ord = [];
for i = 1:Ntruemodels
    idx = find(all(bsxfun(@eq, truemodels(i,:), models), 2),1);
    if ~isempty(idx); ord = [ord, idx]; end    
end
ord = [ord, setdiff(1:Nfitmodels,ord)];

% Build model fit metric for each individual fake subject
score = NaN(Nfitmodels,Nfaketot);
for j = 1:Nfitmodels       
    for i = 1:min(Nfaketot,numel(mmfit{ord(j)}))
       mm = mmfit{ord(j)}{i};
       if ~isempty(mm)
            score(j, i) = mm.metrics.(metric);
       end
   end
end

% Winning model per subject
[maxscore,bestmodel] = max(signum*score,[],1);

% Create model recovery matrix
recomatrix = zeros(Ntruemodels, Nfitmodels);
n = zeros(Ntruemodels, Nfitmodels);

for i = 1:numel(bestmodel)
    if any(isfinite(score(:,i)))
        recomatrix(subjtruemodel(i), bestmodel(i)) = recomatrix(subjtruemodel(i), bestmodel(i)) + 1;
%        recomatrix(subjtruemodel(i), bestmodel(i)) = recomatrix(subjtruemodel(i), bestmodel(i)) + maxscore(i);
        n(subjtruemodel(i), bestmodel(i)) = n(subjtruemodel(i), bestmodel(i)) + 1;
    end
end
recomatrix = bsxfun(@rdivide, recomatrix, sum(recomatrix, 2));

% Parameter recovery error (for matching models only)
for i = 1:Ntruemodels; err{i} = []; end
for i = 1:Nfaketot
    D = fakedata{i};
    if isempty(D); continue; end
    model = find(all(bsxfun(@eq, models, D.truemodel), 2));
    mfit = mmfit{model}{D.dataid};
    if isempty(mfit); continue; end
    err{model}(end+1, :) = mfit.maptheta - D.trueparams;
end

% Plot model recovery matrix

for i = 1:Nfitmodels
    for j = 1:Ntruemodels
        p = recomatrix(j,i);        
        patch([i i i+1 i+1]-0.5, [-j -j-1 -j-1 -j]+0.5, p, 'EdgeColor', 'none');
        if p < 0.01; continue; end        
        if p > 0.5; col = [0 0 0]; else col = 0.6*[1 1 1]; end
        text(i,-j,num2str(p,'%.2f'),'FontSize',axesfontsize,'HorizontalAlignment','center','Color',col);
    end
end
axis([0.5,Ntruemodels+0.5,-Nfitmodels-0.5,-0.5]);
if Nfitmodels == Ntruemodels; axis square; end

set(gca,'XTick',1:Ntruemodels,'YTick',-Nfitmodels:-1,'FontSize',axesfontsize);
if ~isempty(modelnames)
    set(gca,'XTickLabel',modelnames);
    set(gca,'YTickLabel',modelnames(end:-1:1));
end
xticklabel_rotate([],45);

xlabel('Fitted models','FontSize',fontsize);
ylabel('Generating models','FontSize',fontsize);
title('Fraction recovered','FontSize',fontsize);

grid off;
colormap(gray);
set(gca,'TickDir','out');
set(gcf,'Color','w');


end
