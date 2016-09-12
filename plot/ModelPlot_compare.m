function ModelPlot_compare(modelsummary,metric,baseline,action,idx,modelidx)
%MODELPLOT_COMPARE Plot results of model comparison.

if nargin < 3 || isempty(baseline); baseline = 0; end
if nargin < 4 || isempty(action); action = 'mean'; end
if nargin < 5; idx = []; end
if nargin < 6; modelidx = []; end

bootflag = 1;   % Use bootstrap for error bars

% It works with a MBAG too
if isstruct(modelsummary) && isfield(modelsummary,'bag')
    modelsummary = ModelWork_summary(modelsummary);
end

colormap(gray);

modelnumber = 1:size(modelsummary.models,1);
modelnames = modelsummary.modelnames;
metricstring = upper(metric);

if ischar(baseline)
    baseline = find(strcmp(baseline,modelsummary.modelnames));
    if isempty(baseline) || ~isscalar(baseline)
        error('BASELINE model name should be present and unique in the model table.');
    end
end

if baseline == 0
    y = modelsummary.(metric);    
else
    y = bsxfun(@minus,modelsummary.(metric),modelsummary.(metric)(:,baseline));
    y(:,baseline) = [];
    basename = modelnames{baseline};
    modelnames(baseline) = [];
end

% Take only specified datasets
if ~isempty(idx); y = y(idx,:); end
if ~isempty(modelidx)
    if iscell(modelidx)
        newmodelidx = [];
        for i = 1:numel(modelidx)
            if strcmp(modelidx{i},basename); continue; end  % Skip baseline
            tmp = find(strcmp(modelidx{i},modelnames));
            if isempty(tmp) || ~isscalar(tmp)
                error(['Cannot find a unique match for model ' modelidx{i} ' in model table.']);
            end
            newmodelidx = [newmodelidx,tmp];
        end
        modelidx = newmodelidx;
    end
    y = y(:,modelidx);
    modelnames = modelnames(modelidx)
end

[Nsubjs,Nmodels] = size(y);

if bootflag
    nboot = 2e4;
    means = zeros(nboot,Nmodels);
    for i=1:nboot
       bdata = y(ceil(rand(Nsubjs,1)*Nsubjs),:);
       means(i,:) = mean(bdata);
    end
    alpha = 5/Nmodels;
    CI = prctile(means,[alpha/2,100-alpha/2]);
    p = CI(1,:) < 0 & CI(2,:) > 0;
    yerr = std(means);    
else
    yerr = stderr(y,1);
end

switch lower(action(1))
    case 's'    % SUM
        y = sum(y,1);
        groupflag = 1;
    case 'm'    % MEAN
        y = mean(y,1);
        groupflag = 1;
end


if groupflag
    col = 0.7*[1 1 1];
    
    bar(y,'FaceColor',col,'EdgeColor','none'); hold on;
    y_notp = y;
    y_notp(~p) = 0;
    bar(y_notp,'FaceColor','c','EdgeColor','none');    
    for i = 1:Nmodels
        plot(i*[1 1], y(i)+[-1,1]*yerr(i),'k','LineWidth',1);
    end
    % xlabel('Models');
    ylabel(['\Delta' metricstring])
    set(gca,'XTick',1:numel(modelnames),'XTickLabel',modelnames);
    
    text(0.05,1,['n = ' num2str(Nsubjs)],'Units','Normalized');
    xticklabel_rotate([],90);
else

    bar(y);
    xlabel('Subject');
    ylabel(['\Delta' metricstring]);
    legend(modelnames{:}); 
end

if baseline == 0
    title([metricstring]);
else
    title([metricstring ' (wrt ' basename ' model)']);
end

set(gcf,'Color','w');
set(gca,'TickDir','out');
box off;






end