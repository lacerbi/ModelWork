function ModelPlot_compare(modelsummary,metric,baseline,action,idx)
%MODELPLOT_COMPARE Plot results of model comparison.

if nargin < 3 || isempty(baseline); baseline = 0; end
if nargin < 4 || isempty(action); action = 'mean'; end
if nargin < 5; idx = []; end

% It works with a MBAG too
if isstruct(modelsummary) && isfield(modelsummary,'bag')
    modelsummary = ModelWork_summary(modelsummary);
end

colormap(gray);

modelnumber = 1:size(modelsummary.models,1);
metricstring = upper(metric);

if baseline == 0
    y = modelsummary.(metric);
    modelnames = modelsummary.modelnames;
else
    y = bsxfun(@minus,modelsummary.(metric),modelsummary.(metric)(:,baseline));
    y(:,baseline) = [];
    basename = modelsummary.modelnames{baseline};
    modelnames = modelsummary.modelnames(setdiff(modelnumber,baseline));
end

% Take only specified datasets
if ~isempty(idx); y = y(idx,:); end

[nsubjs,nmodels] = size(y);

yerr = stderr(y,1);

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
    for i = 1:nmodels
        plot(i*[1 1], y(i)+[-1,1]*yerr(i),'k','LineWidth',1);
    end
    xlabel('Models');
    ylabel(['\Delta' metricstring])
    set(gca,'XTickLabel',modelnames);
    
    text(0.05,1,['n = ' num2str(nsubjs)],'Units','Normalized');
    
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