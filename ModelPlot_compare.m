function ModelPlot_compare(modelsummary,metric,baseline)
%MODELPLOT_COMPARE Plot results of model comparison.

if nargin < 3 || isempty(baseline); baseline = 0; end

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

bar(y);

set(gcf,'Color','w');
set(gca,'TickDir','out');
box off;
xlabel('Subject');
ylabel(['\Delta' metricstring]);
legend(modelnames{:}); 
if baseline == 0
    title([metricstring]);
else
    title([metricstring ' (wrt ' basename ' model)']);
end






end