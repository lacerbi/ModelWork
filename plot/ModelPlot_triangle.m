function ModelPlot_triangle(mfit)
%MODELPLOT_TRIANGLE Triangle plot of model posterior.

X = mfit.sampling.samples;

if ~isempty(mfit.mp)
    names = mfit.mp.params;
    for i = 1:numel(names); names{i}(names{i} == '_') = '-'; end
    
    lb = mfit.mp.bounds.RLB;
    ub = mfit.mp.bounds.RUB;    
    sample_lb = min(X);
    sample_ub = max(X);
    
    bounds = [min(sample_lb,lb); max(sample_ub,ub)];
    
else
    names = [];
    bounds = [];
end

truths = mfit.maptheta;

cornerplot(X,names,truths,bounds)
set(gcf,'Color','w');

end
