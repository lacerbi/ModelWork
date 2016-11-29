function mbag = ModelBag_convert(mbag,project)
%MODELBAG_CONVERT Covert obsolete model bag to new format.
%   MBAG = MODELBAG_CONVERT(MBAG,PROJECT) convert the model bag MBAG (in 
%   pre-2015 format) to the current model bag format, for project PROJECT.

mbag.project = project;
for i = 1:numel(mbag.bag)
    m = mbag.bag{i};
    m.dataid = m.nid;
    m = rmfield(m,'nid');
    
    metrics = {'aic','aicc','bic','H','Herr','maploglike','marginallike','mapmarginallike','maploglike_sd'};
    for j = 1:numel(metrics)
        m.metrics.(metrics{j}) = m.(metrics{j});
    end
    m = rmfield(m, metrics);    
    
    sampling = {'nsamplesperchain','nchains','nstoredsamples','smpl','loglikes','sumstats'};
    for j = 1:numel(sampling)
        m.sampling.(sampling{j}) = m.(sampling{j});
    end
    m = rmfield(m, sampling);    
    
    mbag.bag{i} = m;
end



end