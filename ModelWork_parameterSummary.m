function paramstat = ModelWork_parameterSummary(mbag,modelsummary,paramname,bms,fullthetanumber,modellist)
% MODELWORK_PARAMETERMEAN Compute mean and SD of parameter across subjects.

if ~exist('bms', 'var'); bms = []; end
if ~exist('modellist', 'var') || isempty(modellist)
    if isempty(bms)
        modellist = modelsummary.models;
    else
        modellist = bms.models;
    end
end
if ~exist('fullthetanumber', 'var') || isempty(fullthetanumber); fullthetanumber = 1; end

subjs = modelsummary.nid;
parvalue = NaN(size(modellist,1), length(subjs));

for iModel = 1:size(modellist,1)
    for iSubj = 1:length(subjs)
        model = modellist(iModel,:);
        mfit = ModelBag_get(mbag,subjs(iSubj),model,{modelsummary.cnd});
        if length(mfit.mp.fulltheta) >= fullthetanumber && ...
                isfield(mfit.mp.fulltheta{fullthetanumber}, paramname)            
            parvalue(iModel, iSubj) = mfit.mp.fulltheta{fullthetanumber}.(paramname);
            
            iPos = find(all(bsxfun(@eq, model, bms.models), 2),1);
            if isempty(bms)
                postp(iModel, iSubj) = 1;
            else
                postp(iModel, iSubj) = bms.g(subjs(iSubj), iPos);
            end
        end
    end
end

paramstat.wp = nansum(postp, 1); % Weighted contribution of each subject
paramstat.means_ss = nansum(parvalue.*postp, 1)./paramstat.wp;
paramstat.mean = nansum(paramstat.means_ss.*paramstat.wp)/nansum(paramstat.wp);

% Compute standard error of the weighted mean by bootstrap
nsmpl = 10^3;
temp = postp; temp(isnan(temp)) = 0;
if isempty(bms); temp = bsxfun(@rdivide, post, sum(temp,1)); end

% Should be able to do this analytically?
for iSubj = 1:length(subjs)
    absentp = max(0, 1 - nansum(temp(:,iSubj)));
    multinp = [temp(:,iSubj)' absentp]/sum([temp(:,iSubj)' absentp]);
    [~,iPos] = nanmax(mnrnd(1, ones(nsmpl,1)*multinp), [], 2);
    temppar = [parvalue(:, iSubj); NaN];
    bootparvalue(iSubj,:) = temppar(iPos)';    
end
paramstat.serr = nanmean(stderr(bootparvalue,[],1));
paramstat.sd = sqrt(nanmean(nanvar(bootparvalue,[],1)));
% paramstat.serr = sqrt(nanmean(stderr(bootparvalue).^2));

end