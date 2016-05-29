function mfit = ModelBag_mergefits(prefix,mfit,m)
% MODELBAG_MERGEMFITS merge two model fit structures.
%
%   MFIT = MODELBAG_MERGEFITS(MFIT,M) merge the base model fit MFIT with
%   model fit M for program PREFIX. Returns the merged MFIT.

% The following elements need to match
equalfields = {'nid', 'model', 'cnd', 'nsamplesperchain'};
for i = 1:length(equalfields)
    f1 = mfit.(equalfields{i});
    f2 = m.(equalfields{i});
    if ~all(f1(:) == f2(:)); error(['Cannot merge: Model fits differ in field ''' equalfields{i} '''.']); end    
end

nstoredsamples = size(mfit.smpl, 1);

if mfit.nsamplesperchain < m.nsamplesperchain && 0
    maxsamples = mfit.nsamplesperchain*m.nchains;   
    step = size(m.smpl, 1)/maxsamples;
    smplmask = round(step:step:size(m.smpl, 1));

    % Thin samples and associated log likelihood
    m.smpl = m.smpl(smplmask, :);
    m.loglikes = m.loglikes(smplmask, :);                        
end

mfit.smpl = [mfit.smpl; m.smpl];
mfit.loglikes = [mfit.loglikes; m.loglikes];
mfit.lastsmplperchain = [mfit.lastsmplperchain; m.lastsmplperchain];

n1 = mfit.nsamplesperchain*mfit.nchains;
n2 = m.nsamplesperchain*m.nchains;
ntot = n1 + n2;

% Summary sampling statistics
if isfield(mfit, 'sumstats') && ~isempty(mfit.sumstats)
    sumstats.smplmean1 = (mfit.sumstats.smplmean1*n1 + m.sumstats.smplmean1*n2)/ntot;
    sumstats.loglikesmean1 = (mfit.sumstats.loglikesmean1*n1 + m.sumstats.loglikesmean1*n2)/ntot;
    sumstats.loglikesmean1exp = (mfit.sumstats.loglikesmean1exp*n1 + m.sumstats.loglikesmean1exp*n2)/ntot;
    sumstats.loglikessqsum1 = mfit.sumstats.loglikessqsum1 + m.sumstats.loglikessqsum1;
    sumstats.n = mfit.sumstats.n + m.sumstats.n;
    mfit.sumstats = sumstats;
else
    mfit.sumstats = [];
end

mfit.nchains = mfit.nchains + m.nchains;

% Update model statistics
hessianflag = ~isempty(mfit.H);
mfit = ModelWork_modelStats(prefix,mfit,hessianflag);

end

