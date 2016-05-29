function [nLogL,extras] = ModelWork_like(project,X,mp,infostruct,theta,logpriorflag,trialloglikesflag,randomizeflag)
% MODELWORK_LIKE Negative log likelihood of the parameters.

if nargin < 3 % Input compatible with other functions
    if isnumeric(project)                % CMAES: y = fun(x,data)
        theta = project;
        mp = X.mp;
        infostruct = X.infostruct;
        if isfield(X,'project'); project = X.project; else project = X.prefix; end        
        X = X.X;
    else                                % MCS: y = fun(data,x)
        theta = X;
        X = project.X;
        mp = project.mp;
        infostruct = project.infostruct;
        if isfield(project,'project'); project = project.project; else project = project.prefix; end        
    end
end
if nargin < 6 || isempty(logpriorflag); logpriorflag = 0; end
if nargin < 7 || isempty(trialloglikesflag); trialloglikesflag = 0; end
if nargin < 8 || isempty(randomizeflag); randomizeflag = 0; end

INFPENALTY = -log(1e-6);

theta = theta(:)';

% Update the number of function calls
mp.funcalls = mp.funcalls + 1;

% Test boundary consistency for parameters
if any(theta < mp.bounds.LB) || any(theta > mp.bounds.UB)
    if isfield(mp, 'nData'); nData = mp.nData;
    elseif isfield(mp.extras, 'trialloglikes'); nData = length(mp.extras.trialloglikes);
    end
    if trialloglikesflag; nLogL = INFPENALTY*ones(1, nData); else nLogL = INFPENALTY*nData; end
    return;
end

setupModelFun = str2func([project '_setupModel']);
modelfitFun = str2func([project '_modelFitBits']);
logpriorFun = str2func([project '_logPrior']);

% Set global model parameters and check their consistency
[mp, outflag] = setupModelFun(mp,theta);
if outflag
    if isfield(mp, 'nData'); nData = mp.nData;
    elseif isfield(mp.extras, 'loglikes'); nData = length(loglikes);
    end
    if trialloglikesflag; nLogL = INFPENALTY*ones(1, nData); else nLogL = INFPENALTY*nData; end
    return; 
end

datalikeFun = modelfitFun('define', [], [], []);

if nargout > 1; flags = [1,randomizeflag]; else flags = [0,randomizeflag]; end
[loglike,extras] = datalikeFun(X,mp,infostruct,flags);            
mp.extras = extras;

if isnan(loglike)
    loglike
end

if trialloglikesflag && ...
    isfield(extras,'trialloglikes') && ~isempty(extras.trialloglikes)
        nLogL = -extras.trialloglikes;
else
    if logpriorflag; nLogL = -loglike - logpriorFun(mp,infostruct);
    else nLogL = -loglike; end
end

end