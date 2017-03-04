function [nLogL,extras] = ModelWork_like(varargin)
%MODELWORK_LIKE Negative log likelihood of the parameters.
%   NLOGL = MODELWORK_LIKE(PROJECT,X,MP,INFOSTRUCT,THETA) returns the
%   negative log likelihood of dataset X with model parameter structure MP
%   and extra structure INFOSTRUCT, evaluate at parameter THETA, for project
%   PROJECT.
%
%   NLOGL = MODELWORK_LIKE(PROJECT,MFIT,THETA) returns the negative log 
%   likelihood for model in MFIT evaluated at parameter THETA. If THETA is
%   empty, the log likelihood is evaluated at the maximum-likelihood solution.
%
%   NLOGL = MODELWORK_LIKE(MFIT,THETA) or
%   NLOGL = MODELWORK_LIKE(THETA,MFIT) is a compact way to call the
%   function, compatible with external programs such as CMAES and MCS.
%   MFIT needs to have an additional field MFIT.project which contains the
%   project name.
%
%   NLOGL = MODELWORK_LIKE(...,LOGPRIORFLAG,TRIALLOGLIKESFLAG,RANDOMIZEFLAG) 
%   sets the flags (all zero by default). LOGPRIORFLAG adds the log prior.
%   TRIALLOGLIKESFLAG returns log likelihood per trial. RANDOMIZEFLAG adds
%   randomization to log likelihood computation.
%
%   [NLOGL,EXTRAS] = MODELWORK_LIKE(...)  returns additional information
%   in structure EXTRAS.
%
%   See also MODELWORK_SETUPMODEL.

if nargin == 0 && nargout == 0; help ModelWork_like; return; end

% Initialize variables
project = [];
mfit = [];
X = [];
mp = [];
infostruct = [];
theta = [];
logpriorflag = [];
trialloglikesflag = [];
randomizeflag = [];

% Standard input format: PROJECT,MFIT,THETA or PROJECT,X,MP,INFOSTRUCT,THETA
if ischar(varargin{1})  % project
    project = varargin{1};
    % First input is a MFIT struct
    if all(isfield(varargin{2},{'X','mp','infostruct','maptheta'}))
        mfit = varargin{2};
        if nargin > 2
            theta = varargin{3};
            if ~isempty(theta) && (size(theta,2) ~= size(mfit.maptheta,2))
                error('THETA does not have the right dimension.');
            end
        end
        idx = 4;
    else
        X = varargin{2};
        mp = varargin{3};
        infostruct = varargin{4};
        theta = varargin{5};
        idx = 6;
    end
    % Read additional flags
    if nargin >= idx; logpriorflag = varargin{idx}; end
    if nargin >= idx+1; trialloglikesflag = varargin{idx+1}; end
    if nargin >= idx+2; randomizeflag = varargin{idx+2}; end
    
% Alternative input format: MFIT,THETA or THETA,MFIT
elseif nargin < 3 % Input compatible with other functions
    if isnumeric(varargin{1})               % CMAES: y = fun(x,data)
        theta = varargin{1};
        mfit = varargin{2};
    else                                    % MCS: y = fun(data,x)
        mfit = varargin{1};
        theta = varargin{2};
    end
end

if ~isempty(mfit)
    X = mfit.X;
    mp = mfit.mp;
    infostruct = mfit.infostruct;
    if isempty(project)
        if isfield(mfit,'project'); project = mfit.project;
        elseif isfield(mfit,'prefix'); project = mfit.prefix;
        else
            error('MFIT structure should have a ''project'' field if calling MODELWORK_LIKE with only two arguments.');
        end
    end
    if isempty(theta); theta = mfit.maptheta; end
end

if isempty(logpriorflag); logpriorflag = 0; end
if isempty(trialloglikesflag); trialloglikesflag = 0; end
if isempty(randomizeflag); randomizeflag = 0; end

INFPENALTY = -log(1e-6);

if size(theta,1) > 1 && size(theta,2) > 1 && nargout == 1
    % Vectorized call (EXTRAS not accepted)
    N = size(theta,1);    
    tmp = ModelWork_like(project,X,mp,infostruct,theta(1,:),logpriorflag,trialloglikesflag,randomizeflag);
    nLogL = NaN(N,numel(tmp));
    nLogL(1,:) = tmp;
    for i = 2:N
        nLogL(i,:) = ModelWork_like(project,X,mp,infostruct,theta(i,:),logpriorflag,trialloglikesflag,randomizeflag);
    end
    return;
end

theta = theta(:)';
extras = [];

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