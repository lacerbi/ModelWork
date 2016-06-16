function [options,models,dataids,cnd] = Dummy_batchrunInit(data,type,options)
% DUMMY_BATCHRUNINIT initialize variables for batch run.

if nargin < 3; options = []; end

debug = 0;

% Get additional type parameters
type = type(1);

% Continue previous sampling if exists
options = setoptions(options,'loadstartx',1,1);

% By default do not compute Hessian unless optimizing (see below)
options = setoptions(options,'hessianflag',0,1);

options = setoptions(options,'optimizationmethod','fmincon',1);
options = setoptions(options,'FitGMM',0,1);

% Get host name
[~,hostname] = system('hostname');
options = setoptions(options,'hostname',strtrim(hostname),1);

% Number of datasets
if isfield(data,'data'); nDatasets = length(data.data);
else nDatasets = length(data); end

% Subjects mask
DATAIDS = [(1:nDatasets)',zeros(nDatasets,1)];

% Default number of samples for MCMC
NSAMPLES = 1000;

% Default optimization steps before starting MCMC
MAXFUNEVALS = 500;

% Optimization steps when optimizing only
NITER_OPTIMIZATION = 500;

% Number of restarts for optimization
nOptimizationRestarts = 3;

if debug
    nOptimizationRestarts = 10;
    NSAMPLES = [10,10];
    MAXFUNEVALS = 10;
end

options = setoptions(options,'nstarts',nOptimizationRestarts,1);
options = setoptions(options,'optfevals',MAXFUNEVALS,1);
options = setoptions(options,'nsamples',NSAMPLES,1);

dataids = DATAIDS;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT MODELS

models = [];
standardmodels = [ 1; ...   % No lapse
    2];                     % Lapse

switch type
    case 0; % DEBUG    

        if nc == 0 % Debug unimodal data     
            options = setoptions(options,'nc',2,1);
            options = setoptions(options,'samplingtemperature',2,1);
            options = setoptions(options,'nsamples',10,1);
            options = setoptions(options,'nstoredsamples',15,1);
            options = setoptions(options,'optfevals',20,1);
            dataids = dataids(1,:);
            models = standardmodels;
            groupcnd = 1;
        end
        
%--------------------------------------------------------------------------        
% BISENSORY ESTIMATION DATA FITS
    
   case {1} % Bisensory standard models
       
       models = standardmodels;
       groupcnd = 1;
       options.jobname = 'base';              
        
end

% Set speed test values
%if any(groupcnd >= 5)
%    options = setoptions(options,'speedtest',10,1);  % Bimodal data
%else
%    options = setoptions(options,'speedtest',0.15,1);  % Unimodal data
%end

% Be verbose
options = setoptions(options,'display','all',1);

% Set conditions
for i = 1:nDatasets; cnd{i} = {groupcnd}; end

% Optimization run
if type > 0 && options.nsamples == 0
    options = setoptions(options,'optfevals',NITER_OPTIMIZATION,1);
    options = setoptions(options,'nsamples',0,1);
    options = setoptions(options,'samplingtemperature',1,1);
    options = setoptions(options,'hessianflag',0,1); % Too expensive - do not compute Hessian
end
        
end    