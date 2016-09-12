function [mfit,exitflag] = ModelWork_modelFit(project,data,mfit,method,options,varargin)
% MODELWORK_MODELFIT Fit single model to data.
%
%   MFIT = MODELWORK_MODELFIT(PROJECT,DATA,MFIT,METHOD,OPTIONS) fits the 
%   dataset DATA for project PROJECT starting from model structure MFIT
%   and with method METHOD, according to the struct OPTIONS. METHOD can
%   be either 'opt' for optimization or 'sam' for sampling.
% 
%   MFIT = MODELWORK_MODELFIT(...,NAME,VALUE) sets the 
%   field NAME of the options struct to a given VALUE (ovverriding settings
%   in OPTIONS).
%
%   [MFIT,EXITFLAG] = MODELWORK_MODELFIT(...) returns an EXITFLAG that 
%   describes the exit condition of MODELWORK_MODELFIT. A value of EXITFLAG 
%   of 0 means a correct execution; a value of 1 indicates an error.

% By Luigi Acerbi <luigi.acerbi@gmail.com>
% Last update: 22/08/2015

if nargin < 1
    help ModelWork_modelFit;
    return;
end

% Load options (MODELFIT can be called independently of BATCHEVAL)
optlist = ModelWork_defaults(project);
options = parseopts(varargin,optlist,options);

% Store model fit method
options = setoptions(options,'fitmethod',method,0);

% Data structure with individual datasets
if isfield(data,'data'); data = data.data; end
dataone = data{options.dataid(1)};

exitflag = 0; % Everything is fine

if options.clearfunctions
    clear functions; % Clear persistent variables (does it by default)
end

switch lower(method(1:3))
    case 'opt'              % Optimizing
        samplingflag = 0;
    case {'sam','smp'}      % Sampling 
        samplingflag = 1;
end

% Verbosity and log file
if isempty(options.outfile); fout = 1; % Standard output
else fout = fopen(options.outfile, 'a'); end

% Random seed (use a fixed seed for debugging purposes)
defaultSeed = options.replica + prod(max(options.model,1)) + 10*options.type + 1000*samplingflag;
options = setoptions(options,'seed',defaultSeed,1);
try RandStream.setGlobalStream(RandStream.create('mt19937ar','seed',options.seed));
catch; RandStream.setDefaultStream(RandStream.create('mt19937ar','seed',options.seed)); end
writelog(fout,'seed',options);

% Define project-dependent functions
modelfitFun = str2func([project '_modelFitBits']);
setupModelFun = str2func([project '_setupModel']);
logpriorFun = str2func([project '_logPrior']);

% Define model- and data- dependent variables
[datalikeFun,infostruct,options] = modelfitFun('define',dataone,options.model,options);
infostruct.cnd = options.cnd;

% Check that the model string is valid
[~,exitflag] = setupModelFun([],[],options.model,infostruct);
if exitflag; return; end

% Initialize model-fit structure
if isempty(mfit)
    [mfit,loadedstartx] = initmodelfit(options,fout); 
    if ~isempty(loadedstartx)
        options.startx = loadedstartx;
        % Do not perform optimization if sampling from loaded point
        if samplingflag; options.optfevals = 0; end
    end
else
    % writelog(fout,'continue',options);
    
    % With MCMC skip burn-in and optimization if starting from a previous sample
    %if samplingflag
        % options.optfevals = 0;
        % Start from the end of the last chains
    % If performing optimization only, start from previous points
    %else
    %    options.startx = options.mfit.optimization.x;        
    %end    
    
end

% Pre-process dataset according to model
[mfit,infostruct] = modelfitFun('preprocessdata',dataone,mfit,options,infostruct);

% Initialize model parameters
[mp,exitflag] = setupModelFun([],[],options.model,infostruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

maptheta = mfit.maptheta; fvalmin = Inf;
if ~isempty(mfit.metrics) && ~isempty(mfit.metrics.maploglike)
    fvalmin = -mfit.metrics.maploglike;
end

% First iteration from midpoint (by default only for first replica)
if isempty(options.firstSetting); firstSetting = (options.replica == 1);
else firstSetting = options.firstSetting; end

% Speed test
if ~isempty(options.speedtest) && options.speedtest > 0
    nvars = mp.ntotparams;
    for i = 1:5
        startx = mp.bounds.RLB + rand(1,nvars).*(mp.bounds.RUB-mp.bounds.RLB);
        t = tic;
        ModelWork_like(project,mfit.X,mp,infostruct,startx,1,0,0); 
        testtime(i) = toc(t);
    end
    avgtime = mean(testtime(2:end));
    if avgtime < options.speedtest
        writelog(fout,'speedsuccess',options,avgtime);
    else
        writelog(fout,'speedfail',options,avgtime);
        error('ModelWork:SpeedFail', 'Failed speed test.');
    end
end

switch lower(method(1:3))
    
    case 'opt'
        % Multi-start optimization for MAP estimate (or MLE with flat priors)
        [maptheta,fvalmin,optexitflag,optoutput] = ...
            ModelWork_map(mfit,mp,infostruct,fout,datalikeFun,setupModelFun,logpriorFun,firstSetting,options);
        mp.funcalls = mp.funcalls + sum(optoutput.funccount);
        
        mfit.optimization.chains = options.replica;
        mfit.optimization.output{1} = optoutput;    

    case {'sam','smp'}
        % Starting points for sampling
        sample0 = []; loglikes = []; weights = [];
        if ~isempty(mfit.sampling) && ~isempty(mfit.sampling.samples)
            sample0 = mfit.sampling.samples;
        elseif ~isempty(mfit.optimization) && ~isempty(mfit.optimization.output)
            for k = 1:numel(mfit.optimization.output)
                sample0 = [sample0; mfit.optimization.output{k}.x];
                loglikes = [loglikes; -mfit.optimization.output{k}.fval(:)];
            end
            weights = exp(loglikes - max(loglikes)); % Weights assigned to starting points
        end
        
        % If passing an external fit, use samples as starting points but
        % otherwise sample from scratch
        if strcmpi(options.mbag,'used'); mfit.sampling = []; end
        
        [sampling,smplexitflag,smploutput] = ...
            ModelWork_sample(mfit,mp,infostruct,fout,datalikeFun,setupModelFun,logpriorFun,sample0,weights,firstSetting,options);    
        mfit.sampling = sampling;
        mp.funcalls = mp.funcalls + sum(smploutput.funccount);
end

mp.computation = 'normal'; % Switch to normal computational precision

% Compute MAP from optimization
mfit.maptheta = maptheta(:)';
mfit.metrics.maploglike = -fvalmin;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING DIAGNOSTICS AND PARAMETER ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute sampled chains and model summary statistics
mfit.mp = mp; mfit.infostruct = infostruct;
mfit.uptodate = 0;
if options.hessianflag; writelog(fout,'hessian',options); end
mfit = ModelWork_modelStats(project,mfit,options.hessianflag,[],[],options);

% Print sampling results
if samplingflag
    mfit.sampling
    writelog(fout,'samplingresults',options,mp.funcalls,mfit.metrics.maploglike,mfit.metrics,mfit.sampling.meantheta); 
end

% Close file
if fout > 1; fclose(fout); end

% Remove temporary fields from model struct
mfit = ModelWork_dropFields(mfit,options);

end % MODELWORK_MODELFIT

%--------------------------------------------------------------------------
function [mfit,loadedstartx] = initmodelfit(options,fout)
%INITMODEL Initialize model-fit structure.

mfit = [];
mFields = {'type','dataid','model','cnd','datafile','nData','X','mp',...
    'infostruct','optimization','sampling','metrics','maptheta'};
for iField = 1:length(mFields); mfit.(mFields{iField}) = []; end

mfit.type = options.type;
mfit.dataid = rmtrailzeros(options.dataid);
mfit.model = rmtrailzeros(options.model);
mfit.cnd = options.cnd;
mfit.datafile = options.datafile;
mfit.uptodate = 0;

% Load starting point for optimization or sampling
loadedstartx = [];
if options.loadstartx && exist('startx.mat', 'file')
    load('startx.mat', 'startxmbag');
    mfittemp = ModelBag_get(startxmbag,options.dataid,options.model,options.cnd);

    if ~isempty(mfittemp)
        writelog(fout,'startx',options);

        if ~isempty(mfittemp.sampling) && ~isempty(mfittemp.sampling.samples)
            % If sampling data are available, start from last sample
            loadedstartx = mfittemp.sampling.samples(end,:);
        elseif ~isempty(mfittemp.maptheta) && all(isfinite(mfittemp.maptheta))
            % If optimization data are available, start from optimum
            loadedstartx = mfittemp.maptheta;
        end
    end
end

end

%--------------------------------------------------------------------------
function writelog(fout,state,options,varargin)
%WRITELOG Write entry in log file.

if isempty(fout); return; end
if ~strcmpi(options.display, 'all'); return; end

modelstring = options.modelstring;
dataidstring = options.dataidstring;

if strncmpi(state,'speed',5)
    if options.speedtest == Inf
        lstring = 'no limit';
    else
        lstring = ['limit = ' num2str(options.speedtest,'%.2g') ' s'];
    end
end

switch lower(state)
    
    case 'seed'
        fprintf(fout, '-----------------------------------------------------------\n');
        fprintf(fout, '%s: Random seed initialized to %d.\n', datestr(now), options.seed);    
    
    case 'startx'
        fprintf(fout, '-----------------------------------------------------------\n');
        fprintf(fout, '%s: Loading starting points for sampling from file.\n', datestr(now));
        
    case 'continue'
        fprintf(fout, '-----------------------------------------------------------\n');
        fprintf(fout, '%s: Continue run and append samples to previous chains.\n', datestr(now));        
        
    case 'speedsuccess'
        testtime = varargin{1};
        fprintf(fout, '%s: Speed test PASSED in %g s (%s).\n', datestr(now), testtime, lstring);
        
    case 'speedfail'
        testtime = varargin{1};
        fprintf(fout, '%s: Speed test FAILED in %g s (%s).\n', datestr(now), testtime, lstring);
                
    case 'hessian'        
        fprintf(fout, '%s: Calculating Hessian for model %s, %s, cnd %s, replica %d.\n', ...
            datestr(now), modelstring, dataidstring, numarray2str(options.cnd), options.replica); 
        
    case 'samplingresults'
        funcalls = varargin{1};
        loglike = varargin{2};
        metrics = varargin{3};
        robusttheta = varargin{4};
        if ~isempty(metrics) && isfield(metrics,'dic') && ~isempty(metrics.dic)
            fprintf(fout, '%s: Sampling: fun calls %d, loglike %g, DIC %g.\n', datestr(now), funcalls, loglike, metrics.dic);            
        else
            fprintf(fout, '%s: Sampling: fun calls %d, loglike %g.\n', datestr(now), funcalls, loglike);
        end
        fprintf(fout, '\tRobust mean params: %s.\n', numarray2str(robusttheta,'%.4g'));        
end

end

%--------------------------------------------------------------------------
function ModelWork_modelFitOne()
%MODELWORK_MODELFITONE Dummy function for retro-compatibility.
    function mlogpdf()
    end
end

function cleanmodelfit()
%CLEANMODELFIT Dummy function for retro-compatibility.

end

