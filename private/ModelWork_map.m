function [x,fval,exitflag,output] = ModelWork_map(mfit,mp,infostruct,fout,datalikeFun,setupModelFun,logpriorFun,firstSetting,options)
%MODELWORK_MAP Model fitting via maximum-a-posteriori (MAP) estimation.
% (MLE estimation is same as MAP with a uniform prior on all parameters)
%

% Import parameters lower and upper bounds
LB = mp.bounds.LB; UB = mp.bounds.UB;
RLB = mp.bounds.RLB; RUB = mp.bounds.RUB;
nvars = size(LB,2);

if ~isempty(fout); displ = 'iter'; else displ = 'off'; end

% Apply parameter prior?
logpriorflag = options.logpriorflag;

nData = mfit.nData;
funccount = 0;

% Initialize direct search (use fminsearch for low dimensions)
if isempty(options.optimizationmethod)            
    if nvars <= 5; options.optimizationmethod = 'fminsearch';
    else options.optimizationmethod = 'fmincon'; end
end
if strcmp(options.optimizationmethod,'fmincon'); localoptions.Algorithm = 'sqp';
else localoptions = []; end

% Default parameters for patternsearch
psoptions.Cache = 'on'; psoptions.CacheTol = 1e-6; 
psoptions.InitialMeshSize = 0.125; psoptions.TolMesh = 1e-5;

MaxFunEvals = options.optfevals;
if isempty(options.optrefinefevals); RefineFunEvals = MaxFunEvals;
else RefineFunEvals = options.optrefinefevals; end
RefineFunEvals = min(RefineFunEvals,100*nvars); % Fcn evals for refinement

% Number of starting points
if options.nstarts == 0
    nstarts = 1;
else
    if isempty(options.nsobol); options.nsobol = options.nstarts*20; end
    nstarts = max(1,[options.nsobol,options.nstarts,round(options.nstarts/10),1]);
end

% Multi-start optimization options
optoptions.InitRange = [RLB; RUB];
optoptions.InitialPoints = options.startx;
optoptions.MidpointStart = firstSetting;
optoptions.SobolInit = 'on';
optoptions.SobolSeed = (options.replica-1)*nstarts(1);
optoptions.Display = displ;
optoptions.MaxFunEvals = [1,MaxFunEvals,RefineFunEvals,RefineFunEvals];
optoptions.MaxIter = MaxFunEvals;
optoptions.TolFun = 1e-4;
optoptions.OptimOptions = {[], localoptions, psoptions, localoptions};
% optoptions.XScale = mp.bounds.SCALE*2;
optoptions.FvalScale = 1;
optoptions.RescaleVars = 'off';
optoptions.Method = {'feval',options.optimizationmethod,'patternsearch',options.optimizationmethod};
optoptions.BPSUseCacheEpochs = 0;

% Save information
optfilename = [options.fullfilename(1:end-4) '_opt.tmp'];
optoptions.LoadFile = optfilename;
optoptions.SaveFile = optfilename;
% Save every ~1h30, jitter to avoid simultaneous save from all processes
optoptions.SaveTime = (0.9 + 0.2*rand())*options.savetime;

optoptions.OutputFcn = @(x,optimValues,state) outputFcn(x,optimValues,state,fout,options.modelstring,options.dataidstring,options.cnd,options.replica,firstSetting);
vararginarray = {{'coarse',logpriorflag},{'normal',logpriorflag},{'normal',logpriorflag},{'precise',logpriorflag}};

% Get trial masks if present
trainmask = [];
testmask = [];
if isfield(mfit.X,'foldmasks') && ~isempty(mfit.X.foldmasks) && isfield(mfit.X.foldmasks,'train')
    trainmask = mfit.X.foldmasks.train;
    % Add prior to training set
    if logpriorflag; trainmask = [trainmask, ones(size(trainmask,1),1)]; end
end
if isfield(mfit.X,'foldmasks') && ~isempty(mfit.X.foldmasks) && isfield(mfit.X.foldmasks,'test')
    testmask = mfit.X.foldmasks.test;
    % Exclude prior from test set
    if ~isempty(testmask) && logpriorflag; testmask = [testmask, zeros(size(testmask,1),1)]; end
end

if ~isempty(trainmask)
    fminfoldfun([],@(theta,precision,logpriorflag) nLogPosterior(theta,precision,logpriorflag,0),trainmask,testmask,1);
    optfun = @fminfoldfun;
else
    optfun = @(theta,precision,logpriorflag) nLogPosterior(theta,precision,logpriorflag,1);
end

% Run multi-start optimization
[~,~,exitflag,output] = fminmulti(optfun,LB,UB,nstarts,optoptions,vararginarray);    

if ~isempty(trainmask)
    tolfun = 0.1;
    optimizer = @(fun,x0) bads(fun,x0,LB,UB,RLB,RUB);    
    [fold.xs,fold.fvals,fold.ftests,fold.exitflag,fold.output] = fminfoldrun(optimizer,[],tolfun,'precise',logpriorflag);
    output.funccount = sum(fold.output.funcCount);    
else
    output.funccount = sum(output.funcalls);    
end


% Remove coarsely evaluated points
idx = 1:nstarts(1);
output.x(idx,:) = [];
output.startx(idx,:) = [];
output.nruns = output.nruns - nstarts(1);
output.fval(idx) = [];
output.funcalls(idx) = [];
output.exitflag(idx) = [];

% Keep best result with maximum precision (last of array)    
output.xmin = output.x(end,:);
output.fvalmin = output.fval(end);

% See if there is a better minimum for the full fold
if ~isempty(trainmask)
    idx_full = find(all(trainmask,2),1);
    if ~isempty(idx_full)
        if fold.fvals(idx_full) < output.fvalmin
            output.fvalmin = fold.fvals(idx_full);
            output.xmin = fold.xs(idx_full,:);
        end
    end
    output.fold = fold;
else
    output.fold = [];
end

x = output.xmin;
fval = output.fvalmin;
output = rmfield(output,'output');  % Save some memory

    %----------------------------------------------------------------------
    function nLL = nLogPosterior(theta,precision,logpriorflag,sumflag)
    % NLOGPOSTERIOR Negative (unnormalized) log posterior of the parameters
        
        %if nargin < 2; precision = []; end
        %if nargin < 3; logpriorflag = []; end
        
        if nargin < 4 || isempty(sumflag); sumflag = true; end
        
        INFPENALTY = -log(1e-6);

        theta = theta(:)';

        % Update the number of function calls
        funccount = funccount + 1;

        % Update computational precision
        if ~isempty(precision); mp.computation = precision; end
                
        % Test boundary consistency for parameters
        if any(theta < LB | theta > UB); nLL = INFPENALTY*nData; return; end
        
        % Set global model parameters and check their consistency
        [mp,exitflag] = setupModelFun(mp,theta);
        if exitflag; nLL = INFPENALTY*nData; return; end

        [loglike,extras] = datalikeFun(mfit.X,mp,infostruct);            
        mp.extras = extras;

        % Apply prior (use flat prior to yield a ML solution)
        if logpriorflag
            logPrior = logpriorFun(mp,infostruct);
        else
            logPrior = 0;
        end
        
        if sumflag
            nLL = -sum(loglike) - logPrior;
        else
            nLL = [-extras.trialloglikes(:)', -logPrior];
        end
    end
end

%--------------------------------------------------------------------------
function stop = outputFcn(x,optimValues,state,fout,modelstring,dataidstring,cnd,replica,firstSetting)
%OUTPUTFCN Log optimization progress to file

switch state
    case 'init'
        if firstSetting
            fprintf(fout, '%s: Starting iter 1 from non-random reasonable x0.\n', datestr(now));
        end
        
    case 'iter'
        if optimValues.iteration == 0 && strcmp(optimValues.procedure,'feval')
            fprintf(fout, '%s: Calculating initial points for epoch 1, model %s, %s, cnd %s, replica %d.\n', ...
                datestr(now), modelstring, dataidstring, numarray2str(cnd,'%d'),replica);
        elseif optimValues.iteration > 0 && ~strcmp(optimValues.procedure,'feval')
            if optimValues.funccount == 0
                fprintf(fout, '%s: Calculating epoch %d, iter %d, for model %s, %s, cnd %s, replica %d.\n', ...
                    datestr(now), optimValues.epoch, optimValues.iteration, modelstring, dataidstring, numarray2str(cnd,'%d'),replica);
                fprintf(fout, '\tStarting params: %s\n', numarray2str(x,'%.4g'));
            else
                fprintf(fout, '%s: Epoch %d, iter %d: fun calls %d, loglike %g.\n', ...
                    datestr(now), optimValues.epoch, optimValues.iteration, optimValues.funccount, -optimValues.fval);
                fprintf(fout, '\tMAP params: %s\n', numarray2str(x,'%.4g'));
            end            
        end
                
    case 'done'
        fprintf(fout, '%s: Optimization result: loglike %g.\n', datestr(now), -optimValues.fval);
        fprintf(fout, '\tMAP params: %s\n', numarray2str(x,'%.4g'));        
end

stop = 0;

end