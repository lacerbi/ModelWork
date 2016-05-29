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

% Multi-start optimization options
optoptions.InitRange = [RLB; RUB];
optoptions.InitialPoints = options.startx;
optoptions.MidpointStart = firstSetting;
optoptions.Display = displ;
optoptions.MaxFunEvals = [1,MaxFunEvals,RefineFunEvals,RefineFunEvals];
optoptions.MaxIter = MaxFunEvals;
optoptions.TolFun = 1e-4;
optoptions.OptimOptions = {[], localoptions, psoptions, localoptions};
% optoptions.XScale = mp.bounds.SCALE*2;
optoptions.FvalScale = 1;
optoptions.RescaleVars = 'off';
optoptions.Method = {'feval',options.optimizationmethod,'patternsearch',options.optimizationmethod};
optoptions.OutputFcn = @(x,optimValues,state) outputFcn(x,optimValues,state,fout,options.modelstring,options.dataidstring,options.cnd,options.replica,firstSetting);
if options.nstarts == 0
    nstarts = 1;
else
    nstarts = max(1,[options.nstarts*20,options.nstarts,round(options.nstarts/10),1]);
end
vararginarray = {{'coarse',logpriorflag},{'normal',logpriorflag},{'normal',logpriorflag},{'precise',logpriorflag}};

% Run multi-start optimization
[~,~,exitflag,output] = fminmulti(@nLogPosterior,LB,UB,nstarts,optoptions,vararginarray);

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
x = output.xmin;
fval = output.fvalmin;
output = rmfield(output,'output');  % Save some memory
output.funccount = funccount;

    %----------------------------------------------------------------------
    function nLL = nLogPosterior(theta,precision,logpriorflag)
    % NLOGPOSTERIOR Negative (unnormalized) log posterior of the parameters
        
        %if nargin < 2; precision = []; end
        %if nargin < 3; logpriorflag = []; end
        
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
        nLL = -sum(loglike) - logPrior;
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