function [sampling,exitflag,output] = ModelWork_sample(mfit,mp,infostruct,fout,datalikeFun,setupModelFun,logpriorFun,x0,weights,firstSetting,options)
%MODELWORK_SAMPLE Model fitting via sampling.

% replica = options.replica;

% Import parameters lower and upper bounds
LB = mp.bounds.LB; UB = mp.bounds.UB;
RLB = mp.bounds.RLB; RUB = mp.bounds.RUB;
nvars = size(LB,2);

funccount = 0;

if ~isempty(fout); displ = 'iter'; else displ = 'off'; end

if isfield(mfit,'sampling') && ~isempty(mfit.sampling)
    % Get sampling struct
    sampling = mfit.sampling;
else        
    % Data-dependent sampling information
    sampling.samples = [];
    sampling.loglikes = [];
    sampling.logpriors = [];
    sampling.nchains = 1;
    sampling.funccount = 0;
    sampling.nsamplestot = 0;
    sampling.sumstats = [];
end

% Check that all observables have the same length
len1 = size(sampling.samples,1);
assert((size(sampling.loglikes,1) == len1) && (size(sampling.logpriors,1) == len1 || isempty(sampling.logpriors)), ...
    'Stored samples and log likelihoods/priors do not match in size.');

% Compute thinning based on requested samples and storage capacity
% (there is no reason for thinning aside of space requirements)
if options.nsamples <= options.maxstoredsamples
    nsamples = options.nsamples;
    thin = 1;   % No thinning, take all samples
else
    thin = floor(options.nsamples/options.maxstoredsamples);
    nsamples = ceil(options.nsamples/thin);
end

% Choose starting point for sampling if not specified
if isempty(x0)
    if isempty(options.startx)
        if firstSetting
            x0 = 0.5*(RUB + RLB);
        else
            x0 = RLB + rand(1,nvars).*(RUB-RLB);
        end
    else
        x0 = options.startx;
    end
end

% Write log to file
writelog(fout,'start',options,nsamples,thin,options.maxstoredsamples);

% unknownsamplingmethod = 0;

% try
    % Different sampling methods
    switch lower(options.samplingmethod)

        case 'adss' % Adaptive Directional Slice Sampling

            samplepdf = @(thx) LogPosterior(thx,1,[],0);
            logpriorpdf = @(thx) LogPosterior(thx,1,[],1);
            smploptions.LogPrior = logpriorpdf;
            if ~isfield(sampling,'burnin') || isempty(sampling.burnin); sampling.burnin = nsamples; end
            smploptions.Burnin = sampling.burnin;
            smploptions.Thin = thin;
            if strcmpi(displ,'iter'); smploptions.Display = 'notify';
            else smploptions.Display = 'off'; end
            samplefilename = [options.fullfilename(1:end-4) '_sample.tmp'];
            smploptions.LoadFile = samplefilename;
            smploptions.SaveFile = samplefilename;
            % Save every ~1h30, jitter to avoid simultaneous save from all processes
            smploptions.SaveTime = (0.9 + 0.2*rand())*options.savetime;
            K = 40 + 5*nvars;
            
            % Randomize starting points according to weights
            if ~isempty(weights)
                idx = mnrnd_private(repmat(weights(:)',[K,1]));
                x0 = x0(idx,:);
            end

            [samples,loglikes,exitflag,output] = adss(samplepdf,x0,K,nsamples,mp.bounds.SCALE,LB,UB,smploptions);

            logpriors = output.logpriors;
            funccount = output.funccount;

%             case 'dramsample'
%                 dram.model.ssfun = @(x,d) 2*LogPosterior(x);
%                 dram.options.method = 'dram';
%                 dram.options.waitbar = 0;
%                 dram.options.nsimu = nsamples;
%                 dram.options.burnintime = nburnin;
%                 dram.options.verbosity = 0;
%                 dram.options.ntry = 4;
%                 dram.options.savefile = ['dram-outfile-' num2str(options.processid) '-' num2str(g) '.mat'];
%                 dram.options.saveiter = options.dramsaveiter;
%                 dram.options.runtimemax = options.runtimemax - toc(options.starttime);
%                 if ~isempty(options.mfit) && options.appendchains
%                     options.dramqcov = options.mfit.dramresults.qcov;
%                 end                        
%                 if isfield(options, 'dramqcov') && ~isempty(options.dramqcov)
%                     dramqcov = options.dramqcov;
%                 else
%                     dramqcov = 2.4^2/D*diag((mp.bounds.SCALE).^2);
%                 end                        
%                 dram.options.qcov = dramqcov;
%                 % Set parameters
%                 for ii = 1:D
%                     dram.params{ii} = {mp.params{ii}, ml{g}(ii), mp.bounds.LB(ii), mp.bounds.UB(ii)};
%                 end
%                 timeout = 0;
%                 [results,samples,~,loglikes2,~,timeout] = mcmcrunplus(dram.model, [], dram.params, dram.options);
%                 % In case of timeout, simply quit without saving
%                 if timeout; exitflag = 1; return; end
%                 samples = samples';
%                 loglikes = -loglikes2/2;
% 
%                 % Avoid copying function handles otherwise MATLAB
%                 % behaves in a funny way
%                 results.ssfun = [];
%                 results.priorfun = [];
% 
%                 mfit.dramresults = results;
%                 % results
% 
%             case 'kebabsample'
%                 [samples loglikes] = kebab_sample(nsamples, nburnin, 0, @(thx) LogPosterior(thx,1), ml{g}, mp.bounds.SCALE, 0, 0);
%                 
%             case 'maxsample'
%                 if ~isempty(options.mfit) && options.appendchains
%                     [samples loglikes] = maxsample(nsamples, @(thx) LogPosterior(thx), options.mfit.smpl, mp.bounds.SCALE/D, options.maxsamplealpha, options.mfit.loglikes);
%                     offset = size(options.mfit.smpl,1);
%                     samples = samples(offset + (1:nsamples), :);
%                     loglikes = loglikes(offset + (1:nsamples));                        
%                 else
%                     [samples loglikes] = maxsample(nsamples, @(thx) LogPosterior(thx), ml{g}', mp.bounds.SCALE/D, options.maxsamplealpha);
%                 end
%                 samples = samples';
%                 
%             case 'slicesample'
%                 [samples loglikes] = slice_sample(nsamples, nburnin, 0, @(thx) LogPosterior(thx,1), ml{g}, mp.bounds.SCALE, 0, 0);

        otherwise
            %unknownsamplingmethod = 1;
            error(['Unknown sampling method: ''' options.samplingmethod '''']);

    end

    % Check that returned observables match in length
    len2 = size(samples,1);
    assert((size(loglikes,1) == len2) && (size(logpriors,1) == len2 || isempty(logpriors)), ...
        'Computed samples and log likelihoods/priors do not match in size.');

    % Concatenate samples to existing chain
    n1 = sampling.nsamplestot;
    n2 = options.nsamples;
    nmax = options.maxstoredsamples;

    [sampling.samples,nc] = concatchains(sampling.samples,n1,samples,n2,nmax);
    sampling.loglikes = concatchains(sampling.loglikes,n1,loglikes,n2,nmax);
    sampling.logpriors = concatchains(sampling.logpriors,n1,logpriors,n2,nmax);
    sampling.nsamplestot = nc;
    sampling.funccount = sampling.funccount + funccount;

%catch err
%    
%    writelog(fout,'error',options,err.message);
%    if unknownsamplingmethod
%        error(['Unknown sampling method: ''' options.samplingmethod '''.']); 
%    else
%        error(['Error while sampling: ' err.message]);
%    end
%    
%end

% Write log to file
writelog(fout,'finish',options,funccount);

output.funccount = funccount;   % Return number of function evaluations

return;
    

    %----------------------------------------------------------------------
    function LL = LogPosterior(theta,trialloglikesflag,precision,doprior)
    % LOGPOSTERIOR Unnormalized log posterior of the parameters.

        persistent oldtheta;
    
        if nargin < 2; trialloglikesflag = 0; end
        if nargin < 3; precision = []; end
        if nargin < 4; doprior = 0; end
                
        theta = theta(:)';
        
        % Call with new parameter vector
        if isempty(oldtheta) || any(theta ~= oldtheta)
            oldtheta = theta;

            % Update the number of function calls
            funccount = funccount + 1;

            % Update computational precision
            if ~isempty(precision); mp.computation = precision; end

            % Test boundary consistency for parameters
            if any(theta < LB) || any(theta > UB)
                if trialloglikesflag; LL = -Inf*ones(1,nData);
                else LL = Inf; end
                return;
            end

            % Set global model parameters and check their consistency
            [mp, exitflag] = setupModelFun(mp,theta);
            if exitflag
                if trialloglikesflag; LL = -Inf*ones(1,nData);
                else LL = -Inf*nData; end
                return;
            end
        end
        
        % Output is either log prior or log likelihood
        if doprior
            LL = logpriorFun(mp,infostruct);
        else
            [loglike,extras] = datalikeFun(mfit.X,mp,infostruct);            
            mp.extras = extras;

            if trialloglikesflag && isfield(extras,'trialloglikes')
                LL = extras.trialloglikes;
            else
                LL = loglike;
            end
        end
    end    
end

%--------------------------------------------------------------------------
function writelog(fout,state,options,varargin)
%WRITELOG Write entry in log file.

if isempty(fout); return; end
if ~strcmpi(options.display, 'all'); return; end

switch lower(state)
    
    case 'start'
        nsamples = varargin{1};
        thin = varargin{2};
        nstoredsamples = varargin{3};
        fprintf(fout, '%s: Start sampling (samples=%d,thin=%d,stored=%d) with method ''%s''.\n', ...
            datestr(now), nsamples, thin, nstoredsamples, options.samplingmethod);
        
    case 'finish'
        funccount = varargin{1};
        fprintf(fout, '%s: Finished sampling after %d function evaluations.\n', ...
            datestr(now), funccount);
        
    case 'error'
        message = varargin{1};
        fprintf(fout, '%s: ERROR while sampling: %s.\n', ...
            datestr(now), message);       
        
end

end

%MNRND_PRIVATE Random sample from unnormalized multinomial distribution.
function r = mnrnd_private(p)
    cdf = cumsum(p,2);
    m = size(p,1);
    samp_k = bsxfun(@gt, cdf, rand(m,1).*sum(p,2));
    idx = diff([zeros(m,1), samp_k],[],2);
    [~,r] = max(idx,[],2);
end
