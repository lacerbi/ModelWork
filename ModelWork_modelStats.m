function [mfit] = ModelWork_modelStats(project,mfit,hessianflag,nnoise,force,options)
% MODELWORK_MODELSTATS compute statistics of a model fit.
%   MFIT = MODELWORK_MODELSTATS(PROJECT,MFIT) updates a series of 
%   statistics and model metrics for project PROJECT and model fit MFIT.
%
%   MFIT = MODELWORK_MODELSTATS(PROJECT,MFIT,1) also computes the Hessian
%   at the MAP solution (default is 0).
%
%   MFIT = MODELWORK_MODELSTATS(PROJECT,MFIT,HESSFLAG,NNOISE) also computes 
%   an estimate of the noise in the log likelihood by taking NNOISE samples
%   (default is 0).
%
%   MFIT = MODELWORK_MODELSTATS(PROJECT,MFIT,HESSFLAG,NNOISE,1) forces
%   a recalculation of all statistics (by default, statistics are 
%   recomputed only if there was a change).

% By Luigi Acerbi <luigi.acerbi@gmail.com>
% Last update: Feb/03/2015

if nargin<3 || isempty(hessianflag); hessianflag = 0; end
if nargin<4 || isempty(nnoise); nnoise = 0; end
if nargin<5 || isempty(force); force = 0; end
if nargin<6; options = []; end

% Check if an update of model statistics is needed
if hessianflag && (~isfield(mfit.metrics,'H') || isempty(mfit.metrics.H)); mfit.uptodate = 0; end
if nnoise > 0 && (~isfield(mfit.metrics,'maploglikeSD') || isnan(mfit.metrics.maploglikeSD)); mfit.uptodate = 0; end
if force; mfit.uptodate = 0; end
if mfit.uptodate; return; end
    
% Clear memory first
clear functions;

% Default options for ModelStats
if isempty(options)
    options.maxstoredsamples = [];  % Maximum number of stored samples
    options.recomputesamplingmetrics = [];      % Force recomputation of log likes for LOO-CV
    options.computemarginallike = [];   % Compute marginal likelihood from MCMC output
end
if isempty(options.maxstoredsamples); options.maxstoredsamples = Inf; end
if isempty(options.recomputesamplingmetrics); options.recomputesamplingmetrics = 0; end
if isempty(options.computemarginallike); options.computemarginallike = 0; end

if isempty(mfit.mp) % If fields are empty, try to load data
    try
        mfit = ModelWork_loadFields(project,mfit);    
        warning(['MFIT has empty fields. Loading data from file ''' mfit.datafile '''.']);
    catch
        error(['MFIT has empty fields and cannot load data from file ''' mfit.datafile '''.']);
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAMPLING DIAGNOSTICS AND PARAMETER ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(mfit,'sampling') && ~isempty(mfit.sampling) && isfield(mfit.sampling,'loglikes')
    mfit.sampling.logliks = mfit.sampling.loglikes;
    mfit.sampling = rmfield(mfit.sampling,'loglikes');
    fprintf('Found obsolete field LOGLIKES in SAMPLING. Correcting to new format.\n');
end
    
samplingflag = ~isempty(mfit.sampling) && ...
    isfield(mfit.sampling,'samples') && size(mfit.sampling.samples,1) > 0;

samplingmetricsflag = samplingflag && ...
    ( size(mfit.sampling.logliks, 2) > 1 || options.recomputesamplingmetrics );

% EXPERIMENTAL FEATURE (do not fully trust results)
marginallikeflag = samplingflag && options.computemarginallike;

% Minus log likelihood function handle
nllfun = @(x,mp,logpriorflag,randomizeflag) ModelWork_like(project,mfit.X,mp,mfit.infostruct,x,logpriorflag,0,randomizeflag);

% Copy number of data points
if isfield(mfit,'nData'); mfit.mp.nData = mfit.nData; end

% Sampling
if samplingflag   
    sampling = mfit.sampling;

    % When sampling you need to return and store for each sample both the 
    % log likelihood and the log prior, and do stuff here correctly
    
    % Compute sampled chains summary statistics
    if ~isfield(mfit, 'sumstats') || isempty(mfit.sumstats)
        sampling.sumstats = samplingStats(sampling);
    end
    
    % Check if log prior is present
    if isempty(sampling.logpriors); logpriors = 0; else logpriors = sampling.logpriors; end
    
    % Update MAP if sampling found a better log likelihood value
    logpost = sum(sampling.logliks, 2) + logpriors;
    [maxlogpost,idx] = max(logpost,[],1);
    if maxlogpost > mfit.metrics.maploglike
        mfit.maptheta = sampling.samples(idx,:); 
        mfit.metrics.maploglike = maxlogpost;
    end
    
    % Diagnostics on sampled chains (Gelman and Rubin's scale statistic R)
    % It assumes all chains have the same length
    if sampling.nchains > 1
        nsamplesperchain = floor(size(sampling.samples, 1)/sampling.nchains);
        X = zeros(nsamplesperchain, size(sampling.samples, 2), sampling.nchains);
        for k = 1:sampling.nchains; X(:,:,k) = sampling.samples((1:nsamplesperchain) + (k-1)*nsamplesperchain, :); end
    else
        % If there is only one chain, split it in two
        halflen = floor(size(sampling.samples,1)/2);
        X(:,:,1) = sampling.samples(1:halflen,:);
        X(:,:,2) = sampling.samples((1:halflen) + halflen,:);
    end
    if size(X, 1) > 5
        [R,neff,~,~,~,tau] = psrf(X);
    else
        nvars = size(sampling.samples,2);
        R = NaN(1,nvars);
        neff = NaN(1,nvars);
        tau = NaN(1,nvars);
    end
    sampling.sumstats.R = R;
    sampling.sumstats.neff = neff;
    sampling.sumstats.tau = tau;
    
    % Multivariate effective sample size (mESS) analysis
    nvars = size(sampling.samples, 2);
    if size(sampling.samples, 1) > sampling.nchains*nvars
        X = [];
        nsamplesperchain = floor(size(sampling.samples, 1)/sampling.nchains);
        for k = 1:sampling.nchains
            X{k} = sampling.samples((1:nsamplesperchain) + (k-1)*nsamplesperchain, :);
        end
        mESS = multiESS(X,[],'lESS',10,200);
        sampling.sumstats.mESS_chain = mESS;
    end    
    
    % Compute mean parameter estimates
    sampling.meantheta = mean(sampling.samples,1);
    sampling.robusttheta = trimmean(sampling.samples, 20, 1);
    
    % Compute DIC
    pointloglike = -nllfun(sampling.meantheta,mfit.mp,0,0);
    mfit.metrics.dic = -4*sum(sampling.sumstats.loglikesmean1, 2) + 2*pointloglike;    
    
    % Compute WAIC and PSIS-LOO cross validation score
    if samplingmetricsflag
        [sampling,mfit.metrics] = trialSampleMetrics(sampling,project,mfit,options);
    end
    
    % Approximate marginal likelihood via weighted harmonic mean (EXPERIMENTAL)
    if marginallikeflag
        [sampling,mfit.metrics] = computeMarginalLike(sampling,project,mfit,options);
    end

    % Do not keep full log data table
    sampling.logliks = sum(sampling.logliks, 2);

    % Thin data if samples are more than allowed number of stored samples
    if size(sampling.samples, 1) > options.maxstoredsamples
        storedchains(1:sampling.nchains-1) = floor(options.maxstoredsamples/sampling.nchains);
        storedchains(sampling.nchains) = options.maxstoredsamples - sum(storedchains(1:sampling.nchains-1));

        for k = 1:sampling.nchains
            nsamplesperchain = floor(size(sampling.samples, 1)/sampling.nchains);
            idx = (1:nsamplesperchain) + (k-1)*nsamplesperchain;
            thin = round(linspace(nsamplesperchain/storedchains(k), nsamplesperchain, storedchains(k)));
            idx_thin = (1:storedchains(k)) + sum(storedchains(1:k-1));

            chainsmpl = sampling.samples(idx, :);
            newsmpl(idx_thin, :) = chainsmpl(thin, :);

            chainloglikes = sampling.logliks(idx, :);
            newloglikes(idx_thin, :) = chainloglikes(thin, :);

            if ~isempty(sampling.logpriors)
                chainlogpriors = sampling.logpriors(idx, :);
                newlogpriors(idx_thin, :) = chainlogpriors(thin, :);
            end
        end

        sampling.samples = newsmpl;
        sampling.logliks = newloglikes;
        if ~isempty(sampling.logpriors); sampling.logpriors = newlogpriors; end
    end
    
    mfit.sampling = sampling;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL FITTING AND MODEL COMPARISON
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Save model parameter structure
setupModelFun = str2func([project '_setupModel']);

% Temporarily update the model parameter structure to MAP
mfit.mp.computation = 'precise';
[mfit.mp, exitflag] = setupModelFun(mfit.mp, mfit.maptheta);

% Compute MLE/MAP and Hessian
mfit.metrics.maploglike = -nllfun(mfit.maptheta',mfit.mp,1,0);
% Compute marginal likelihood and Hessian
mfit.metrics.mapmarginallike = mfit.metrics.maploglike;

% Rough estimate of noise in the computation of the log likelihood
if nnoise > 0
    ll = zeros(1,nnoise);
    for i = 1:nnoise
        clear functions;
        ll(i) = -nllfun(mfit.maptheta',mfit.mp,0,1);
    end
    mfit.metrics.maploglikeSD = std(ll);
else
    mfit.metrics.maploglikeSD = NaN;
end    

% Compute BIC, AIC and AICc (corrected AIC)
n = mfit.nData;
k = length(mfit.maptheta);
mfit.metrics.bic = -2 * mfit.metrics.maploglike + k * log(n);
mfit.metrics.aic = -2 * mfit.metrics.maploglike + k * 2;
mfit.metrics.aicc = mfit.metrics.aic + 2 * k * (k + 1) / (n - k - 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTATION OF THE HESSIAN AND MARGINAL LIKELIHOOD (VIA LAPLACE'S METHOD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute Hessian if requested (necessary for marginal likelihood)
if hessianflag    
    mfit.mp.computation = 'hessian';
    [mfit.mp,exitflag] = setupModelFun(mfit.mp, mfit.maptheta);
    mlogpdfhess = @(x) ModelWork_like(project,mfit.X,mfit.mp,mfit.infostruct,x,1,[],0);
    % Prevent computations exactly on a parameter boundary
    dtheta = max(mfit.mp.bounds.SCALE/10000, 1e-4);
    hesstheta = min(max(mfit.maptheta, mfit.mp.bounds.LB + dtheta), mfit.mp.bounds.UB - dtheta);         
    [mfit.metrics.H, mfit.metrics.Herr] = hessian(mlogpdfhess, hesstheta');
else
    mfit.metrics.H = []; mfit.metrics.Herr = [];    
end

% Calculate marginal likelihood via Laplace (if Hessian is available)
if isempty(mfit.metrics.H)
    mfit.metrics.marginallike = [];
else
    mfit.metrics.marginallike = mfit.metrics.mapmarginallike + 0.5*length(mfit.maptheta)*log(2*pi) - 0.5*log(abs(det(mfit.metrics.H)));
end

% Order fields
if ~isempty(mfit.metrics); mfit.metrics = orderfields(mfit.metrics); end
if ~isempty(mfit.sampling); mfit.sampling = orderfields(mfit.sampling); end

mfit.mp.computation = []; % Go back to normal

% Update model parameter structure to posterior mean if sampling (deprecated)
%if samplingflag
%    [mfit.mp,exitflag] = setupModelFun(mfit.mp, mfit.sampling.meantheta);
%end

mfit.uptodate = 1;

end

%--------------------------------------------------------------------------
function [sampling,metrics] = trialSampleMetrics(sampling,project,mfit,options)
%TRIALSAMPLEMETRICS Compute trial-based sampling metrics (WAIC, PSIS-LOO CV).

metrics = mfit.metrics;
extras = [];

% If needed, recompute the per-trial log likelihood
if size(sampling.logliks,2) == 1
    nllfun = @(th_) ModelWork_like(project,mfit.X,mfit.mp,mfit.infostruct,th_,0,1,0);
    [nlls(1,:),extras] = nllfun(sampling.samples(1,:));
    if size(nlls,2) > 1
        startmsg = 'Recomputing trial-based log likelihoods... Sample ';
        filename = [options.fullfilename(1:end-4) '_batch.tmp'];        
        % Jitter save to avoid simultaneous save from all processes
        SaveTime = (0.9 + 0.2*rand())*options.savetime;
        temp = fbatcheval(nllfun,sampling.samples(2:end,:),1,filename,SaveTime,startmsg);
        nlls(2:size(sampling.samples,1),:) = temp;        
        % Remove empty trials
        idx = all(nlls == 0,1);
        nlls(:,idx) = [];        
        sampling.logliks = -nlls;
        if size(sampling.logliks,2) ~= mfit.nData
            warning(['Declared number of trials in model fit (n=' num2str(mfit.nData) ') does not match log likelihood (n=' num2str(size(sampling.logliks,2)) ').']);
        end
    end
end

% Compute trial-based metrics if trial log-likelihoods are available
if size(sampling.logliks,2) > 1
    if isfield(options,'binnedloglik') && ~isempty(options.binnedloglik) ...
            && options.binnedloglik
        % Log likelihood is binned, need to unpack
        if isempty(extras); [~,extras] = nllfun(sampling.samples(1,:)); end
        error('Unpacking of log likelihood not supported.');
    else
        logliks = sampling.logliks;
    end
    
    % Recompute summary statistics
    sampling.sumstats = samplingStats(sampling);
    
    % Compute WAIC1 and 2 (see Bayesian Data Analysis, 3rd edition)
    waic1 = 2*sum(log(sampling.sumstats.loglikesmean1exp) - sampling.sumstats.loglikesmean1, 2);
    sum1sq = sampling.sumstats.loglikesmean1.^2.*sampling.sumstats.n;
    waic2 = sum((sampling.sumstats.loglikessqsum1 - sum1sq)./(sampling.sumstats.n-1), 2);
    lppd = sum(log(sampling.sumstats.loglikesmean1exp), 2);
    mfit.metrics.waic1 = -2*lppd + 2*waic1;
    mfit.metrics.waic2 = -2*lppd + 2*waic2;

    [loo,loos,ks] = psisloo(logliks);
    metrics.loocv = loo;        
    sampling.sumstats.loos = loos;
    sampling.sumstats.ks = ks;
end

end
%--------------------------------------------------------------------------
function [sampling,metrics] = computeMarginalLike(sampling,project,mfit,options)

metrics = mfit.metrics;

vboptions.Nstarts = 3;
vboptions.Display = 'off';
X = sampling.samples';
Y = sum(sampling.logliks,2);

logpriors = sampling.logpriors;
% If all log priors are empty or zero, assume uniform priors
if isempty(logpriors) || all(logpriors == 0)
    logpriors = -sum(log(mfit.mp.bounds.UB - mfit.mp.bounds.LB));
end
Y = Y + logpriors;

vbmodel = vbgmmfit(X,[],[],vboptions);

% Marginal likelihood approximation methods
methods = {'whmg','whmu','rlr'};

for m = 1:numel(methods)    
    [logZ,slogZ,fracZ,neff] = vbgmmmarglike(vbmodel,X,Y,methods{m});
    % Pick estimate with lowest uncertainty
    [slogZ,idx] = min(slogZ);
    logZ = logZ(idx);
    sampling.marginallike.(methods{m}).logZ = logZ;
    sampling.marginallike.(methods{m}).slogZ = slogZ;
    sampling.marginallike.(methods{m}).fracZ = fracZ;
    sampling.marginallike.(methods{m}).neff = neff;
    metrics.(['marginallike_' methods{m}]) = logZ;
end

% This occupies a lot of memory, do not save
% sampling.vbmodel = vbmodel;

end

%--------------------------------------------------------------------------
function sumstats = samplingStats(sampling)
%SAMPLINGSTATS Compute sampling summary statistics.

if isfield(sampling,'sumstats'); sumstats = sampling.sumstats; end
logliks = sum(sampling.logliks,2);
lz = max(logliks);
sumstats.smplmean1 = nanmean(sampling.samples, 1);
sumstats.loglikesmean1 = nanmean(logliks, 1);
sumstats.loglikesmean1exp = exp(lz)*nanmean(exp(logliks-lz), 1);
sumstats.loglikessqsum1 = nansum(logliks.^2, 1);
sumstats.n = sum(~isnan(logliks), 1);

end