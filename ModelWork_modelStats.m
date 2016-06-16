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
    options.recomputeloo = [];      % Force recomputation of log likes for LOO-CV
end
if isempty(options.maxstoredsamples); options.maxstoredsamples = Inf; end
if isempty(options.recomputeloo); options.recomputeloo = 0; end

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

samplingflag = ~isempty(mfit.sampling) && ...
    isfield(mfit.sampling,'samples') && size(mfit.sampling.samples,1) > 0;

loocvflag = samplingflag && ...
    ( size(mfit.sampling.logliks, 2) > 1 || options.recomputeloo );

% Minus log likelihood function handle
nllpdf = @(x,mp,logpriorflag,randomizeflag) ModelWork_like(project,mfit.X,mp,mfit.infostruct,x,logpriorflag,0,randomizeflag);

% Copy number of data points
if isfield(mfit,'nData'); mfit.mp.nData = mfit.nData; end

% Sampling
if samplingflag   
    sampling = mfit.sampling;

    % When sampling you need to return and store for each sample both the 
    % log likelihood and the log prior, and do stuff here correctly
    
    % Compute sampled chains summary statistics
    if ~isfield(mfit, 'sumstats') || isempty(mfit.sumstats)
        logliks = sum(sampling.logliks,2);
        lz = max(logliks);
        sampling.sumstats.smplmean1 = nanmean(sampling.samples, 1);
        sampling.sumstats.loglikesmean1 = nanmean(logliks, 1);
        sampling.sumstats.loglikesmean1exp = exp(lz)*nanmean(exp(logliks-lz), 1);
        sampling.sumstats.loglikessqsum1 = nansum(logliks.^2, 1);
        sampling.sumstats.n = sum(~isnan(logliks), 1);
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
    
    % Compute mean parameter estimates
    sampling.meantheta = mean(sampling.samples,1);
    sampling.robusttheta = trimmean(sampling.samples, 20, 1);
    
    % Compute DIC and WAIC
    pointloglike = -nllpdf(sampling.meantheta,mfit.mp,0,0);
    mfit.metrics.dic = -4*sum(sampling.sumstats.loglikesmean1, 2) + 2*pointloglike;    

    % Compute WAIC1 and 2 (see Bayesian Data Analysis, 3rd edition)
    waic1 = 2*sum(log(sampling.sumstats.loglikesmean1exp) - sampling.sumstats.loglikesmean1, 2);
    sum1sq = sampling.sumstats.loglikesmean1.^2.*sampling.sumstats.n;
    waic2 = sum((sampling.sumstats.loglikessqsum1 - sum1sq)./(sampling.sumstats.n-1), 2);
    lppd = sum(log(sampling.sumstats.loglikesmean1exp), 2);
    mfit.metrics.waic1 = -2*lppd + 2*waic1;
    mfit.metrics.waic2 = -2*lppd + 2*waic2;
    
    % Compute PSIS-LOO
    if loocvflag
        % If needed, recompute the per-trial log likelihood
        if size(sampling.logliks,2) == 1
            loocvfun = @(th_) ModelWork_like(project,mfit.X,mfit.mp,mfit.infostruct,th_,0,1,0);
            logliks(1,:) = loocvfun(sampling.samples(1,:));
            if size(logliks,2) > 1
                logliks = repmat(logliks, [size(sampling.samples,1),1]);        
                for iSamp = 2:size(sampling.samples,1)
                    logliks(iSamp,:) = loocvfun(sampling.samples(iSamp,:));
                end
                sampling.logliks = logliks;
            end
        end
        
        % Compute PSIS-LOO CV score if trial log-likelihoods are available
        if size(sampling.logliks,2) > 1
            [loo,loos,ks] = psisloo(sampling.logliks);
            mfit.metrics.loocv = loo;        
            sampling.sumstats.loos = loos;
            sampling.sumstats.ks = ks;
        end
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
mfit.metrics.maploglike = -nllpdf(mfit.maptheta',mfit.mp,1,0);
% Compute marginal likelihood and Hessian
mfit.metrics.mapmarginallike = mfit.metrics.maploglike;

% Rough estimate of noise in the computation of the log likelihood
if nnoise > 0
    ll = zeros(1,nnoise);
    for i = 1:nnoise
        clear functions;
        ll(i) = -nllpdf(mfit.maptheta',mfit.mp,0,1);
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

% Calculate marginal likelihood (if Hessian is available)
if isempty(mfit.metrics.H)
    mfit.metrics.marginallike = [];
else
    mfit.metrics.marginallike = mfit.metrics.mapmarginallike + 0.5*length(mfit.maptheta)*log(2*pi) - 0.5*log(abs(det(mfit.metrics.H)));
end

% Order fields
if ~isempty(mfit.metrics); mfit.metrics = orderfields(mfit.metrics); end
if ~isempty(mfit.sampling); mfit.sampling = orderfields(mfit.sampling); end

mfit.mp.computation = []; % Go back to normal

% If sampling, update model parameter structure to the posterior mean
if samplingflag
    [mfit.mp,exitflag] = setupModelFun(mfit.mp, mfit.sampling.meantheta);
end

mfit.uptodate = 1;

end