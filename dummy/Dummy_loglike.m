% DUMMY_DATALIKE compute minus log likelihood of a dataset.
% 
%   LOGLIKE = VESTBMS_DATALIKE([],MP,INFOSTRUCT) returns log likelihood 
%   of model-parameter struct MP with extra info INFOSTRUCT.
%
%   [LOGLIKE,EXTRAS] = VESTBMS_DATALIKE(...) also returns a series of 
%   additional structures in EXTRAS.
%
%  By Luigi Acerbi <luigi.acerbi@gmail.com>

function [loglike,extras] = Dummy_loglike(X,mp,infostruct,flags)

if nargin < 4 || isempty(flags); flags = [0 0]; end
if size(flags,2) < 2; flags(2) = 0; end

debug = flags(1);
randomize = flags(2);

extras = [];

% Persistent variables
% persistent oldparams;
% persistent oldloglikes;
% persistent oldtrialloglikes;
% persistent starttrials;
% if isempty(oldparams)
%     % Prepare variables to store parameters between function calls
%     for iicnd = 1:4; oldparams{iicnd} = zeros(1, 10); end
%     for iicnd = 5:7; oldparams{iicnd} = zeros(1, 20); end
%     oldloglikes = zeros(1, 7);
%     
%     % Count expected number of trial types per condition
%     ntrials = zeros(1, 7);
%     for iicnd = 1:4; ntrials(iicnd) = numel(X.unibins{iicnd}(:)); end
%     for iicnd = 5:7
%         for iTask = 1:numel(X.bimbins{iicnd-4})
%             ntrials(iicnd) = ntrials(iicnd) + numel(X.bimbins{iicnd-4}{iTask}(:));
%         end
%     end
%     starttrials = [0 cumsum(ntrials)]+1;    
%     oldtrialloglikes = NaN(1, sum(ntrials));
% end

cnd = mp.cnd;
model = mp.model;

loglikes = zeros(1, max(cnd));

for iicnd = 1:length(cnd)        
    fulltheta = mp.fulltheta{iicnd};   
    
    mu = fulltheta.mu;
    sigma = fulltheta.sigma;
    lambda = fulltheta.lambda;
    
    % Compute psychometric function
    like = psychomodel(mu,sigma,lambda,X(:,1),X(:,2));    
    trialloglikes = log(like);
    loglikes(iicnd) = sum(trialloglikes);
end

%oldloglikes = loglikes;
%oldtrialloglikes = trialloglikes;

if any(isnan(loglikes) | isinf(loglikes) | ~isreal(loglikes))
    loglikes
    pause
end

loglike = sum(loglikes);
if debug; extras.struct = extras; end
if ~isempty(trialloglikes); extras.trialloglikes = trialloglikes; end

end