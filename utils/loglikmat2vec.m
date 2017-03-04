function ll_vec = loglikmat2vec(ll_mat,X)
%LOGLIKMAT2VEC Convert a log likelihood matrix to a array of trials.
%   LL_VEC = LOGLIKMAT2VEC(LL_MAT,X) converts log likelihood matrix LL_MAT
%   into an array of log likelihoods per trial given data matrix X.
%   LL_MAT and X are matrices with the same dimensions. Each element of X
%   specifies the number of counts within each bin or trial type.

if any(size(ll_mat) ~= size(X))
    error('Matrices LL_MAT and X should have the same size.');
end

ll_mat = ll_mat(:);
X = X(:);

nTrials = sum(X);   % Total number of trials
ll_vec = zeros(nTrials,1);

% Loop over nonempty trials
idx = 1;
for i = find(X)'
    n = X(i);
    ll_vec(idx:idx+n-1) = ll_mat(i);
    idx = idx + n;
end

end