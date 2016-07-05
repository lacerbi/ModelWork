function loglikvec = loglikmat2vec(llmat,X)
%LOGLIKMAT2VEC Convert a log likelihood matrix to a vector of trials.

llmat = llmat(:);
X = X(:);

nTrials = sum(X);
loglikvec = zeros(nTrials,1);

% Loop over nonempty trials
idx = 1;
for i = find(X)'
    n = X(i);
    loglikvec(idx:idx+n-1) = llmat(i);
    idx = idx + n;
end

end