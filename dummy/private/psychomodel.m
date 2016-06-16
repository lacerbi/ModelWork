function [like,pright] = psychomodel(mu,sigma,lambda,X,R)
%PSYCHOMODEL Standard psychometric function for left/right.

MIN_P = 1e-6;

z = (X - mu)./sigma;
respRight = R == 1;
respLeft = R == 0;

pright = 0.5*erfc(-z/sqrt(2));

% Add lapse
if lambda > 0
    pright =  lambda/2 + (1 - lambda) .* pright;
end

% Compute likelihood
like = pright.*respRight + (1-pright).*respLeft;
like = (1-MIN_P) * like + MIN_P;

end