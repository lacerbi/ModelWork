% DUMMY_CREATEDATASET
function data = Dummy_createDataset()

data = [];
nSubjs = 10;
nTrials = 250;

xbounds = [-3,3];

mu = randn(1,nSubjs);
sigma = exp(randn(1,nSubjs));
lambda = -0.05*log(rand(1,nSubjs)); % Exponential

for i = 1:nSubjs
    data{i}.mu = mu(i);
    data{i}.sigma = sigma(i);
    data{i}.lambda = lambda(i);
    
    X = rand(nTrials,1)*diff(xbounds) + xbounds(1);
    
    [~,pright] = psychomodel(mu(i),sigma(i),lambda(i),X,zeros(nTrials,1));
    R = rand(nTrials,1) < pright;
    
    data{i}.X = [X, R];
end

save('Dummy_data.mat','data');

end