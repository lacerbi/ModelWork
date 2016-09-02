function ModelWork_loglikeCheck(project,mfit,N)
%MODELWORK_LOGLIKECHECK Check validity of computation and data generation.
%   MODELWORK_LOGLIKECHECK(PROJECT,MFIT) tests the computation of the log 
%   likelihood and fake data generation for project PROJECT and provided 
%   model structure MFIT.

if nargin < 3 || isempty(N); N = 1000; end

TolErr = sqrt(eps);

gendataFun = str2func([project '_gendata']);
analyticsFun = str2func([project '_analytics']);
modelFitFun = str2func([project '_modelFitBits']);
mfit.sampling = [];     % Use only maximum-likelihood
genall = gendataFun(N,mfit);

[Ntrials,Ncols] = size(genall{1});
genmat = zeros(Ntrials,Ncols,N);
ll = zeros(1,N);

options.dataid = mfit.dataid;
options.cnd = mfit.cnd;

% Regenerate true data
gentrue = gendataFun(1,mfit,'r');
[ll_true,mfit_true] = computeLL(project,gentrue,mfit,options,analyticsFun,modelFitFun);

%truedata = analyticsFun(gentrue);
%[mfit_true,infostruct] = ...
%    modelFitFun('preprocessdata',truedata{1},mfit,options,mfit.infostruct);
%mfit_true.infostruct = infostruct;
%clear functions;
%ll_true = -ModelWork_like(project,mfit_true,mfit.maptheta);

err = abs(ll_true - mfit.metrics.maploglike);
if err > TolErr
    error(['Recomputed log likelihood (LL=' num2str(ll_true,'%.4f') ') does not match stored value (LL=' num2str(mfit.metrics.maploglike,'%.4f') ').']);
end

infodata = [];

for i = 1:N;
    % genmat(:,:,i) = genall{i};
    [ll(i),mfit_curr] = computeLL(project,{genall{i}},mfit,options,analyticsFun,modelFitFun);
    
    %tmp = analyticsFun({genmat(:,:,i)});
    %gendata{i} = tmp{1};
    
    %[tmpfit,infostruct] = ...
    %    modelFitFun('preprocessdata',gendata{i},mfit,options,mfit.infostruct);
    %X = tmpfit.X;
    %clear functions;
    %ll(i) = -ModelWork_like(project,X,mfit.mp,infostruct,mfit.maptheta);
    
    [infodata,ll_est] = VestBMS_loglikeCheck(mfit_curr,mfit_true,infodata);
    
    %if ~isempty(X.bimbins{1}{2}); task = 2; else task = 3; end
    %for iNoise = 1:3
    %    p(:,:,iNoise,i) = X.bimbins{iNoise}{task};
    %    P(:,:,iNoise) = mfit_true.X.bimbins{iNoise}{task};
    %end
    %pmean = mean(p,4);
    %pmean = bsxfun(@rdivide, pmean, pmean(:,1,:)+pmean(:,2,:));
    %pmean = lambda*0.5 + (1-lambda)*pmean;
    
    %LL0 = nansum(log(pmean(:)).*P(:));
    
    fprintf('%3d# True LL: %.2f. Estimated LL: %.2f. Fake data LL: %.2f ± %.2f.\n', ...
        i, ll_true, ll_est, mean(ll(1:i)), std(ll(1:i)));
end

genmat = mean(genmat,3);

end

%--------------------------------------------------------------------------
function [ll,mfit] = computeLL(project,datamat,mfit,options,analyticsFun,modelFitFun)
%COMPUTELL Compute log likelihood of data matrix.

data = analyticsFun(datamat);
[mfit,infostruct] = ...
    modelFitFun('preprocessdata',data{1},mfit,options,mfit.infostruct);
mfit.infostruct = infostruct;
clear functions;
ll = -ModelWork_like(project,mfit,mfit.maptheta);

end
