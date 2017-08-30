function fakedata = ModelWork_generateFakeData(project,mbag,nid,models,cnd,N)
% MODELWORK_GENERATEFAKEDATA Generates fake datasets.
%
% Common usage:
% fakedata = ModelWork_generateFakeData(mbag,modelsummary.nid,modelsummary.models,modelsummary.cnd)

warning('You may want to use MODELWORK_GENDATA.');

setupModelFun = str2func([project '_setupModel']);
gendataFun = str2func([project '_gendata']);
analyticsFun = str2func([project '_analytics']);

fakedata = [];

% Loop over models
for im = 1:size(models, 1)
    display(['Generating datasets for model #' num2str(im) '...'])
    
    if size(cnd, 1) == size(models, 1)
        modelcnd = cnd(im, :);
    else
        modelcnd = cnd;
    end
    model = models(im, :);
    mfit = ModelBag_get(mbag, nid, model, {modelcnd});
        
    thetastar = zeros(length(mfit), size(mfit{1}.maptheta, 2));

    % Take MAP solution
    for i = 1:length(mfit); thetastar(i, :) = mfit{i}.maptheta; end
    
    % Compute mean and covariance matrix of fitted parameters
    mustar = mean(thetastar);
    sigmastar = cov(thetastar);
    
    for k = 1:N
        id = mod(k-1, length(mfit)) + 1;
        X = mfit{id}.X;
        mp = mfit{id}.mp;
        
        while 1
            try

                % Generate valid theta
                while 1
                    thetarnd = mvnrnd(mustar, sigmastar);
                    if all(thetarnd >= mp.bounds.LB & thetarnd <= mp.bounds.UB); break; end            
                end

                % Update parameter structure
                [mp, outflag] = setupModelFun(mp, thetarnd, model, mfit{id}.infostruct);

                % Generate single dataset
                X = gendataFun(1, X, mp);
                tempdata = analyticsFun(X,[],[],1);
                % tempdata = tempdata{1};
                tempdata{1}.id = mfit{id}.nid;
                tempdata{1}.truemodel = model;
                tempdata{1}.truetheta = thetarnd;
                fakedata = [fakedata, tempdata];
                break;
            catch mexc
                mexc
                some_error = 1;
                warning('MATLAB:ModelWork:generateFakeData', ...
                    ['Error while generating dataset #' num2str(k) ' of model #' num2str(im) '. Retrying...']);
            end
        end
    end
             
end

% gendata = CueBMS_addUnimodalFit(gendata, mfituni);

% Update datasets' id (needs to be unique), keep true one
for i = 1:length(fakedata)
    fakedata{i}.id
    fakedata{i}.trueid = fakedata{i}.id;
    fakedata{i}.id = i;        
end

end