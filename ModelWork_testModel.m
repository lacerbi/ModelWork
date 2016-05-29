
function mfit = ModelWork_testModel(project,data,type,modeln,datan,varargin)

clear functions;

if ischar(data); data = load(data); end 

optlist = ModelWork_defaults(project);
options = parseoptions(varargin,optlist);

% Initialize project options according to job TYPE
defaultsFun = str2func([project '_defaults']);
batchrunInitFun = str2func([project '_batchrunInit']);
options = setoptions(options,'type',type,1);
options = setoptions(options,'chain',1,1);
options = setoptions(options,'processid',1,1);
[options,models,dataids,cnds] = batchrunInitFun(data,type,options);

models
dataids

% Selected model, dataset (subject), and condition
model = models(modeln,:);
dataid = dataids(datan,:);
cnd = cnds{datan}{1};

options = setoptions(options,'model',model,1);
options = setoptions(options,'dataid',dataid,1);
options = setoptions(options,'cnd',cnd,1);

% Get formatted strings for model and data
[modelstring,dataidstring] = defaultsFun('strings',options.model,options.dataid);
options = setoptions(options,'modelstring',modelstring,1);
options = setoptions(options,'dataidstring',dataidstring,1);

% Sampling or optimizing?
options.nsamples = 0;
options.niter = 1;
options.samplingflag = options.nsamples > 0;
options = setoptions(options,'nstarts',0,1);

% Compute one step of minimization per epoch
[mfit,exitflag] = ModelWork_modelFit(project,data,options);

% Recover model fit fields
mfit = ModelWork_loadFields(project,mfit,data);

end