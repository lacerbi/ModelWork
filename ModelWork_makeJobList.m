function ModelWork_makeJobList(project,data,type,replicas,nprocs,varargin)
%MODELWORK_MAKEJOBLIST Create text files with job list.
%
%   MODELWORK_MAKEJOBLIST(PROJECT,DATA,TYPE,REPLICAS,NPROCS) creates job
%   list for project PROJECT, dataset DATA, type TYPE and replica REPLICAS.
%   DATA can either be a data struct or a string file name. 
%   REPLICAS is a numeric array of replicas.
%   NPROCS is a numeric scalar.
%
%   MODELWORK_MAKEJOBLIST(...,NAME,VALUE) sets the field NAME of the 
%   options struct to a given VALUE.
%
%   See also MODELWORK_BATCHEVAL.

if nargin < 1; help ModelWork_makeJobList; return; end

if nargin < 5; nprocs = []; end

verbose = 1;

batchrunInitFun = str2func([project '_batchrunInit']);
costFun = str2func([project '_modelFitBits']);

if verbose; fprintf('Loading data structure...\n'); end

% Load data structure
data = ModelWork_loaddatafile(project,data,type);

if verbose; fprintf('Loading default options...\n'); end

%% Default options
optlist = ModelWork_defaults(project);
options = parseoptions(varargin,optlist);

if verbose; fprintf('Loading model-dependent instructions...\n'); end

% Model-dependent instructions
[options,models,dataids,cnd] = batchrunInitFun(data,type,options);

if verbose; fprintf('Creating job list...\n'); end

% Create job list
jobs = CreateJobsList(costFun,models,dataids,cnd,replicas);

if verbose; fprintf('Removing completed jobs from list...\n'); end

% Remove completed jobs from list
defaultsFun = str2func([project '_defaults']);
removejobs = false(1,numel(jobs));
for iJob = 1:length(jobs)
    thisjob = jobs{iJob};    
    [options.modelstring,options.dataidstring] = ...
        defaultsFun('strings',thisjob.model,thisjob.dataid); 
    fileinfo = ModelWork_fileinfo(options,thisjob);    
    if exist(fileinfo.fullfilename,'file')
        temp = load(fileinfo.fullfilename,'jobdone');
        if temp.jobdone; removejobs(iJob) = true; end
    end
end
jobs(removejobs) = [];

% If not specified, assign one processor per job
if isempty(nprocs) || ~isfinite(nprocs); nprocs = numel(jobs); end

if verbose; fprintf('Assigning jobs to processors...'); end

for iProc = 1:nprocs
    if verbose; fprintf('%d..', iProc); end
    
    jobproc = SplitJobs(iProc,nprocs,jobs);
    jobfilename = [num2str(iProc) '.job'];
    fout = fopen(jobfilename,'w');
    for iJob = 1:length(jobproc)
        thisjob = jobproc{iJob};
        model = packuint(thisjob.model);
        dataid = packuint(thisjob.dataid);
        cnd = thisjob.cnd;
        replica = thisjob.replica;
        fprintf(fout,'%d %s %s %s %s\n',type,model,dataid,numarray2str(cnd),numarray2str(replica));
    end
    fclose(fout);    
end
if verbose; fprintf('\n'); end

end

%--------------------------------------------------------------------------
function jobs = SplitJobs(id,nprocs,jobs)
% SPLITJOBS Split a jobs list among processes 

costJobs = NaN(1,length(jobs));
for i = 1:length(jobs)
    if isfield(jobs{i},'cost'); costJobs(i) = jobs{i}.cost; end        
end

% Compute all lists
loclist = [];
costProcess = zeros(1, nprocs); % Current cost for each process
for i = 1:nprocs; loclist{i} = []; end
while 1
   if all(isnan(costJobs)); break; end
   [~, processindex] = min(costProcess); % Pick the process with less cost
   [cost, listindex] = max(costJobs); % Pick the job with highest cost
   loclist{processindex}{end+1} =  jobs{listindex};
   costJobs(listindex) = NaN;
   jobs{listindex} = [];       
   costProcess(processindex) = costProcess(processindex) + cost;
end

jobs = loclist{id};

% Display locally chosen models
% length(jobs)
% for i = 1:length(jobs); display(jobs{i}); end
    
end

%--------------------------------------------------------------------------
function jobs = CreateJobsList(costFun,models,dataids,cnd,replicas)
% CREATEJOBSLIST Create a list of independent jobs to run. 
% Each job is a combination of model, subject, conditions and replica.

jobs = [];

for iModel = 1:size(models, 1);
    job.model = models(iModel,:);
    for iData = 1:size(dataids, 1);            
        job.dataid = dataids(iData,:);
        cnd1 = cnd{job.dataid(1)};    % Subject-dependent cell array of cnds
        if ~iscell(cnd1); cnd1 = {cnd1}; end
        for iCnd = 1:length(cnd1)
            job.cnd = cnd1{iCnd};
            job.cost = costFun('cost', job.model, job.cnd);
            for iChain = 1:length(replicas)
                job.replica = replicas(iChain);
                jobs{end+1} = job;
            end
        end
    end
end        
    
end