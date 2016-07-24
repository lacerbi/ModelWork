function ModelWork_batchEval(project,data,jobs,varargin)
%MODELWORK_BATCHEVAL Run a batch model evaluation.
%
%   MODELWORK_BATCHEVAL(PROJECT,DATA,JOBS) runs batch model evaluation for 
%   project PROJECT, dataset DATA and job array JOBS. 
%   DATA can either be a data struct or a string file name. 
%   JOBS can be either a cell array of jobs structs or a string filename for a
%   job text file created with MODELWORK_MAKEJOBLIST.
%
%   MODELWORK_BATCHEVAL(...,NAME,VALUE) sets the field NAME of the options 
%   struct to a given VALUE.
%
%   See also MODELWORK_MAKEJOBLIST, MODELWORK_MODELFIT.

% By Luigi Acerbi <luigi.acerbi@gmail.com>
% Last update: 22/08/2015

starttime = tic; % Start measuring time

if nargin < 1; help ModelWork_batchEval; return; end
if isempty(project); error('You need to specify the project name PROJECT.'); end

% Read job list from file
if ischar(jobs)
    fin = fopen(jobs,'r');
    jobs = [];
    joblist = textscan(fin,'%d %s %s %s %s');
    for i = 1:length(joblist{1})
        jobs{i}.type = joblist{1}(i);
        jobs{i}.model = unpackuint(joblist{2}{i});
        jobs{i}.dataid = unpackuint(joblist{3}{i});
        jobs{i}.cnd = fullstr2num(joblist{4}{i});
        jobs{i}.replica = fullstr2num(joblist{5}{i});
    end
    fclose(fin);
end

% Check JOBS format as cell array
if isstruct(jobs) && length(jobs) > 1; 
    error('Struct arrays not supported for JOBS. Use cell arrays instead.');
end
if ~iscell(jobs); jobs = {jobs}; end
if ~all(isfield(jobs{1},{'type','model','dataid','cnd','replica'}))
    error('JOBS is not a valid job file or a cell array of job structs.');
end

% Load default options (general and project-specific)
optlist = ModelWork_defaults(project);
options = parseoptions(varargin,optlist);
    
% Run each job sequentially
for iProc = 1:length(jobs)    
    batchEval_one(project,data,jobs{iProc},starttime,options);
end

end

%--------------------------------------------------------------------------
function batchEval_one(project,data,job,starttime,options)
% BATCHEVAL_ONE Run a single job.

jobdone = 0;    % Job is not done yet
fitstep = 0;    % Step of fitting (0: none; 1: optimization; 2:sample)
exitflag = 0;   % No problems so far

% Retrieve process id
if isempty(options.procid)
    procid = 0;
else
    procid = options.procid;
end

% Basic functions
defaultsFun = str2func([project '_defaults']);
batchrunInitFun = str2func([project '_batchrunInit']);

% Load data structure
[data,datafile] = ModelWork_loaddatafile(project,data,job.type);
if ~isempty(datafile); setoptions(options,'datafile',datafile); end

MAXLENGTH = 32;             % Add trailing zeros to data/model vectors
job.dataid = [job.dataid,zeros(1,MAXLENGTH-length(job.dataid))];
job.model = [job.model,zeros(1,MAXLENGTH-length(job.model))];

% Initialize project options according to job TYPE
options = setoptions(options,'type',job.type,1);
options = setoptions(options,'model',job.model,1);
options = setoptions(options,'dataid',job.dataid,1);
options = setoptions(options,'cnd',job.cnd,1);
options = setoptions(options,'replica',job.replica,1);
options = setoptions(options,'procid',procid,1);
options = batchrunInitFun(data,job.type,options);

% Get formatted strings for model and data
[modelstring,dataidstring] = defaultsFun('strings',options.model,options.dataid);
options = setoptions(options,'modelstring',modelstring,1);
options = setoptions(options,'dataidstring',dataidstring,1);

% Sampling or optimizing?
% options.samplingflag = options.nsamples > 0;

% Load file and path information
fileinfo = ModelWork_fileinfo(options,job);

% Create sub-directories if needed
if ~exist(fileinfo.basedir1,'dir'); mkdir(fileinfo.basedir1); end
if ~exist([fileinfo.basedir1 filesep() fileinfo.basedir2],'dir'); mkdir([fileinfo.basedir1 filesep() fileinfo.basedir2]); end

% Print full logs on file (print to video on my laptop)
fileinfo.mylaptop = options.mylaptop;
options = setoptions(options,'stdout',strncmp(fileinfo.hostname,fileinfo.mylaptop,length(fileinfo.mylaptop)),1);
fileinfo.outfile = []; fileinfo.successfile = [];
fileinfo.failfile = ['fail-' num2str(procid) '.log'];
if ~options.stdout
    fileinfo.outfile = ['diary-' num2str(procid) '.log'];
    fileinfo.successfile = ['success-' num2str(procid) '.log'];
end

options = setoptions(options,'outfile',fileinfo.outfile,1);
options = setoptions(options,'filename',fileinfo.filename,1);
options = setoptions(options,'fullfilename',fileinfo.fullfilename,1);

if exist(fileinfo.fullfilename,'file')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RECOVERY JOB
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load run variables from file
    prevrun = load(fileinfo.fullfilename);
    % options = prevrun.options;
    % job = prevrun.job;
    mbag = prevrun.mbag;
    jobdone = prevrun.jobdone;
    fitstep = prevrun.fitstep;
        
    if jobdone && ~options.continueflag
        error(['File ' fileinfo.fullfilename ' already exists and job done.']);
    end
    
    % Check that loaded file matches with current job
    assert(all(prevrun.job.type == job.type) ...
        && all(prevrun.job.model == job.model) ...
        && all(prevrun.job.dataid == job.dataid) ...
        && all(prevrun.job.cnd == job.cnd) ...
        && (prevrun.job.replica == job.replica), ...
        'Mismatch between job in recovered run and current job.');
    
    % Recompute the 'scratch' field
    %options2 = batchrunInitFun(data,job.type,options);
    %if isfield(options2, 'scratch'); options.scratch = options2.scratch; end
    
    % Number of samples may change in a re-iteration
    %if jobdone && options.continueflag
    %    options.niter = options2.niter;
    %    options.nsamples = options2.nsamples;
    %end
    %clear options2;
    
    options.continueflag = 1;
    
    writelog(fileinfo.outfile,'recover',procid,job,options);

    clear prevrun;
    
% Otherwise, start new run
else
    mbag = [];              % Empty model bag
    options = setoptions(options,'starttime',starttime);        
    writelog(fileinfo.outfile,'init',procid,job,options,options.continueflag);
    options.continueflag = 0;       % Start from scratch   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BATCH RUN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear functions; % Clear persistent variables

try
    
    if ~isempty(options.mbag)
        mfit = ModelBag_get(options.mbag,job.dataid,job.model,job.cnd);
        % Get rid of MBAG in OPTIONS, but keep track if it was used
        if ~isempty(mfit)
            options.mbag = 'used';
        else
            options.mbag = [];
        end
    elseif options.continueflag && ~isempty(mbag)
        mfit = ModelBag_get(mbag,job.dataid,job.model,job.cnd);
    else
        mfit = [];
    end

    if isempty(options.optfevals)
        error('OPTIONS.optfevals is empty. Set OPTIONS.optfevals = 0 if you want to skip optimization.');
    end
    
    % 1: Multi-start optimization
    if options.optfevals > 0 && fitstep == 0
        [mfit,exitflag] = ModelWork_modelFit(project,data,mfit,'optimize',options);
        fitstep = 1;
        if ~isempty(mfit) && exitflag == 0        
            mbag = ModelBag_add(mbag,mfit,'merge',project);
        end
        save(fileinfo.fullfilename);
    end
        
    % 2: Sample
    if options.nsamples > 0
        [mfit,exitflag] = ModelWork_modelFit(project,data,mfit,'sample',options);
        if ~isempty(mfit) && exitflag == 0        
            mbag = ModelBag_add(mbag,mfit,'replace',project);
        end
    end
    
    % 3: Recompute model statistics
    if options.optfevals == 0 && options.nsamples == 0
        mfit = ModelWork_modelStats(project,mfit,options.hessianflag,[],1,options);
        if ~isempty(mfit)
            mbag = ModelBag_add(mbag,mfit,'replace',project);
        end
    end
    
    %if ~isempty(mfit) && exitflag == 0        
    %    mbag = ModelBag_add(mbag,mfit,'merge',project);
    %end

catch fitexception
    jobpath = path;
    jobcd = cd;
    clear data; % Do not save data

    save(fileinfo.fullfilename);

    writelog(fileinfo.outfile,'fail_short',procid,job,options);    
    writelog(fileinfo.failfile,'fail_long',procid,job,options,fitexception);
    
    dataidstring = ['S' num2str(job.dataid(1)) '-' packuint(job.dataid(2:end))];
    error(['Error in model ' packuint(job.model) ', ' dataidstring ', cnd ' numarray2str(job.cnd) ', replica ' num2str(job.replica) '.']);
    exitflag = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOB INTERRUPTED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exitflag > 0
    writelog(fileinfo.outfile,'interrupt',procid,job,options);    
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOB DONE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jobdone = 1;
clear ans data mfit batchrunInitFun defaultsFun;
if isfield(options,'scratch'); options = rmfield(options, 'scratch'); end;

% MATLAB may fail while saving and leave a corrupt file -- at least avoid
% that it ruins both the old and the new file
if exist(fileinfo.fullfilename,'file')
    [pathstr,name,ext] = fileparts(fileinfo.fullfilename);
    tempfile = fullfile(pathstr,[name '.old']);
    % tempfile = ['.' filesep pathstr filesep name filesep 'old']
    movefile(fileinfo.fullfilename,tempfile);
    for iAttempt = 1:10
        temp = [];
        try
            save(fileinfo.fullfilename);
            temp = load(fileinfo.fullfilename); % See if it worked
            if ~isempty(temp); break; end
        catch
            % Did not work
            writelog(fileinfo.successfile,'savefail',procid,job,options,iAttempt);
        end
    end
    if isempty(temp); error('Cannot save current model, aborting.'); end
    % Remove temporary old file
    if exist(tempfile,'file'); delete(tempfile); end
else
    save(fileinfo.fullfilename);
end

writelog(fileinfo.outfile,'done',procid,job,options);
writelog(fileinfo.successfile,'success',procid,job,options);

end   % MODELWORK_BATCHEVAL


%--------------------------------------------------------------------------
function writelog(outfile,state,iProc,job,options,varargin)
%WRITELOG Write entry in log file.

if isempty(outfile); return; end
modelstring = options.modelstring;
dataidstring = options.dataidstring;

fout = fopen(outfile, 'a');
switch lower(state)
    case 'init'    
        fprintf(fout, '\n\n-----------------------------------------------------------\n');
        fprintf(fout, '%s: Start model fitting for thread #%d.\n', datestr(now), iProc);
        fprintf(fout, '%s: Model %s, D%s, cnd %s, replica %d.\n', datestr(now), modelstring, dataidstring, numarray2str(job.cnd), job.replica);
        continueflag = varargin{1};
        if continueflag
            fprintf(fout, '%s: Could not find checkpoint file. Starting from scratch.\n', datestr(now));
        end
        
    case 'recover'
        fprintf(fout, '\n\n-----------------------------------------------------------\n');
        fprintf(fout, '%s: Recover model fitting for thread #%d.\n', datestr(now), iProc);
        fprintf(fout, '%s: Model %s, D%s, cnd %s, replica %d.\n', datestr(now), modelstring, dataidstring, numarray2str(job.cnd), job.replica);                
        
    case 'fail_short'
        fprintf(fout, 'ERROR in Model %s, D%s, cnd %s, replica %d.\n', modelstring, dataidstring, numarray2str(job.cnd), job.replica);
        
    case 'fail_long'
        fitexception = varargin{1};
        fprintf(fout, '%s: Job ID #%d failed.\n', datestr(now), iProc);
        fprintf(fout, 'ERROR in Model %s, D%s, cnd %s, replica %d.\n', modelstring, dataidstring, numarray2str(job.cnd), job.replica);
        fprintf(fout, 'File ''%s'', function ''%s'', line %d.\n', fitexception.stack(1).file, fitexception.stack(1).name, fitexception.stack(1).line);
        fprintf(fout, '%s\n%s\n', fitexception.identifier, fitexception.message);        
        
    case 'interrupt'
        fprintf(fout, '%s: Job interrupted (may be able to recover) for thread #%d.\n', datestr(now), iProc);
        fprintf(fout, '-----------------------------------------------------------\n');

    case 'done'
        fprintf(fout, '%s: End model fitting for thread #%d.\n', datestr(now), iProc);
        fprintf(fout, '-----------------------------------------------------------\n');
                
    case 'success'
        fprintf(fout, '%s: Job ID #%d was successfully completed.\n', datestr(now), iProc);
        fprintf(fout, '%s: Model %s, D%s, cnd %s, replica %d.\n', datestr(now), modelstring, dataidstring, numarray2str(job.cnd), job.replica);                

    case 'savefail'
        iAttempt = varargin{1};
        fprintf(fout, '%s: Error while saving job ID #%d at completion (attempt #%d).\n', datestr(now), iProc, iAttempt);
        fprintf(fout, '%s: Model %s, D%s, cnd %s, replica %d.\n', datestr(now), modelstring, dataidstring, numarray2str(job.cnd), job.replica);                
        fprintf(fout, 'Trying again...\n');
        
end     
fclose(fout);
end

%--------------------------------------------------------------------------
function num = fullstr2num(str)
% FULLSTR2NUM Convert string to number, everything that is not numbers or 
% points is converted to whitespace.

if isnumeric(str); num = str; return; end
str(str < 48 | str > 57) = ' ';
num = str2num(str);
end