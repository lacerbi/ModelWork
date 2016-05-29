function fileinfo = ModelWork_fileinfo(options,job)
%MODELWORK_FILEINFO Return common file and storage path information.
%
%   FILEINFO = MODELWORK_FILEINFO(OPTIONS,JOB) returns a struct FILEINFO
%   with common file and storage path information for options struct
%   OPTIONS and job struct JOB.

% Sampling or optimizing?
if ~isfield(options,'samplingflag') || isempty(options.samplingflag)
    options.samplingflag = options.nsamples > 0;
end
% if options.samplingflag; suffix = 'smp'; else suffix = 'opt'; end
suffix = 'fit';

fileinfo.basedir1 = [options.jobname '@C' numarray2str(job.cnd,[],'','','')];
fileinfo.basedir2 = ['M' options.modelstring];
fileinfo.filename = [options.dataidstring '@' suffix num2str(job.replica) '.mat'];
fileinfo.fullfilename = [fileinfo.basedir1 filesep() fileinfo.basedir2 filesep() fileinfo.filename];
if ~isfield(options,'hostname') || isempty(options.hostname)
    [~,fileinfo.hostname] = system('hostname');
else
    fileinfo.hostname = options.hostname;
end