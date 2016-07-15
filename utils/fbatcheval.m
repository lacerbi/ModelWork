function Y = fbatcheval(fun,X,fout,filename,cptime,startmsg)
%FBATCHEVAL Batch evaluation with checkpoints of expensive function.
%   Y = FBATCHEVAL(FUN,X) evaluates function FUN on a number of samples X,
%   where X is a matrix with one sample per row. FUN accepts as input an
%   array of the same size as a row of X and returns either a scalar or a
%   row array. Y is a matrix such that Y(i,:) is FUN(X(i,:)).
%
%   Y = FBATCHEVAL(FUN,X,FOUT) writes log on file FOUT (default FOUT = 1,
%   standard output).
%
%   Y = FBATCHEVAL(FUN,X,FOUT,FILENAME) uses file FILENAME to periodically
%   save checkpoints.
%
%   Y = FBATCHEVAL(FUN,X,FOUT,FILENAME,CPTIME) checkpoints every CPTIME
%   seconds (default CPTIME=1000).
%
%   Y = FBATCHEVAL(FUN,X,FOUT,FILENAME,CPTIME,STARTMSG) writes message
%   STARTMSG at beginning of execution.
%

%  Date: 16/07/2016
%  Author: Luigi Acerbi
%  Email: luigi.acerbi@gmail.com

if nargin < 3 || isempty(fout); fout = 1; end   % Standard output
if nargin < 4; filename = []; end
if nargin < 5 || isempty(cptime); cptime = 1000; end
if nargin < 6 || isempty(startmsg); startmsg = 'Computing function... Sample '; end

% Evaluate first sample, get output size
Ytemp = fun(X(1,:));
k = numel(Ytemp);

n = size(X,1);
Y = NaN(n,k);
Y(1,:) = Ytemp(:);

if ~isempty(filename) && exist(filename,'file')
    %% Load interrupted execution from file
    
    % Store current options
    fout_new = fout; filename_new = filename;
    cptime_new = cptime; startmsg_new = startmsg;
    
    load(filename,'-mat');  % Load run from recovery file
    i0 = iSamp;                     % Start from recovered iteration
    if fout > 0; fprintf('Loading sampling from file ''%s''.\n', filename); end

    fout = fout_new; filename = filename_new;
    cptime = cptime_new; startmsg = startmsg_new;
    clear fout_new filename_new cptime_new startmsg_new;
else
    i0 = 2;
end

if fout > 0; fprintf(fout, startmsg); end

lastsave = tic;

for iSamp = i0:n
    
    % Save current progress to file
    if ~isempty(filename) && (toc(lastsave) > cptime || iSamp == n)
        if fout > 0; fprintf('\nSaving temp file ''%s''...\n', filename); end
        
        if ~exist(filename,'file')
            save(filename);
        else

            % File already exists, move it to temporary copy
            [pathstr,name,ext] = fileparts(filename);
            tempfile = fullfile(pathstr,[name '.old_']);
            movefile(filename,tempfile);
            for iAttempt = 1:3
                temp = [];
                try
                    save(filename);
                    temp = load(filename,'-mat'); % See if it worked
                    if ~isempty(temp); break; end
                catch
                    % Did not work
                end
            end
            if isempty(temp)
                movefile(tempfile,filename);
                error(['Cannot save file ''' filename '''.']);
            end

            % Remove temporary old file
            if exist(tempfile,'file'); delete(tempfile); end
        end
        
        lastsave = tic;
    end
    
    if mod(iSamp,100) == 0 && fout > 0; fprintf(fout,'%d..', iSamp); end
    Ytemp = fun(X(iSamp,:));
    Y(iSamp,:) = Ytemp(:);
end

if fout > 0; fprintf(fout, '\n'); end

end
