function optlist = ModelWork_defaults(project)
%MODELWORK_DEFAULTS Return list of default options.
%
%   OPTLIST = MODELWORK_DEFAULTS() returns the standard default options
%   list. OPTLIST is a struct array (see PARSEOPTIONS for description).
%
%   OPTLIST = MODELWORK_DEFAULTS(PROJECT) also loads the default options
%   list for a given project.
%
%   See also PARSEOPTIONS.

if nargin < 1; project = []; end

%% INPUT/OUTPUT
optlist(1) = struct('name','stdout','type','flag','default',[]);
optlist(end+1) = struct('name','hostname','type','string','default',[]);
optlist(end+1) = struct('name','mylaptop','type','string','default','Dyson');
optlist(end+1) = struct('name', 'display', 'type', 'string', 'default', 'all');
optlist(end+1) = struct('name', 'outfile', 'type', 'string', 'default', []);

%% RECOVERY AND TROUBLESHOOTING
optlist(end+1) = struct('name', 'continueflag', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'mbag', 'type', 'string', 'default', []);
optlist(end+1) = struct('name', 'appendchains', 'type', 'flag', 'default', 0);
optlist(end+1) = struct('name', 'skipcompletedjobs', 'type', 'matrix', 'default', 1);
optlist(end+1) = struct('name', 'debug', 'type', 'flag', 'default', 0);
optlist(end+1) = struct('name', 'clearfunctions', 'type', 'matrix', 'default', 1);

%% CURRENT RUN DEFINITION
optlist(end+1) = struct('name', 'experimentName', 'type', 'string', 'default', []);
optlist(end+1) = struct('name', 'datafile', 'type', 'string', 'default', []);
optlist(end+1) = struct('name', 'type', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'model', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'procid', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'dataid', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'chain', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'cnd', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'modelstring', 'type', 'string', 'default', []);
optlist(end+1) = struct('name', 'dataidstring', 'type', 'string', 'default', []);
% optlist(end+1) = struct('name', 'mfit', 'type', 'matrix', 'default', []);

%% GENERAL MODEL FITTING
optlist(end+1) = struct('name', 'firstSetting', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'speedtest', 'type', 'matrix', 'default', Inf);
optlist(end+1) = struct('name', 'nstarts', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'seed', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'startx', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'loadstartx', 'type', 'flag', 'default', 0);
optlist(end+1) = struct('name', 'savetime', 'type', 'matrix', 'default', 3600);
optlist(end+1) = struct('name', 'binnedloglik', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'logpriorflag', 'type', 'flag', 'default', 0);

%% OPTIMIZATION
optlist(end+1) = struct('name', 'optimizationmethod', 'type', 'string', 'default', '');
optlist(end+1) = struct('name', 'optfevals', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'optrefinefevals', 'type', 'matrix', 'default', []);
%optlist(end+1) = struct('name', 'niter', 'type', 'matrix', 'default', []);

%% SAMPLING
optlist(end+1) = struct('name', 'samplingmethod', 'type', 'string', 'default', 'eiss');
optlist(end+1) = struct('name', 'nsamples', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'maxstoredsamples', 'type', 'matrix', 'default', 5000);
optlist(end+1) = struct('name', 'burnin', 'type', 'matrix', 'default', []);
optlist(end+1) = struct('name', 'thin', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'samplingflag', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'samplingtemperature', 'type', 'matrix', 'default', 1); % Obsolete       
% optlist(end+1) = struct('name', 'kebabsample', 'type', 'flag', 'default', 0);

%% POSTPROCESSING
optlist(end+1) = struct('name', 'recomputesamplingmetrics', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'computemarginallike', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'hessianflag', 'type', 'matrix', 'default', 0);
optlist(end+1) = struct('name', 'removefunhandles', 'type', 'string', 'default', []);

if ~isempty(project)
    defaultsFun = str2func([project '_defaults']);
    optlist = [optlist, defaultsFun('options')];
end
