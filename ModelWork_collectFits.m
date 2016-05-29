function [mbag,modelsummary] = ModelWork_collectFits(project,name,data)
% MODELWORK_COLLECTFITS Collect model fits in a single structure.
%
%   MBAG = MODELWORK_COLLECTFITS(PROJECT) collect all model fits in this
%   directory, assigning them to project PROJECT.
%
%   MBAG = MODELWORK_COLLECTFITS(PROJECT,NAME) also recursively load model
%   fits from directories that match pattern NAME. Use '*' to recursively
%   get model fits from all subdirectories of the current directory.
%
%   [MBAG,MODELSUMMARY] = MODELWORK_COLLECTFITS(...) computes a model
%   summary struct.
%
% Author: Luigi Acerbi 
% Email: luigi.acerbi@gmail.com
% Date: 10/Oct/2015

if nargin < 2 || isempty(name); name = []; end
if nargin < 3 || isempty(data); data = load([project '_data.mat']); end

mbag = [];

% Loop over mat files in the directory
if name(1) == '*'
    for f = dir('*.mat')'
        fitdata = load(f.name);
        if isfield(fitdata,'mbag')
            fprintf('%s\n', f.name);
            mbag = ModelBag_add(mbag,fitdata.mbag,[],project);
        end
    end
else
    fprintf('\\. (skipping files in current directory)\n');
end

% Recursively check sub-directories
if ~isempty(name)
    for d = dir(name)'
        if ~d.isdir || strcmp(d.name,'.') || strcmp(d.name,'..')
            continue;
        end
        cd(d.name);
        fprintf('\\%s\n', d.name);
        mbag = ModelBag_add(mbag,ModelWork_collectFits(project,'*'),[],project);
        cd('..');        
    end
end

% Compute model summary if requested
if nargout > 1
    modelsummary = ModelWork_summary(mbag);
end

end