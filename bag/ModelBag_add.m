function mbag = ModelBag_add(mbag,m,method,project,data)
% MODELBAG_ADD add a model fit to a model bag.
%   MBAG = MODELBAG_ADD(MBAG,M) adds M to the model bag MBAG. M can be
%   a single fit or another modelbag struct. If MBAG is the empty object,
%   create a new model bag. M can be also a vector of fits.
%
%   MBAG = MODELBAG_ADD(MBAG,M,METHOD) specifies the method to use
%   with duplicate models:
%      'e'  error               give an error message
%      'm'  merge (default)     merge models and chains
%      'r'  replace             replace original model fit in MBAG with M
%      's'  skip                keep original model fit in MBAG
%
%   MBAG = MODELBAG_ADD(MBAG,M,METHOD,PROJECT) passes project name 
%   PROJECT to a new model bag.
%
%   MBAG = MODELBAG_ADD(MBAG,M,METHOD,PROJECT,DATA) passes dataset DATA
%   to a new model bag.

TABLESIZE = 9973;                   % Table size for hash indexing

if nargin < 1
    help ModelBag_add;
    return;
end
if nargin < 3 || isempty(method); method = 'merge'; end
if nargin < 4; project = []; end
if nargin < 5; data = []; end

if length(m) > 1
    for i = 1:length(m)
        mbag = ModelBag_add(mbag,m{i},method,project,data);
    end
    return;
end

if isempty(project) && isfield(mbag,'project')
    project = mbag.project;
end

% Initialize model bag if empty
if isempty(mbag)
    mbag.bag = [];
    mbag.project = project;
    mbag.model = [];
    mbag.index = [];
    mbag.dataid = [];
    mbag.cnd = [];
    mbag.ver = '2015a';
    mbag.tablesize = TABLESIZE;
    mbag.table{mbag.tablesize} = [];
    if isempty(mbag.project); warning('Project name PROJECT not specified.'); end
else
    mbag.project = project;    
end

if isempty(m); return; end

% Check if it is a model bag
if isstruct(m) && isfield(m,'bag')
    for i = 1:length(m.bag)
        mbag = ModelBag_add(mbag,m.bag{i},method,[],data);
    end
    return;
end

% Find model in model bag, otherwise append model
model = m.model;
% Retrocompatibility
if isfield(m,'dataid'); dataid = m.dataid; else dataid = m.nid; end
cnd = m.cnd;

% Find index of given model
index = findmodelhash(mbag,dataid,model,cnd);

% If the model already exists, behavior depends on the flag
if ~isempty(index)
    switch lower(method(1))
        case 'e'; error('Duplicate model fits in MBAG and M.');
        case 's';
        case 'r';  mbag.bag{index} = m;
        case 'm'; mbag.bag{index} = ModelBag_merge(project,mbag.bag{index},m,data);
        case 'x'
            error('EXTENDCHAINS method not supported yet.');
            
            if mbag.bag{index}.nchains ~= m.nchains
                error('Number of chains is different.');
            else
                error('Chain extension not implemented yet.');

                newn = size(m.smpl, 1) + size(mbag.bag{index}.smpl, 1);
                newsmpl = zeros(newn, size(m.smpl, 2)); 
                newloglikes = zeros(newn, size(m.loglikes, 2));
                oldnpc = mbag.bag{index}.nsamplesperchain;
                addnpc = m.nsamplesperchain;

                for ii = 1:m.nchains
                    newsmpl((1:oldnpc) + (ii-1)*(oldnpc+addnpc), :) = mbag.bag{index}.smpl((1:oldnpc) + (ii-1)*oldnpc, :);
                    newsmpl(oldnpc + (1:addnpc) + (ii-1)*(oldnpc+addnpc), :) = m.smpl((1:addnpc) + (ii-1)*addnpc, :);
                    newloglikes((1:oldnpc) + (ii-1)*(oldnpc+addnpc), :) = mbag.bag{index}.loglikes((1:oldnpc) + (ii-1)*oldnpc, :);
                    newloglikes(oldnpc + (1:addnpc) + (ii-1)*(oldnpc+addnpc), :) = m.loglikes((1:addnpc) + (ii-1)*addnpc, :);
                end                        
            end
            mbag.bag{index}.smpl = newsmpl;
            mbag.bag{index}.loglikes = newloglikes;
            mbag.bag{index}.nsamplesperchain = oldnpc + addnpc;
        otherwise; error('Unknown METHOD for dealing with duplicate models.');
    end

else

    % Model does not exist, append at the end
    index = length(mbag.bag) + 1;        
    mbag.bag{index} = m;

    % Hash indexing
    h = modelhash(dataid,model,cnd,mbag.tablesize);
    mbag.table{h} = [mbag.table{h}, index];                        

end
