function varargout = ModelWork_dropFields(mbag,options)
%MODELWORK_DROPFIELDS Drop unnecessary fields of model fit structure.
%
%   MBAG = MODELWORK_DROPFIELDS(MBAG) drop unnecessary fields of model fits
%   in model bag MBAG (usually dropped fields are 'X', 'mp', and 
%   'infostruct').
%
%   MFIT = MODELWORK_LOADFIELDS(MFIT) drop unnecessary fields of model fit
%   structures MFIT. MFIT can be a single model structure or a cell array 
%   of model structures.
%
%   MFIT = MODELWORK_LOADFIELDS(...,OPTIONS) uses additional information
%   in options struct OPTIONS (currently unused).
%
%   Dropped fields can be recovered via MODELWORK_LOADFIELDS.
%
%   See also MODELWORK_LOADFIELDS.

if nargin < 1; help ModelWork_dropFields; return; end
if nargin < 2; options = []; end

if isfield(mbag,'bag')     % MBAG
    mfit = mbag.bag;
else                    % Cell array of MFIT structures
    mfit = mbag;
    mbag = [];
end
    
if ~iscell(mfit); uncellflag = 1; mfit = {mfit};
else uncellflag = 0; end

for i = 1:numel(mfit)
    if isempty(mfit{i}); continue; end
    
    % These fields can be recovered with MODELWORK_LOADFIELDS
    mfit{i}.X = [];
    mfit{i}.infostruct = [];
    mfit{i}.mp = [];

    % Remove temporary scratch field
    if isfield(mfit{i},'scratch')
        mfit{i} = rmfield(mfit{i}, 'scratch');
    end

    % Remove function handles (keep for retrocompatibility)
    if ~isempty(options) && isfield(options,'removefunhandles') && ~isempty(options.removefunhandles) 
        fun_handles = options.removefunhandles;
        for j = 1:numel(mfit.mp.fulltheta)
            for k = 1:numel(fun_handles)
                if isfield(mfit.mp.fulltheta{j}, fun_handles{k})
                    mfit.mp.fulltheta{j}.(fun_handles{k}) = [];
                end
            end        
        end
    end
end

if uncellflag; mfit = mfit{1}; end

if ~isempty(mbag)
    mbag.bag = mfit;
    varargout{1} = mbag;
else
    varargout{1} = mfit;
end
    