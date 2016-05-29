function varargout = ModelWork_loadFields(varargin)
%MODELWORK_LOADFIELDS Load missing fields of model fit structure.
%
%   MBAG = MODELWORK_LOADFIELDS(MBAG,DATA) loads empty fields of model fits
%   in model bag MBAG for datasets DATA, dropped with MODELWORK_DROPFIELDS 
%   (usually the empty fields are 'X', 'mp', and 'infostruct').
%
%   MFIT = MODELWORK_LOADFIELDS(PROJECT,MFIT,DATA) loads missing fields of
%   model fit structures MFIT for project PROJECT and datasets DATA. MFIT
%   can be a single model structure or a cell array of model structures.
%
%   DATA can be left empty, in which case MODELWORK_LOADFIELDS will attempt
%   to load the data from file.
%
%   See also MODELWORK_BATCHEVAL, MODELWORK_DROPFIELDS, MODELWORK_MODELFIT.

if nargin < 1; help ModelWork_loadFields; return; end

data = [];
if isfield(varargin{1},'project')    % MBAG
    mbag = varargin{1}; 
    mfit = mbag.bag;
    project = mbag.project;
    if length(varargin) > 1; data = varargin{2}; end
elseif ischar(varargin{1})          % Cell array of MFIT structs
    mbag = [];
    project = varargin{1};
    mfit = varargin{2};
    if length(varargin) > 2; data = varargin{3}; end
else
    error('Wrong syntax for MODELWORK_LOADFIELDS.');
end

% Model-dependent functions
batchrunInitFun = str2func([project '_batchrunInit']);
modelfitFun = str2func([project '_modelFitBits']);
setupModelFun = str2func([project '_setupModel']);
    
if ~iscell(mfit); uncellflag = 1; mfit = {mfit};
else uncellflag = 0; end

for i = 1:length(mfit)
    if isempty(mfit{i}); continue; end
    
    if isempty(mfit{i}.X) || isempty(mfit{i}.mp) || isempty(mfit{i}.infostruct)
                
        if isempty(data)
            datafile = mfit{i}.datafile;
            if isempty(datafile); datafile = [project '_data.mat']; end
            data = load(datafile); 
        end
        if isfield(data,'data'); data = data.data; end 
        
        mm = mfit{i};
        optlist = ModelWork_defaults(project);
        options = parseoptions([],optlist);        
        
        fields = {'type','model','dataid','cnd'};
        for iField = 1:length(fields)
            options = setoptions(options,fields{iField},mfit{i}.(fields{iField}));
        end
        options = batchrunInitFun(data,mm.type,options);    
        [~,infostruct,options] = modelfitFun('define',data{options.dataid(1)},options.model,options);
        infostruct.cnd = options.cnd;
        [mm,infostruct] = modelfitFun('preprocessdata',data{options.dataid(1)},mm,options,infostruct);
        [mp,exitflag] = setupModelFun([],mm.maptheta,options.model,infostruct);
        mp.nData = mm.nData;
        mm.infostruct = infostruct;
        mm.mp = mp;
        mfit{i} = mm;
    end
    
end

if uncellflag; mfit = mfit{1}; end

if ~isempty(mbag)
    mbag.bag = mfit;
    varargout{1} = mbag;
else
    varargout{1} = mfit;
end
    