function varargout = Dummy_defaults(command,varargin)
%PROJECT_DEFAULTS Get default variables for a given project.
%
%   OPTLIST = PROJECT_DEFAULTS('options') returns a project-dependent 
%   options list. OPTLIST is a struct array with the following fields:
%      'name'       name of the option
%      'type'       'flag','matrix','string' (flags can take values 0 or 1)
%      'default'    default value
%
%   [MODELSTRING,DATAIDSTRING] = PROJECT_DEFAULTS('strings',MODEL,DATAID) 
%   returns the project-dependent model string and data id string for
%   a given MODEL and DATAID.
%
%   See also MODELWORK_DEFAULTS, PARSEOPTIONS.

switch lower(command)
    case 'options'
        optlist(1) = struct('name', 'dummyflagoption', 'type', 'flag', 'default', false);
        optlist(end+1) = struct('name', 'dummymatrixoption', 'type', 'matrix', 'default', 1);
        optlist(end+1) = struct('name', 'dummystringoption', 'type', 'string', 'default', 'dummiest');
        
        varargout{1} = optlist;
        
    case 'strings'        
        model = varargin{1};
        dataid = varargin{2};
        
        modelstring = packuint(model);
        suffix = packuint(dataid(2:end));
        if isempty(suffix)
            dataidstring = ['S' num2str(dataid(1))];                
        else
            dataidstring = ['S' num2str(dataid(1)) '-' packuint(dataid(2:end))];
        end
        
        varargout{1} = modelstring;
        varargout{2} = dataidstring;          
end