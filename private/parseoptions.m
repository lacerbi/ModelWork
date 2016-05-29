function options = parseoptions(optstring,optlist,oldoptions)
%PARSEOPTIONS parses a list of strings and returns an options struct.
%
%   OPTIONS = PARSEOPTIONS(OPTSTRING,OPTLIST) parses the list of strings
%   and values OPTSTRING (typically a VARARGIN) using the options in 
%   OPTLIST (see below) and returns a struct OPTIONS.
%
%   OPTLIST is a struct array with the following fields:
%      'name'       name of the option
%      'type'       'flag','matrix','string' (flags can take values 0 or 1)
%      'default'    default value
%
%   OPTIONS = PARSEOPTIONS(OPTSTRING,OPTLIST,OLDOPTIONS) uses the options
%   in OLDOPTIONS instead of the defaults.

if nargin < 3; oldoptions = []; end

if isempty(oldoptions)
    options = struct();
else
    options = oldoptions;
end

% Check parsed struct for options
for i = 1:length(optstring)
    if i > length(optstring); break; end
    if isempty(optstring{i}); continue; end
    
    match = 0;
    % Check for parameter in argument list (start from end of the list)
    for j = length(optlist):-1:1
        opt = optlist(j);
        if strcmpi(optstring{i}, opt.name)
          switch lower(opt.type)
              case 'flag'
                  options = setoptions(options,opt.name,1);
              case 'matrix'
                  temp = cell2mat(optstring(i+1));
                  if ischar(temp); temp = str2double(temp); end
                  options = setoptions(options,opt.name,temp);
                  optstring{i+1} = [];
              case 'string'
                  options = setoptions(options,opt.name,optstring{i+1});
                  optstring{i+1} = [];                        
          end
          match = 1;
          break;
        end
    end
    
   if ~match
        error(['Error in parsing options: unknown option ''', optstring{i}, '''.']);
   end     
end

% Assign default values to empty fields
for j = 1:length(optlist)
    name = optlist(j).name;
    flag = 1;
    if ~isempty(oldoptions) && isfield(oldoptions,name)
        defvalue = oldoptions.(name);
        if isfield(oldoptions,'userset_') && ...
                any(strcmp(name,oldoptions.userset_)); flag = 0; end
    else
        defvalue = optlist(j).default;
    end
    options = setoptions(options,name,defvalue,flag);
end