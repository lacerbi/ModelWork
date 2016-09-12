function options = parseoptions(options,varargin)
%PARSEOPTIONS Parse options either as struct or variable arguments in name/value format.
%   OPTIONS = PARSEOPTIONS(OPTIONS,'PROPERTY1',VALUE1,'PROPERTY2',VALUE2,...)
%   sets the fields propertyX in default OPTIONS structure to valueX. Input 
%   field names are not case sensitive. 
%
%   OPTIONS = PARSEOPTIONS(OPTIONS,NEWOPTS) assigns fields in struct NEWOPTS
%   to OPTIONS. Matching fields are not case sensitive for OPTIONS.
%
%   OPTIONS = PARSEOPTIONS(OPTIONS,NEWOPTS,'PROPERTY1',VALUE1,'PROPERTY2',VALUE2,...) 
%   first assigns values from struct NEWOPTS, and then name/value pairs.
%

%   Author: Luigi Acerbi
%   Email:  luigi.acerbi@gmail.com
%   Date:   Sep/08/2016

if nargin < 1; help parseoptions; return; end

if isempty(options)
    error('parseOptions:emptyOptions','Default OPTIONS struct should be nonempty.');
end

if isempty(varargin)
    return;
end

deff = fields(options)';

if isstruct(varargin{1})                    % Input as NEWOPTS struct
    newopts = varargin{1};
    for f = fields(newopts)'
        idx = find(strcmpi(f{:}, deff),1);
        if isempty(idx)
            error('parseOptions:unknownProperty', ...
                ['Unknown property ''' f{:} ''' in NEWOPTS.']);
        else
            options.(deff{idx}) = newopts.(f{:});
        end
    end
    varargin(1) = [];
end

if ~isempty(varargin)
    if ischar(varargin{1})                      % Input in name/value format
        % check for correct number of inputs
        if mod(numel(varargin),2) == 1
            error('parseOptions:wrongInputFormat', ...
                'Name and value input arguments must come in pairs.');
        end

        % parse arguments
        for i = 1:2:numel(varargin)
            if ischar(varargin{i})
                idx = find(strcmpi(varargin{i}, deff),1);
                if isempty(idx)
                    error('parseOptions:unknownProperty', ...
                        ['Unknown property name ''' varargin{i} '''.']);
                else
                    options.(deff{idx}) = varargin{i+1};
                end            
            else
                error('parseOptions:wrongInputFormat', ...
                    'Name and value input arguments must come in pairs.');
            end
        end
    else
        error('parseOptions:wrongInputFormat', ...
                'Input should come as a NEWOPTS struct and/or with name/value pairs.');
    end
end

end