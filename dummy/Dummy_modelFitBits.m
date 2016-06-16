

function varargout = Dummy_modelFitBits(command, varargin)

switch lower(command)

    % function [datalikeFun, infostruct] = ModelDependentDefine(dataone)
    % MODELDEPENDENTDEFINE Define model-dependent variables.
    case 'define'
        dataone = varargin{1};
        model = varargin{2};
        options = varargin{3};

        % Define data function
        datalikeFun = @Dummy_loglike;        
        varargout{1} = datalikeFun;
        
        if ~isempty(options)             
            cnd = options.cnd;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Specific data for model-parameters structure
            %infostruct.experimentName = options.experimentName;
            %if isempty(infostruct.experimentName)
            %    error('Expriment name not specified.');
            %end
            
            % 'Temperature' during slice sampling (multiplies window step)
            infostruct.samplingtemperature = options.samplingtemperature;

            % Convert model number to model vector
            % if isscalar(model); model = CueBMS_getmodelnumber(model); end
            
            varargout{2} = infostruct;
            varargout{3} = options;
        end
                
    % function modelstruct = ModelDependentInit(dataone, modelstruct, options)
    %  Model-dependent preprocessing of data structure
    case 'preprocessdata'
        
        modelstruct = varargin{2};
        options = varargin{3};
        infostruct = varargin{4};
        
        modelstruct.X = varargin{1}.X;                    
        modelstruct.nData = size(modelstruct.X,1);
        
    %----------------------------------------------------------------------
    % DATASET FLAGS: Datasets keep/remove trials
    % Flag 1-7: Available
    %----------------------------------------------------------------------
        
        options.dataid = [options.dataid 0];
                
        varargout{1} = modelstruct;
        varargout{2} = infostruct;
        
    % function cost = ModelDependentCost(dataone, model)
    % MODELCOST Model-dependent cost.
    case 'cost'
        
        model = varargin{1};
        cnd = varargin{2};
               
        % Compute approximate model computational cost
        % priornmix = [1 1 2 2];
        % lapsecost = [1 1.25 1.5 1.25 1.5 1.25 1.25];
        cost = ones(1, length(cnd));
        for j = 1:length(cnd)
           % ...        
        end   
        
        varargout{1} = sum(cost);
        
end

end