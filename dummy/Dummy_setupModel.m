% DUMMY_SETUPMODEL initialize and update model parameter structure.
%
% MP = CUEBMS_SETUPMODEL([], THETA, MODEL, INFOSTRUCT) returns a
% model-parameter structure MP initialized with parameter vector THETA 
% according to the model specified in model vector MODEL. INFOSTRUCT is a
% structure whose fields contain additional initialization information.
% THETA can be empty, in which case the parameter structure is only 
% partially initialized.
%
% MP = CUEBMS_SETUPMODEL(MP, THETA) updates existing model-parameter 
% structure MP with global parameter values THETA.
%
% [MP, OUTFLAG] = CUEBMS_SETUPMODEL(...) returns the flag OUTFLAG 
% which takes the value 1 if there were errors, 0 otherwise.

% The MODEL array contains, in order:
%--------------------------------------------------------------------------
% MODEL(1) Lapse:
% 1 No lapse; 2 Lapse (1 param)

function [mp,exitflag] = Dummy_setupModel(mp,theta,model,infostruct)

if ~exist('theta', 'var'); theta = []; end
if ~exist('model', 'var'); model = []; end
if ~exist('infostruct', 'var'); infostruct = []; end

% Initialize model
if isempty(mp)
    [mp, exitflag] = initModel(model, infostruct);
    if ~isempty(theta) && exitflag == 0
        [mp, exitflag] = updateModel(mp,theta);
    end
else
    [mp, exitflag] = updateModel(mp,theta);
end

end

% INITMODEL Initialize model.
function [mp, outflag] = initModel(model, infostruct)
    mp = [];
    outflag = 0;
    
    % Take information from infostruct
    cnd = infostruct.cnd;
                    
    % Number of function calls in the optimization
    mp.funcalls = 0;
    
    mp.model = model;
    mp.cnd = cnd;
    mp.ncnd = length(cnd);
        
    % For advanced sampling parameters - is kappa going to start from zero?
    mp.kappaminusoneflag = 1;

    % Fixed parameters
    % if ~isfield(infostruct, 'fixedtheta'); infostruct.fixedtheta = []; end
    % mp.fixedtheta = infostruct.fixedtheta;
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setup model parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Bounds are maximal lower/upper bounds, reasonable lower/upper bounds,
    % and reasonable width of the posterior (for slice sampling)

    pbounds = [];
    params = [];
    
    % Number of parameter families
    mp.nparams = zeros(1, 2);
    for i = 1:size(mp.nparams, 2); pbounds{i} = []; params{i} = []; end
    
    % Parameter shared across conditions
    if ~isfield(infostruct, 'sharedparams'); infostruct.sharedparams = []; end
    if isempty(infostruct.sharedparams);
       % By default, all parameters are shared
       infostruct.sharedparams = ones(1, size(mp.nparams, 2));
    end
    mp.sharedparams = infostruct.sharedparams;
    %if ~mp.sharedparams(5)
    %    error('Sensory mapping parameters need to be shared between sessions.');
    %end
    
    mp.theta = [];
    
    mubound = [-3,3,-1,1,1];
    sigmabound = log([0.05 10, 0.5 5, exp(1)]);
    lambdabound = log([1e-6 0.5, 1e-6 0.1, exp(3)]);

    % Visual noise variance
    % Set default values
    for icnd = 1:mp.ncnd
        mp.nparams(1) = mp.nparams(1) + 2;
        params{1} = {'mu','sigma'};
        pbounds{1} = [mubound; sigmabound];        
    end
    
    % Lapse model
    for icnd = 1:mp.ncnd
        mp.fulltheta{icnd}.lambda = 0; 
    end    
    switch model(1)
        case 1 % No lapse
        case 2 % Lapse (1-param)
            mp.nparams(2) = 1;
            pbounds{2} = lambdabound;
            params{2} = {'lambda'};
        otherwise; error('Unsupported lapse model.');
    end 
    
    %----------------------------------------------------------------------
    % Finalize parameters

    % Unwrap parameters
    pboundslist = []; paramslist = [];
    for i = 1:size(mp.nparams, 2)
        if ~isempty(pbounds{i})
            pboundslist = [pboundslist; pbounds{i}];
            for j = 1:length(params{i}); paramslist{end+1} = params{i}{j}; end
        end
    end

    mp.nparams = mp.nparams + mp.nparams.*(1-mp.sharedparams)*(mp.ncnd-1);
    
    for icnd = 2:mp.ncnd
        for ii = 1:size(mp.sharedparams, 2)
            if ~mp.sharedparams(ii)
                pboundslist = [pboundslist; pbounds{ii}];
                for j = 1:length(params{ii}); paramslist{end+1} = [params{ii}{j} '#' num2str(icnd)]; end
            end
        end
    end
        
    % Total number of parameters
    mp.ntotparams = sum(mp.nparams, 2);
    
    % for icnd = 1:mp.cnd; mp.fulltheta{icnd} = zeros(1, mp.nparams); end
    
    % Parameters extreme lower bounds and upper bounds
    mp.bounds.LB = zeros(1, mp.ntotparams); mp.bounds.UB = zeros(1, mp.ntotparams); 
    
    % Parameters reasonable lower bounds and upper bounds
    mp.bounds.RLB = zeros(1, mp.ntotparams); mp.bounds.RUB = zeros(1, mp.ntotparams);
        
    for k = 1:mp.ntotparams
        mp.bounds.LB(k) = pboundslist(k, 1);
        mp.bounds.UB(k) = pboundslist(k, 2);
        mp.bounds.RLB(k) = pboundslist(k, 3);
        mp.bounds.RUB(k) = pboundslist(k, 4);
        mp.bounds.SCALE(k) = pboundslist(k, 5);                
    end    
    
    mp.paramstree = params;
    mp.params = paramslist;
    
end

% UPDATEMODEL Adjust global model parameters based on current
% model and given parameter vector. Returns a flag OUTFLAG which is
% one if everything went okay or zero otherwise.
function [mp,exitflag] = updateModel(mp,theta)

    mp.theta = theta;
    model = mp.model;

    % Set XGRID values for different types of computation (precision)
    if ~isfield(mp, 'computation') || isempty(mp.computation)
        mp.computation = 1;  % Standard level of precision
    end
        
    exitflag = 0;
    actparam = 1;    

   % Information for log priors 

    for icnd = 1:mp.ncnd
        thiscnd = mp.cnd(icnd);

        % Mu and sigma
        updateparams(icnd, 1, {'id','exp'});
                
        % Lapse model
        switch model(1)
            case 1 % No lapse
            case 2 % Lapse (1-param)
                updateparams(icnd, 2, {'exp'});
        end
                
        % mp.fulltheta{1}
    end
    
   
    return;
    
    % UPDATEPARAMS Generic parameter update.
    function updateparams(icnd, iParam, oplist)        
        if mp.sharedparams(iParam) && icnd > 1
            for kk = 1:length(mp.paramstree{iParam})
                mp.fulltheta{icnd}.(mp.paramstree{iParam}{kk}) = mp.fulltheta{1}.(mp.paramstree{iParam}{kk});
            end
        else
            for kk = 1:length(mp.paramstree{iParam})
                th = theta(actparam - 1 + kk);
                switch oplist{kk}
                    case 'id' % Identity
                    case 'exp' % Exponential (for log representation)
                        th = exp(th);
                    case 'inv' % Inverse representation
                        th = 1/th - 1;
                end
                mp.fulltheta{icnd}.(mp.paramstree{iParam}{kk}) = th;
            end
            actparam = actparam + length(mp.paramstree{iParam});            
        end
    end
    
end