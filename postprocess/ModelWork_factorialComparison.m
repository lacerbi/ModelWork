function [bms,fac] = ModelWork_factorialComparison(modelsummary,metric,modelnames,modelbase,prior,factors,factornames,hg,masks,factorfixed)
%MODELWORK_FACTORIALCOMPARISON Performa and plot factorial model comparison.

if nargin < 8; hg = []; end
if nargin < 9; masks = []; end
if nargin < 10 || isempty(factorfixed); factorfixed = 0; end

if ~iscell(factors); factors = {factors}; end
if ~iscell(hg); hg = {hg}; end
Nfgroups = numel(factors);

if isempty(masks); masks{1} = ones(1,size(factors{1},2)); end
if ~iscell(masks); masks = {masks}; end

Nmasks = numel(masks);
if isscalar(factorfixed); factorfixed = factorfixed * ones(Nmasks,Nfgroups); end

for m = 1:Nmasks
    hg_panel = hg{m};
    
    prior_one = prior(min(m,end),:);
    
    if isempty(hg_panel)
        grid = [ones(3,Nfgroups), 2:Nfgroups+1];
        hg = plotify(grid,'gutter',[0.05,0.1],'margins',[0.05 0.05,0.1 0.05]);
    end

    if hg_panel(1) ~= 0
        axes(hg_panel(1));
        plotfirst = 1;
    else
        plotfirst = 0;
    end
    bms = plotcomparison(modelsummary,metric,modelnames,modelbase,prior_one,1,[],[],plotfirst);
    if plotfirst
        h = get(gca,'xlabel');
        set(h,'FontSize',12);    
    end

    for k = 1:Nfgroups
        axes(hg_panel(k+1));
        labelxticks = (m == Nmasks);
        [~,tab] = plotcomparison(modelsummary,metric,modelnames,modelbase,prior_one,bms, ...
            factornames{k},bsxfun(@times,factors{k},masks{m}),1,labelxticks,factorfixed(m,k));
        
        fac{m,k}.tab = tab;
        fac{m,k}.factors = factors{k};
        fac{m,k}.factornames = factornames{k};

        set(gca,'FontSize',10,'Ytick',[0 0.5 1]);
        if k == 1;
            ylabel('$\tilde{p}$','Interpreter','LaTeX','FontSize',12)
        else
            ylabel('');
        end
        if k > 1; set(gca,'YTickLabel',''); end
    end
end

end

%--------------------------------------------------------------------------
function [bms,tab] = plotcomparison(modelsummary,metric,modelnames,bestmodel,prior,BMS,factornames,factors,plotflag,labelxticks,factorfixed)

if nargin < 5 || isempty(prior); prior = 1; end
if nargin < 6 || isempty(BMS); BMS = 0; end
if nargin < 7; factornames = []; end
if nargin < 8; factors = []; end
if nargin < 9; plotflag = 1; end
if nargin < 10; labelxticks = []; end
if nargin < 11; factorfixed = 0; end

if (isnumeric(BMS) || islogical(BMS)) && BMS == 0
    ModelPlot_compare(modelsummary,lower(metric),bestmodel,'mean',[],modelnames);
    bms = [];
else        % Bayesian Model Selection for group studies
    if ~isempty(factors)    % Remove models not present in factors
        prior = prior .* (sum(factors,1) > 0);
    end
    M = sum(prior > 0);
    fprintf('Number of effective models: %d.\n', M);
    alpha0 = prior./sum(prior) * M; %/ sqrt(M);
    if (isnumeric(BMS) || islogical(BMS)) && BMS == 1
        plottype = 'bars';
        BMSplot = 'post';
    else
        plottype = 'factors';
        BMSplot = 'pxp';
        alpha0 = alpha0;
    end
    
    [models,tab,bms] = ModelWork_plotModelComparison(modelsummary,lower(metric),alpha0, ...
        'PlotType',plottype,'ModelList',modelnames,'SSorder',1,'BMSorder',0, ...
        'Factors',factors,'FactorNames',factornames,'BMSplot',BMSplot,'Plot',plotflag, ...
        'LabelXTicks', labelxticks, 'FactorFixed', factorfixed);    
end

end