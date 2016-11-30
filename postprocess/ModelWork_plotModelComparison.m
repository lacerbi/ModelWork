function [models,tab,bms] = ModelWork_plotModelComparison(modelsummary,metric,bms,varargin)
% MODELWORK_PLOTMODELCOMPARISON Plot of summary results of model comparison.

if nargin < 2 || isempty(metric); metric = 'aicc'; end
if nargin < 3; bms = []; end

% Default options
options.plottype = 'square';
options.modellist = [];
options.scorebounds = 50;
options.bmsalgorithm = 'V';     % Default BMS method is variational inference
options.bmsorder = 0;
options.bmsplot = 'post';       % BMS plot: 'post', 'xp', 'pxp'
options.ssorder = 0;
options.deviance = [];          % Is the metric specified in deviance units?
options.factors = [];
options.maxmodels = 20;         % Maximum # models shown
options.nsamples = 1e6;         % Samples used in BMS
options.threshold = 5;          % Score threshold

% Parse variable arguments
options = parseoptions(options,varargin{:});

plottype = options.plottype;
modellist = options.modellist;
scorebounds = options.scorebounds;
ssorder = options.ssorder;
deviance = options.deviance;
bmsalgorithm = options.bmsalgorithm;
bmsorder = options.bmsorder;
bmsplot = options.bmsplot;
factors = options.factors;
nsamples = options.nsamples;
threshold = options.threshold;

if isempty(modellist); modellist = modelsummary.models; end

% By default perform simple model comparison
if isempty(bms); bms = 0; end

% Specify metrics that are defined in terms of deviance
deviancemets = {'aic','aicc','bic','dic','dic1','dic2','waic','waic1','waic2'};
if isempty(deviance); deviance = any(strcmpi(metric,deviancemets)); end
if deviance; devmult = -2; else devmult = 1; end

if plottype == 3; bmsorder = 1; end

fontsize = 14;
axesfontsize = 10;
MaxModels = options.maxmodels;

height = 1000;
bms_precomputed = 0;

set(gcf, 'Color', 'w');
box off;

lme = modelsummary.(metric);
models = modelsummary.models;

% Correction for Inf and NaNs in the marginal likelihood
if strcmpi(metric, 'mlike') || strcmpi(metric, 'marginallike')
    f = isnan(lme) | isinf(lme);
    if any(f(:))
        warning(['Found ' num2str(sum(f(:))) ' invalid entries (Inf''s or NaN''s) in marginal likelihood table; using BIC for those.']);
        lme(f) = -0.5*modelsummary.bic(f);
    end
end

%% Extract target models from MODELLIST

% Convert from model names to model vectors
if iscell(modellist) && ischar(modellist{1})
    modellist = getmodelvectors(modellist,modelsummary);
end

Nmodels = size(modellist,1);

pos = zeros(1,Nmodels);
for i = 1:Nmodels
    idx = find(all(bsxfun(@eq, modellist(i,:), models),2));
    if isempty(idx) || ~isscalar(idx)
        error(['Cannot find a unique match for input model #' num2str(i) ' in MODELSUMMARY.']);
    end
    pos(i) = idx;
end

% Get model names
if isfield(modelsummary, 'modelnames') && ~isempty(modelsummary.modelnames)
    modelnames = modelsummary.modelnames;
else
    for i = 1:size(models,1); modelnames{i} = ['M' num2str(i)]; end
end

modelnames = modelnames(pos);
models = models(pos,:);
lme = lme(:,pos);

%% Run model comparison

alpha0 = [];
if isscalar(bms) && isnumeric(bms) && bms == 0
    do_BMS = false;
    alpha0 = 0;
else
    do_BMS = true;  % Perform group BMS method
    if ischar(bms) && strcmpi(bms,'sparse') % Give sparse prior for alpha0
        alpha0 = alpha0 / sqrt(Nmodels);
        % bms_alpha0 = 0.285 * size(lme, 2)^(-0.3) - 0.033;
    elseif isnumeric(bms) || islogical(bms)
        alpha0 = double(bms);
    end
    if isscalar(alpha0)
        alpha0 = alpha0.*ones(1,Nmodels);
    end
    
    % Keep only models with nonzero priors
    pos = alpha0 > 0;
    lme = lme(:,pos);
    modelnames = modelnames(pos);
    models = models(pos,:);
    if ~isempty(factors); factors = factors(:,pos); end
    alpha0 = alpha0(pos);
end

[Ns,Nk] = size(lme);


if isstruct(bms) && isfield(bms, 'alpha0')
    bms_precomputed = 1;
    alpha = bms.alpha;
    exp_r = bms.exp_r;
    xp = bms.xp;
    pxp = bms.pxp;
    g = bms.g;
    bms_smpl = bms.smpl;
    models = bms.models;
    ssorder = bms.ssorder;
    modelnames = bms.modelnames;
    modelstr = modelnames;
    alpha0 = bms.alpha0;
    lme = bms.lme;
    % tab = devmult*log(g);
    tab = log(g);
    if ~strcmpi(metric, bms.metric)
        warning(['Current metric (' upper(metric) ') differs from stored BMS metric (' upper(bms.metric) ').']);
    end
elseif do_BMS
    if sum(alpha0) < 1; str = ' (sparse prior)'; else str = []; end
    display(['Performing Bayesian Model Selection (BMS) with alpha_0 = ' num2str(alpha0) str '.']);
    
    % Convert from deviance units to log likelihood units if needed
    lme = lme/devmult;
            
    switch lower(bmsalgorithm(1))
        case 'v'    % Variational inference
            [alpha,exp_r,xp,pxp,bor,g] = spm_BMS_fast(lme,[],[],[],[],alpha0.*ones(1,Nk));
            bms_smpl = [];
        case 's'    % Sampling
            error('Sampling BMS currently not supported.');
            % [exp_r,xp,g,bms_smpl] = gen_BMS_alpha(tab_adj,[],'sparse',[]);
            [exp_r,xp,g,bms_smpl] = gen_BMS(tab_adj,[],'sparse',[]);
            alpha = exp_r;
        otherwise
            error('Unknown algorithm for Bayesian Model Selection. Use either ''V''ariational inference or ''S''ampling (default ''V'')');
    end
    % tab = devmult*log(g);
    tab = log(g);
end

%% Rearrange results of model comparison

if bms_precomputed
    if ~isempty(factors); factors = factors(:,bms.modelsorder); end    
else
    if bmsorder && do_BMS
        % Order models according to posterior expectation of model probability
        [~, order] = sort(exp_r, 'descend');
    else
        % Order models according to average model score
        [~, order] = sort(nansum(tab,1), 'descend');
    end
    tab = tab(:,order);
    lme = lme(:,order);
    models = models(order,:);
    modelstr = modelnames(order);

    if do_BMS
        bms = [];
        bms.lme = lme;
        bms.alpha0 = alpha0(order);
        bms.alpha = alpha(order);
        bms.exp_r = exp_r(order);
        bms.xp = xp(order);
        bms.pxp = pxp(order);
        bms.g = g(:,order);
        bms.metric = metric;
        bms.smpl = bms_smpl;
        bms.modelnames = modelstr;
        bms.modelsorder = order; % Store reordering
        if ~isempty(bms.smpl)
            bms.smpl.r = bms.smpl.r(:,order);
            if isfield(bms.smpl,'z'); bms.smpl.z = bms.smpl.z(:,order); end
            if isfield(bms.smpl,'u'); bms.smpl.u = bms.smpl.u(:,order); end
        end
        if ~isempty(factors); factors = factors(:,order); end
        alphatot = sum(bms.alpha); % Compute SD of r
        bms.std_r = sqrt(bms.alpha.*(alphatot - bms.alpha)./(alphatot^2*(alphatot+1)));
        % tab = devmult*log(bms.g);
        tab = log(bms.g);
    else
        bms = [];
    end
    
    % tab = -tab*sign(devmult);
    
    % Reorder subjects
    if isscalar(ssorder) && ~ssorder
        ssorder = 1:Ns;
    elseif isscalar(ssorder)
        basetemp = bsxfun(@minus, tab, max(tab, [], 2));
        [bestscore,orderscore] = sort(basetemp,2,'descend');
        temp = orderscore(:, 1)*1e12;
        if size(bestscore, 2) > 1; temp = temp - bestscore(:, 2).*orderscore(:, 2)*1e6; end
        if size(bestscore, 2) > 2; temp = temp - bestscore(:, 3).*orderscore(:, 3); end
        [~,ssorder] = sort(temp, 1, 'ascend');    
    end

    if do_BMS
        if ~bms_precomputed
            bms.lme = bms.lme(ssorder,:);
            bms.g = bms.g(ssorder,:);
        end
        bms.models = models;
        bms.ssorder = ssorder;
    end
    if ~bms_precomputed; tab = tab(ssorder, :); end

end

%% Plot results of model comparison

hold on;
maxscore = scorebounds(1);

% Plot scores
switch lower(plottype(1:3))
    case {'mos','che','squ'} % Square or checkerboard table
        % Remove each individual subject's best model
        tab = bsxfun(@minus, tab, max(tab, [], 2));

        % Cap maximum difference, for visualization
        tab(tab < -maxscore) = -maxscore;

        if Nk > MaxModels
            tab = tab(:, 1:MaxModels);
            Nk = MaxModels;
        end
        
        ss = 1:Ns;
        for i = 1:Ns
            for k = 1:Nk
                patch([i i i+1 i+1], [-k -k-1 -k-1 -k], -tab(i, k), 'EdgeColor', 'none');
            end
        end

        % Best models
        for i = 1:Ns
            [besttab1, bestmodel1] = sort(tab(i, :), 'descend');
            besttab1 = besttab1 - besttab1(1);
            besttab1 = find(besttab1 > -threshold);
            for k = 1:length(besttab1)
                text(i + 0.5, -bestmodel1(besttab1(k))-0.5, num2str(k), 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', axesfontsize, 'Color', 0.8*[1 1 1]);
            end
        end

        colormap;
        ytick = round(linspace(0, maxscore, 6));
        for i = 1:length(ytick); yticklabel{i} = num2str(ytick(i)); end
        yticklabel{end} = ['> ' num2str(ytick(end))];
        h = colorbar('YTick', ytick, 'YTickLabel', yticklabel);
        set(h, 'YDir', 'reverse');
        h2 = get(h, 'Ylabel');
        set(h2, 'Interpreter', 'Tex');
        axis([1 (length(ss)+1), (-Nk-1) -1, -height height]);

        % Axes
        set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial');

        for i = 1:length(ss); xticklabel{i} = num2str(ssorder(i)); end
        set(gca, 'Xtick', 1.5:1:length(ss)+1, 'XtickLabel', xticklabel);

        for i = 1:Nk
            yticklabel{Nk-i+1} = modelstr{i};
        end
        set(gca, 'Ytick', (-Nk:-1)-0.5, 'YtickLabel', yticklabel);

        xlabel('Subject number', 'FontName', 'Arial', 'FontSize', fontsize);

    case 'bar' % Bar graph
        tab = bsxfun(@minus, tab, tab(:, 1));
        tabbase = tab;
        
        % Standard group-average model comparison
        if ~do_BMS
            tabse = stderr(tab);
            tab = nanmean(tab);
                    
            % Compute p-value
            [~,pval] = ttest(tabbase);
            % pval = pval*size(tab, 2)

            % x-axis limit
            if length(scorebounds) > 1; xrlimit = scorebounds(2); else xrlimit = scorebounds(1) + 50; end

        % Hierarchical Bayesian model comparison    
        else            
            tab = bms.exp_r;
            tabse = bms.std_r;
            pval = 1 - bms.pxp;
            xrlimit = 1; % x-axis limit
        end
        
        if Nk > MaxModels
            tab = tab(:, 1:MaxModels);
            Nk = MaxModels;
        end

        % Remove ticks
        % patch([0.1 20 20 0.1], [-size(tab, 2) -size(tab, 2) -0.6 -0.6], [1 1 1], 'EdgeColor', 'none');

        h = barh(-1:-1:-Nk, tab');
        colmap = colormap;
        % colmap = flipud(colormap);
        for k = 1:Nk
            hb(k) = barh(-k, tab(k), 0.8);
            colindex = max(1,1 + floor(min([tab(k)/maxscore, 1])*(size(colmap, 1)-1)));
            set(hb(k), 'FaceColor', colmap(colindex, :));
        end
        
        % Plot errors
        for i = 1:Ns
            for k = 1:Nk
                if isempty(hb(i,k)); continue; end
                x = hb(i,k).XData;
                y = hb(i,k).YData;
                clipping = 'on';
                
                seleft = [[1 1]*tab(i,k) + [0 1]*tabse(i,k), [1 1]*x];
                plot(seleft(1:2), seleft(3:4), 'k', 'LineWidth', 1, 'Clipping', clipping);
                if tab(i,k) > maxscore*0.1 && tab(i,k) < maxscore*0.7; col2 = [0 0 0]; else col2 = 0.8*[1 1 1]; end
                seright = [[1 1]*tab(i,k) + [-1 0]*tabse(i,k), [1 1]*x];
                plot(seright(1:2), seright(3:4), 'Color', col2, 'LineWidth', 1, 'Clipping', clipping);
            end
        end    

        % Plot pvalues
        for k = 1:Nk
            if pval(k) > 0.05 || isnan(pval(k)); continue; end
            if pval(k) < 0.001; ptext = '***'; elseif pval(k) < 0.01; ptext = '**'; else ptext = '*'; end 
            text(xrlimit*0.9, -k, ptext, 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', fontsize);        
%            text(tab(i) + 30, -k, ptext, 'HorizontalAlignment', 'Left', 'FontName', 'Arial', 'FontSize', fontsize);        
        end

        % symb = {'', 'd', 'o'};
        % col = [0 0 0; 0.9 1 1; 0.5 0.8 0.8];
        % Plot nodes
        %if ngroups == 2
        %    for i = 2:3
        %        for j = 1:size(dic, 2)
        %            hnodes(i) = plot(dic(i, j), -j -0.25*(i-2.5), [symb{i} 'k-'], 'MarkerFaceColor', col(i, :), 'LineWidth', 0.5);            
        %        end
        %    end
        %end

        % Best models
        for i = 1:size(tab, 1)
            [bestscore1, bestmodel1] = sort(tab(i, :), 'descend');
            bestscore1 = bestscore1 - bestscore1(1);
            bestscore1 = find(bestscore1 > -threshold);
        end

        colormap;
        % colormap(flipud(colormap));
        axis([0 xrlimit, (-size(tab,2)-0.5) -0.5]);

        % Axes
        set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial', 'Xgrid', 'on');        

        for i = 1:size(tab,2)
            yticklabel{size(tab,2)-i+1} = modelstr{i};
        end
        set(gca, 'Ytick', (-size(tab,2):-1), 'YtickLabel', yticklabel);

        if plottype == 2
            xlabel(['{\Delta}' upper(metric)], 'FontName', 'Arial', 'FontSize', fontsize, 'Interpreter', 'Tex');
        elseif plottype == 3
            xlabel('Posterior EMF', 'FontName', 'Arial', 'FontSize', fontsize, 'Interpreter', 'Tex');            
        end
        
    case 'fac' % Factor graph
        lme = bms.lme;        
        %tab = bsxfun(@minus, tab, tab(:, 1));
        %tabbase = tab;
        
        if isempty(factors)
            error('FACTORS not specified; cannot compute factor comparison.');
        end
        Nfactors = size(factors,1);
        
        % Normalize factors
        factors_norm = bsxfun(@rdivide, factors, sum(factors,1));
        
        horizontalplot = 0;
        
        switch lower(bmsalgorithm)
            case 'v'
                
                % Add alpha levels of each model to each feature
                alpha_fac = sum(bsxfun(@times, bms.alpha, factors),2)';
                alpha0_fac = sum(bsxfun(@times, bms.alpha0, factors),2)';
                        
                % Compute posterior summed over factors
                g_fac = zeros(Ns, Nfactors);
                lme_fac = zeros(Ns, Nfactors);
                for f = 1:Nfactors
                    g_fac(:,f) = sum(bsxfun(@times, bms.g, factors_norm(f,:)), 2);
                    idx = logical(factors_norm(f,:));
                    maxlme = max(lme(:,idx),[],2);
                    lme_fac(:,f) = log(sum(bsxfun(@times,factors_norm(f,idx)/sum(factors_norm(f,idx)),exp(bsxfun(@minus,lme(:,idx),maxlme))),2)) + maxlme;
                end
                
                % Compute protected exceedance probability for each factor
                [xp,bor,pxp] = spm_dirichlet_exceedance_fast(alpha_fac, nsamples, lme_fac, g_fac, alpha0_fac);

%                 tmp = bsxfun(@minus, lme_fac, max(lme_fac,[],2));
%                 bsxfun(@rdivide, exp(tmp), sum(exp(tmp),2))
%                 xp
%                 bor
%                 pxp
                
                % 'Bayesian' p-value
                pval = 1 - pxp;
                xrlimit = 1; % x-axis limit
                
                switch lower(bmsplot)
                    case 'xp'
                        tab = xp;
                        tabse = zeros(size(tab));
                    case 'pxp'
                        tab = pxp;
                        tabse = zeros(size(tab));
                    case 'post'
                        % Compute posterior expected feature frequency and SD
                        tab = alpha_fac./sum(alpha_fac);
                        tabse = sqrt(alpha_fac.*(sum(alpha_fac)-alpha_fac)/(sum(alpha_fac)^2*(sum(alpha_fac)+1)));
                end
                
                nonfeatures = (alpha_fac == 0);
                tab(nonfeatures) = []; tabse(nonfeatures) = [];
                
            case 's' % Sampling                
                error('Sampling method currently not supported.');
                
                exp_r = zeros(size(bms.smpl.r,1), length(modelgroups));
                for iModel = 1:size(models, 1)
                    model = models(iModel,:);
                    flag = zeros(1, length(modelgroups));
                    for iGroup = 1:length(modelgroups)
                        if any(all(bsxfun(@eq, model, modelgroups{iGroup}), 2))
                            flag(iGroup) = 1;
                        end
                    end
                    % Sum contributions of current model
                    exp_r = exp_r + bsxfun(@times,bms.smpl.r(:, iModel),flag)/sum(flag);
                end
                
                % Compute posterior expected feature frequency and SD
                tab = mean(exp_r,1)./sum(mean(exp_r,1));
                tabse = std(exp_r,[],1);
        
                nonfeatures = (tab == 0);
                tab(nonfeatures) = []; tabse(nonfeatures) = [];
                
                % Compute excess probability (i.e., the probability of being the 
                % most likely component)
                [~,iMax] = max(exp_r,[],2);
                for iGroup = 1:length(modelgroups)
                    xp(iGroup) = sum(iMax == iGroup)/size(exp_r,1);
                end
                pval = 1 - xp;
                xrlimit = 1; % x-axis limit
                
            otherwise
                error('Unknown algorithm for Bayesian Model Selection. Use either ''V''ariational inference or ''S''ampling (default ''V'')');
        end
        
        if horizontalplot
            h = barh(-1:-1:-size(tab, 2), tab');
        else
            h = bar(1:size(tab, 2), tab', 'EdgeColor', 'none');
        end
        colmap = colormap;
        for k = 1:size(tab, 2)
            if horizontalplot
                hb(k) = barh(-k, tab(k), 0.8);
            else
                hb(k) = bar(k, tab(k), 0.8);
            end
            % colindex = 1 + floor(min([tab(i)/maxscore, 1])*(size(colmap, 1)-1));        
            % set(hb(i), 'FaceColor', colmap(colindex, :));
            set(hb(k), 'FaceColor', 0.8*[1 1 1], 'EdgeColor', 'none');
        end
        
        % Plot errors
        for i = 1:size(tab,1)
            for k = 1:size(tab, 2)
                if isempty(hb(i, k)); continue; end
                x = hb(i,k).XData;
                y = hb(i,k).YData;
                if horizontalplot
                    plot([1 1]'*tab(i, k) + [0 1]'*tabse(i, k), [1 1]'*y, 'k', 'LineWidth', 0.5, 'Clipping', 'off');
                else
                    plot([1 1]'*k, [1 1]'*tab(i,k) + [0 1]'*tabse(i, k), 'k', 'LineWidth', 0.5, 'Clipping', 'off');                        
                end
                   
                if tab(i,k) > maxscore*0.1 && tab(i, k) < maxscore*0.7; col2 = [0 0 0]; else col2 = 0.8*[1 1 1]; end
                if horizontalplot
                    plot([1 1]'*tab(i, k) + [-1 0]'*tabse(i, k), [1 1]'*y, 'Color', col2, 'LineWidth', 0.5, 'Clipping', 'off');
                else
                    plot([1 1]'*k, [1 1]'*tab(i, k) + [-1 0]'*tabse(i, k), 'Color', col2, 'LineWidth', 0.5, 'Clipping', 'on');                    
                end
            end
        end    

        % Plot pvalues
        for k = 1:size(tab, 2)
            if pval(k) > 0.05 || isnan(pval(k)); continue; end
            if pval(k) < 0.001; ptext = '***'; elseif pval(k) < 0.01; ptext = '**'; else ptext = '*'; end 
            if horizontalplot
                text(xrlimit*0.9, -k, ptext, 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', fontsize);
            else
                text(k, xrlimit*0.95, ptext, 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', fontsize);                
            end
%            text(tab(i) + 30, -i, ptext, 'HorizontalAlignment', 'Left', 'FontName', 'Arial', 'FontSize', fontsize);        
        end

        % colormap;
        % colormap(flipud(colormap));
        if horizontalplot; axis([0 xrlimit, (-size(tab,2)-0.5) -0.5]); else axis([0.5 size(tab,2)+0.5, 0 xrlimit]); end

        % Axes
        set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial','TickDir','out'); %,'TickLength',2*get(gca,'TickLength'));
        % box on;

        if isfield(modelsummary,'groupnames') && ~isempty(modelsummary.groupnames)
            str = [];
            modelsummary.groupnames
            for k = 1:length(nonfeatures)
                if ~nonfeatures(k); str{end+1} = modelsummary.groupnames{k}; end
            end
            modelstr = str;
        else
            for iGroup = 1:Nfactors
                modelstr{iGroup} = ['F#' num2str(iGroup)];
            end
        end
        
        if horizontalplot
            for k = 1:size(tab,2)
                yticklabel{size(tab,2)-k+1} = modelstr{k};
            end
            set(gca, 'Xtick', 0:0.5:1, 'Ytick', (-size(tab,2):-1), 'YtickLabel', yticklabel);
            xlabel('Posterior Model Feature Frequency', 'FontName', 'Arial', 'FontSize', fontsize, 'Interpreter', 'Tex');
        else
            for k = 1:size(tab,2)
                xticklabel{k} = modelstr{k};
            end
            set(gca, 'Ytick', 0:0.5:1, 'Xtick', 1:size(tab,2), 'XtickLabel', xticklabel);
            xticklabel_rotate([],45);
            ylabel('Posterior Model Feature Frequency', 'FontName', 'Arial', 'FontSize', fontsize);
        end
end

end

%--------------------------------------------------------------------------
function models = getmodelvectors(modelnames,modelsummary)
%GETMODELVECTORS Extract model vector from modelsummary given model names.

M = numel(modelnames);
modelidx = zeros(1,M);
for i = 1:M
    idx = find(strcmp(modelnames{i},modelsummary.modelnames));
    if isempty(idx) || ~isscalar(idx)
        error(['Cannot find a unique match for model ' modelnames{i} ' in model table.']);
    end
    modelidx(i) = idx;
end
models = modelsummary.models(modelidx,:);

end

%--------------------------------------------------------------------------
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