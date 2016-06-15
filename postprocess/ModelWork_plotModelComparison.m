function [models,tab,bms] = ModelWork_plotModelComparison(modelsummary,score,plottype,modellist,scorebounds,bms_alpha0,bmsorder,ssorder)
% MODELWORK_PLOTMODELCOMPARISON Plot of summary results of model comparison.

if ~exist('plottype', 'var') || isempty(plottype); plottype = 1; end
if ~exist('modellist', 'var') || isempty(modellist); modellist = modelsummary.models; end
if ~exist('scorebounds', 'var') || isempty(scorebounds); scorebounds = 100; end
if ~exist('bms_alpha0', 'var') || isempty(bms_alpha0); bms_alpha0 = 0; end
if ~exist('bmsorder', 'var') || isempty(bmsorder); bmsorder = 0; end
if ~exist('ssorder', 'var') || isempty(ssorder); ssorder = 0; end

if plottype == 3; bmsorder = 1; end

bmsalgorithm = 2; % Sampling

fontsize = 16;
axesfontsize = 12;

height = 1000;
bms_precomputed = 0;

MAXMODELS = 20; % Maximum number of models shown on the table

% opengl software;
set(gcf, 'Color', 'w');
box off;

tab = modelsummary.(score);
models = modelsummary.models;

% Correction for Inf and NaNs in the marginal likelihood
if strcmpi(score, 'mlike') || strcmpi(score, 'marginallike')
    f = isnan(tab) | isinf(tab);
    if any(f(:))
        warning(['Found ' num2str(sum(f(:))) ' invalid entries (Inf''s or NaN''s) in marginal likelihood table; using BIC for those.']);
        tab(f) = -0.5*modelsummary.bic(f);
    end
    
end

mf = false(size(models, 1), 1);

% If MODELLIST is a cell array, flatten it
if iscell(modellist)
    modelgroups = modellist;
    modellist = modelgroups{1};
    for iGroup = 2:length(modelgroups)
        for iModel = 1:size(modelgroups{iGroup}, 1)
            model = modelgroups{iGroup}(iModel,:);
            if ~any(all(bsxfun(@eq, model, modellist), 2))
                modellist(end+1,:) = model;
            end
        end
    end    
end

for i = 1:size(modellist, 1)
    for j = 1:size(models, 1)
        if all(models(j, :) == modellist(i, :)); mf(j) = true; pos(i) = j; end
    end
end

% Get model names
if isfield(modelsummary, 'modelnames') && ~isempty(modelsummary.modelnames)
    modelnames = modelsummary.modelnames;
else
    for i = 1:size(models, 1); modelnames{i} = ['M' num2str(i)]; end
end 
newmodelnames = [];
for i = 1:length(modelnames)
    if mf(i); newmodelnames{end+1} = modelnames{i}; end
end

models = models(mf, :);
tab = tab(:, mf);

% Give sparse prior over alpha0
if isnumeric(bms_alpha0) && (isinf(bms_alpha0) || isnan(bms_alpha0))
    bms_alpha0 = 0.285 * size(tab, 2)^(-0.3) - 0.033;
end

if isstruct(bms_alpha0) && isfield(bms_alpha0, 'alpha0')
    bms_precomputed = 1;
    alpha = bms_alpha0.alpha;
    exp_r = bms_alpha0.exp_r;
    xp = bms_alpha0.xp;
    g = bms_alpha0.g;
    bms_smpl = bms_alpha0.smpl;
    models = bms_alpha0.models;
    ssorder = bms_alpha0.ssorder;
    newmodelnames = bms_alpha0.modelnames;
    bms_alpha0 = bms_alpha0.alpha0;
    tab = -2*log(g);
elseif bms_alpha0 > 0
    if bms_alpha0 < 1; str = ' (sparse prior)'; else str = []; end
    display(['Performing Bayesian Model Selection (BMS) with alpha_0 = ' num2str(bms_alpha0) str '.']);
    switch lower(score)
        case {'aic', 'aicc', 'bic', 'dic', 'waic1', 'waic2'}
            tab_adj = -0.5*tab;
        otherwise
            tab_adj = tab;
    end
    if bmsalgorithm == 1
        [alpha,exp_r,xp,g] = spm_BMS(tab_adj, [], [], [], [], bms_alpha0);
        bms_smpl = [];
    else
        % [exp_r,xp,g,bms_smpl] = gen_BMS_alpha(tab_adj,[],'sparse',[]);
        [exp_r,xp,g,bms_smpl] = gen_BMS(tab_adj,[],'sparse',[]);
        alpha = exp_r;
    end
    tab = -2*log(g);
    clear tab_adj;
elseif strcmpi(score, 'marginallike')
    tab = -2*tab;
end

if bmsorder && bms_alpha0 > 0
    % Order models according to posterior expectation of model probability
    [~, order] = sort(exp_r, 'descend');    
else
    % Order models according to average model score
    [~, order] = sort(nansum(tab, 1), 'ascend');
end
tab = tab(:, order);
models = models(order, :);

if bms_alpha0 > 0
    bms.alpha0 = bms_alpha0;
    bms.alpha = alpha(order);
    bms.exp_r = exp_r(order);
    bms.xp = xp(order);
    bms.g = g(:,order);
    bms.smpl = bms_smpl;
    if ~isempty(bms.smpl)
        bms.smpl.r = bms.smpl.r(:,order);
        if isfield(bms.smpl,'z'); bms.smpl.z = bms.smpl.z(:,order); end
        if isfield(bms.smpl,'u'); bms.smpl.u = bms.smpl.u(:,order); end
    end
    alphatot = sum(bms.alpha); % Compute SD of r
    bms.std_r = sqrt(bms.alpha.*(alphatot - bms.alpha)./(alphatot^2*(alphatot+1)));
    tab = -2*log(bms.g);
else
    bms = [];
end

% Reorder subjects
if isscalar(ssorder) && ~ssorder
    ssorder = 1:size(tab, 1);
elseif isscalar(ssorder)
    basetemp = bsxfun(@minus, tab, min(tab, [], 2));
    [bestscore,orderscore] = sort(basetemp,2,'ascend');
    temp = orderscore(:, 1)*10^12;
    if size(bestscore, 2) > 1; temp = temp + bestscore(:, 2).*orderscore(:, 2)*10^6; end
    if size(bestscore, 2) > 2; temp = temp + bestscore(:, 3).*orderscore(:, 3); end
    [~,ssorder] = sort(temp, 1, 'ascend');    
end

if bms_alpha0 > 0
    if ~bms_precomputed; bms.g = bms.g(ssorder,:); end
    bms.models = models;
    bms.ssorder = ssorder;
end
if ~bms_precomputed; tab = tab(ssorder, :); end

% Get model strings
newmodelnames
size(modellist)
size(newmodelnames)
size(order)
max(order)
for i = 1:size(modellist,1); modelstr{i} = newmodelnames{order(i)}; end
maxscore = scorebounds(1);

if bms_alpha0 > 0; bms.modelnames = modelstr; end


hold on;
% Plot scores
switch plottype
    case 1 % Patch table
        % Remove each individual subject's best model
        tab = bsxfun(@minus, tab, min(tab, [], 2));

        % Cap maximum difference, for visualization
        tab(tab > maxscore) = maxscore;

        if size(tab, 2) > MAXMODELS; tab = tab(:, 1:MAXMODELS); end
        
        ss = 1:size(tab, 1);
        for i = 1:size(tab, 1)
            for j = 1:size(tab, 2)
                patch([i i i+1 i+1], [-j -j-1 -j-1 -j], tab(i, j), 'EdgeColor', 'none');
            end
        end

        % Best models
        for i = 1:size(tab, 1)
            [besttab1, bestmodel1] = sort(tab(i, :));
            besttab1 = besttab1 - besttab1(1);
            besttab1 = find(besttab1 < 10);
            for j = 1:length(besttab1)
                text(i + 0.5, -bestmodel1(besttab1(j))-0.5, num2str(j), 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', axesfontsize, 'Color', 0.8*[1 1 1]);
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
        axis([1 (length(ss)+1), (-size(tab,2)-1) -1, -height height]);

        % Axes
        set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial');

        for i = 1:length(ss); xticklabel{i} = num2str(ssorder(i)); end
        set(gca, 'Xtick', 1.5:1:length(ss)+1, 'XtickLabel', xticklabel);

        for i = 1:size(tab,2)
            yticklabel{size(tab,2)-i+1} = modelstr{i};
        end
        set(gca, 'Ytick', (-size(tab,2):-1)-0.5, 'YtickLabel', yticklabel);

        xlabel('Subject number', 'FontName', 'Arial', 'FontSize', fontsize);

    case {2, 3} % Bar graph
        tab = bsxfun(@minus, tab, tab(:, 1));
        tabbase = tab;
        
        % Standard group-average model comparison
        if plottype == 2 || isempty(bms)
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
            pval = 1 - bms.xp;
            xrlimit = 1; % x-axis limit
        end
        
        if size(tab, 2) > MAXMODELS; tab = tab(:, 1:MAXMODELS); end

        % Remove ticks
        % patch([0.1 20 20 0.1], [-size(tab, 2) -size(tab, 2) -0.6 -0.6], [1 1 1], 'EdgeColor', 'none');

        h = barh(-1:-1:-size(tab, 2), tab');
        colmap = colormap;
        % colmap = flipud(colormap);
        for i = 1:size(tab, 2)
            hb(i) = barh(-i, tab(i), 0.8);
            colindex = 1 + floor(min([tab(i)/maxscore, 1])*(size(colmap, 1)-1));        
            set(hb(i), 'FaceColor', colmap(colindex, :));
        end
        
        % Plot errors
        for i = 1:size(tab, 1)
            for j = 1:size(tab, 2)
                if isnan(hb(i, j)); continue; end
                y = get(get(hb(i, j), 'children'), 'ydata');
                y = fliplr(mean(y([1 3], :)));
                plot([1 1]'*tab(i, j) + [0 1]'*tabse(i, j), [1 1]'*y, 'k', 'LineWidth', 1, 'Clipping', 'off');
                if tab(i, j) > maxscore*0.1 && tab(i, j) < maxscore*0.7; col2 = [0 0 0]; else col2 = 0.8*[1 1 1]; end
                plot([1 1]'*tab(i, j) + [-1 0]'*tabse(i, j), [1 1]'*y, 'Color', col2, 'LineWidth', 1, 'Clipping', 'off');
            end
        end    

        % Plot pvalues
        for i = 1:size(tab, 2)
            if pval(i) > 0.05 || isnan(pval(i)); continue; end
            if pval(i) < 0.001; ptext = '***'; elseif pval(i) < 0.01; ptext = '**'; else ptext = '*'; end 
            text(xrlimit*0.9, -i, ptext, 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', fontsize);        
%            text(tab(i) + 30, -i, ptext, 'HorizontalAlignment', 'Left', 'FontName', 'Arial', 'FontSize', fontsize);        
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
            [bestscore1, bestmodel1] = sort(tab(i, :));
            bestscore1 = bestscore1 - bestscore1(1);
            bestscore1 = find(bestscore1 < 10);
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
            xlabel(['{\Delta}' upper(score)], 'FontName', 'Arial', 'FontSize', fontsize, 'Interpreter', 'Tex');
        elseif plottype == 3
            xlabel('Posterior EMF', 'FontName', 'Arial', 'FontSize', fontsize, 'Interpreter', 'Tex');            
        end
        
    case 4 % Model features bar graph
        tab = bsxfun(@minus, tab, tab(:, 1));
        tabbase = tab;
        
        horizontalplot = 0;
                
        switch bmsalgorithm
            case 1 % Add alpha levels of each model to each feature
                alpha = zeros(1, length(modelgroups));
                for iModel = 1:size(models, 1)
                    model = models(iModel,:);
                    flag = zeros(1, length(modelgroups));
                    for iGroup = 1:length(modelgroups)
                        if any(all(bsxfun(@eq, model, modelgroups{iGroup}), 2))
                            flag(iGroup) = 1;
                        end
                    end
                    % Sum contributions of current model
                    alpha = alpha + bms.alpha(iModel).*flag/sum(flag);
                end
        
                % Compute posterior expected feature frequency and SD
                tab = alpha./sum(alpha);
                tabse = sqrt(alpha.*(sum(alpha)-alpha)/(sum(alpha)^2*(sum(alpha)+1)));
        
                nonfeatures = (tab == 0);
                tab(nonfeatures) = []; tabse(nonfeatures) = [];

                % Compute excess probability (i.e., the probability of being the 
                % most likely component)
                nsmpl = 10^6;
                [~,iMax] = max(drchrnd(alpha,nsmpl),[],2);
                for iGroup = 1:length(modelgroups)
                    xp(iGroup) = sum(iMax == iGroup)/nsmpl;
                end        
                pval = 1 - xp;
                xrlimit = 1; % x-axis limit
                
            case 2 % Sampling
                
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
        end
        
        if horizontalplot; h = barh(-1:-1:-size(tab, 2), tab'); else h = bar(1:size(tab, 2), tab'); end
        colmap = colormap;
        for i = 1:size(tab, 2)
            if horizontalplot; hb(i) = barh(-i, tab(i), 0.8); else hb(i) = bar(i, tab(i), 0.8); end
            % colindex = 1 + floor(min([tab(i)/maxscore, 1])*(size(colmap, 1)-1));        
            % set(hb(i), 'FaceColor', colmap(colindex, :));
            set(hb(i), 'FaceColor', 'w');
        end
        
        % Plot errors
        for i = 1:size(tab, 1)
            for j = 1:size(tab, 2)
                if isnan(hb(i, j)); continue; end
                y = get(get(hb(i, j), 'children'), 'ydata');
                y = fliplr(mean(y([1 3], :)));
                if horizontalplot
                    plot([1 1]'*tab(i, j) + [0 1]'*tabse(i, j), [1 1]'*y, 'k', 'LineWidth', 0.5, 'Clipping', 'off');
                else
                    plot([1 1]'*j, [1 1]'*tab(i, j) + [0 1]'*tabse(i, j), 'k', 'LineWidth', 0.5, 'Clipping', 'off');                        
                end
                   
                if tab(i, j) > maxscore*0.1 && tab(i, j) < maxscore*0.7; col2 = [0 0 0]; else col2 = 0.8*[1 1 1]; end
                if horizontalplot
                    plot([1 1]'*tab(i, j) + [-1 0]'*tabse(i, j), [1 1]'*y, 'Color', col2, 'LineWidth', 0.5, 'Clipping', 'off');
                else
                    plot([1 1]'*j, [1 1]'*tab(i, j) + [-1 0]'*tabse(i, j), 'Color', col2, 'LineWidth', 0.5, 'Clipping', 'on');                    
                end
            end
        end    

        % Plot pvalues
        for i = 1:size(tab, 2)
            if pval(i) > 0.05 || isnan(pval(i)); continue; end
            if pval(i) < 0.001; ptext = '***'; elseif pval(i) < 0.01; ptext = '**'; else ptext = '*'; end 
            if horizontalplot
                text(xrlimit*0.9, -i, ptext, 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', fontsize);
            else
                text(i, xrlimit*0.95, ptext, 'HorizontalAlignment', 'Center', 'FontName', 'Arial', 'FontSize', fontsize);                
            end
%            text(tab(i) + 30, -i, ptext, 'HorizontalAlignment', 'Left', 'FontName', 'Arial', 'FontSize', fontsize);        
        end

        % colormap;
        % colormap(flipud(colormap));
        if horizontalplot; axis([0 xrlimit, (-size(tab,2)-0.5) -0.5]); else axis([0.5 size(tab,2)+0.5, 0 xrlimit]); end

        % Axes
        set(gca, 'FontSize', axesfontsize, 'FontName', 'Arial','TickDir','out','TickLength',2*get(gca,'TickLength'));
        % box on;

        if isfield(modelsummary,'groupnames');
            str = [];
            modelsummary.groupnames
            for i = 1:length(nonfeatures)
                if ~nonfeatures(i); str{end+1} = modelsummary.groupnames{i}; end
            end
            modelstr = str;
        else
            for iGroup = 1:length(modelgroups)
                modelstr{iGroup} = ['F#' num2str(iGroup)];
            end
        end
        
        if horizontalplot
            for i = 1:size(tab,2)
                yticklabel{size(tab,2)-i+1} = modelstr{i};
            end
            set(gca, 'Xtick', 0:0.2:1, 'Ytick', (-size(tab,2):-1), 'YtickLabel', yticklabel);
            xlabel('Posterior Model Feature Frequency', 'FontName', 'Arial', 'FontSize', fontsize, 'Interpreter', 'Tex');
        else
            for i = 1:size(tab,2)
                xticklabel{i} = modelstr{i};
            end
            set(gca, 'Ytick', 0:0.2:1, 'Xtick', 1:size(tab,2), 'XtickLabel', xticklabel);
            xticklabel_rotate([],45);
            ylabel('Posterior Model Feature Frequency', 'FontName', 'Arial', 'FontSize', fontsize);
        end
end

end

function r = drchrnd(alpha,n)
% DRCHRND Take a sample from a dirichlet distribution.
alpha = alpha(:)';
r = gamrnd(repmat(alpha,n,1),1,[n,size(alpha,2)]);
% r = gamrnd(repmat(alpha,n,1),1,[n,size(alpha,2)]);
r = bsxfun(@rdivide, r, sum(r,2));
end
