function ModelWork_plotChains(mfit, paramn, plottype, thin)
% MODELWORK_PLOTCHAINS Plot sampled chains in model structures MFIT. 
% Can specify a parameter number PARAMN (otherwise plots all parameters in 
% different graphs). PLOTTYPE can be either 'pdf' or 'cdf' (default 'pdf').
%
% Common usage:
% ModelWork_plotChains(ModelBag_get(mbag,1:16,modelsummary.models(1,:),1))

if nargin < 1
    help ModelWork_plotChains;
    return;
end

% Wrap singleton model struct
if isfield(mfit, 'nid'); mfit = {mfit}; end

% Plot all parameters if not specified
if ~exist('paramn', 'var'); paramn = []; end
if isempty(paramn); paramn = [1:size(mfit{1}.smpl, 2) 0]; end

% Plot type
if ~exist('plottype', 'var'); plottype = []; end
if isempty(plottype); plottype = 'pdf'; end
switch lower(plottype)
    case 'pdf'; plottype = 1;
    case 'cdf'; plottype = 2;
    otherwise; error('Unknown plot type');
end

% Thinning factor (by default do not thin)
if ~exist('thin', 'var'); thin = []; end
if isempty(thin); thin = 1; end


% Only one subject -- plot all parameters in a single figure
if length(mfit) == 1
    figure(mfit{1}.nid);

    sep = [0.05, 0.1];
    z = ones(5, 1);
    gg = []; gbottom = []; gleft = 1:4;
    for j = 1:ceil(length(paramn)/4)
        gg = [gg; z*((1:4)+(j-1)*4)];
        gbottom = [gbottom; j*4];
    end
    hg = multigraph(gg', sep, [0.1, 0.1]);
    set(gcf,'color','w');

    for i = 1:length(paramn)
        plotParamChains(mfit, paramn(i), plottype, hg(i), gleft, gbottom, thin);
    end
    
    
else
    for i = 1:length(paramn)
        figure(i);
        
        sep = [0.02, 0.02];
        z = ones(5, 1);
        gg = []; gbottom = []; gleft = 1:4;
        for j = 1:ceil(length(mfit)/4)
            gg = [gg; z*((1:4)+(j-1)*4)];
            gbottom = [gbottom; j*4];
        end
        hg = multigraph(gg', sep, [0.1, 0.1]);
        set(gcf,'color','w');
                
        plotParamChains(mfit, paramn(i), plottype, hg, gleft, gbottom, thin);
        
        
    end
end


end

% Plot a single parameter for all subjects
function plotParamChains(mfit, pn, plottype, hg, gleft, gbottom, thin)

fontsize = 18;
    
xmin = NaN(1, length(mfit));
xmax = NaN(1, length(mfit));

% Take variable bounds
for i = 1:length(mfit)    
    m = mfit{i};    
    if isempty(m); continue; end
    if pn == 0; smpl = sum(m.loglikes, 2); else smpl = m.smpl(:, pn); end
    xmin(i) = quantile(smpl, 0.01);
    xmax(i) = quantile(smpl, 0.99);
%    xmin(i) = min(m.smpl(:, pn));
%    xmax(i) = max(m.smpl(:, pn));
end

xbounds = [quantile(xmin, 0.025), quantile(xmax, 0.975)];

for i = 1:length(mfit)
    axes(hg(i));
    hold on;
        
    m = mfit{i};
    if isempty(m); continue; end

    paramName = plotSubjectParamChains(m, pn, plottype, xbounds, fontsize, thin);
    % if all(i ~= gbottom); set(gca,'XTick', []); end
    if all(i ~= gleft); set(gca,'YTick', []); end    
end

% Parameter name
if length(mfit) == 1
    title(paramName, 'FontSize', fontsize);   
    % set(gca, 'Xtick', xbounds);
else
    sep = [0.02, 0.02];
    hg2 = multigraph(1, sep, [0.08, 0.08]);
    axes(hg2(1));
    title(['Parameter: ' paramName ], 'FontSize', fontsize);
    axis off;
end

end



function paramName = plotSubjectParamChains(m, pn, plottype, xbounds, fontsize, thin)
    maxp = 0;
    if pn == 0
        paramName = 'Log likelihood';
    elseif pn <= length(m.mp.params)
        paramName = m.mp.params{pn};
    else
        paramName = ['Parameter #' num2str(pn)];        
    end    
    
    if pn > size(m.smpl, 2); return; end

    for k = 1:m.nchains;
        nsamplesperchain = size(m.smpl, 1)/m.nchains;
        if pn == 0
            samples = sum(m.loglikes((thin:thin:nsamplesperchain) + (k-1)*nsamplesperchain, :), 2);
        else
            samples = m.smpl((thin:thin:nsamplesperchain) + (k-1)*nsamplesperchain, pn);
        end
        [~, pdf, x, cdf] = kde(samples, 1024);
        col = [1 1 1]*(k-1)/m.nchains;    
    
        if plottype == 1
            plot(x, pdf, 'k', 'LineWidth', 2, 'Color', col);
            maxp = max(max(pdf), maxp);
        else
            plot(x, cdf, 'k', 'LineWidth', 2, 'Color', col);
            maxp = 1;
        end
                
    end
    
    % Plot lower and upper bounds for initial conditions
    if pn > 0
        rminx = m.mp.bounds.RLB(pn);
        rmaxx = m.mp.bounds.RUB(pn);
        plot(rminx*[1 1], [0, maxp], 'k--', 'LineWidth', 1); % Lower bound
        plot(rmaxx*[1 1], [0, maxp], 'k--', 'LineWidth', 2); % Upper bound
    end
    
    h = text(0, 0, num2str(m.nid), 'HorizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', fontsize, 'FontName', 'Arial');
    set(h, 'Position', [0.15 0.85 0]);
    xlabel(''); ylabel(''); title('');
    
    axis([xbounds, 0 maxp]); 

end
