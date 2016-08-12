function ModelPlot_drawPanel(panel)

fontsize = assign(panel, 'fontsize', 16);
axesfontsize = assign(panel, 'axesfontsize', 12);

hold on;
herr = []; % Handle for error bars
herrsize = [];

% Make priority list
for iPlot = 1:numel(panel.plots)
    panel.plots{iPlot}.priority = assign(panel.plots{iPlot},'priority',1);
    priority(iPlot) = panel.plots{iPlot}.priority;
end
priorityList = unique(priority);

for iPriority = 1:numel(priorityList)
    for iPlot = 1:length(panel.plots)
        thisplot = panel.plots{iPlot};
        if thisplot.priority ~= priorityList(iPriority); continue; end

        switch lower(thisplot.type)
            case {'dot','dots'}
                x = thisplot.x(:);
                y = thisplot.y(:);
                yerr = assign(thisplot, 'yerr', []); yerr = yerr(:);
                [x,y,yerr] = removeinvalidpoints(x,y,yerr);

                linecolor = assign(thisplot, 'linecolor', [0 0 0]);
                color = assign(thisplot, 'color', [0 0 0]);
                linewidth = assign(thisplot, 'linewidth', 1);
                if ~isempty(yerr)
                    % herr(end+1) = errorbar(x, y, yerr, 'k', 'LineWidth', linewidth,'LineStyle','none');
                    for i = 1:numel(x)
                        herr(end+1) = plot(x(i)*[1 1], y(i)+yerr(i)*[-1,1], '-', 'Color', linecolor, 'LineWidth', linewidth);
                    end
                    herrsize(end+1) = assign(thisplot, 'errorbar_size', 80);
                end
                markertype = assign(thisplot, 'markertype', 'o');
                markersize = assign(thisplot, 'markersize', 6);

                if ~isempty(y) && ~any(strcmpi(markertype,{'','none'}))
                    plot(x, y, markertype, 'MarkerSize', markersize, 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'LineWidth', linewidth);
                end
                
            case {'errorbar','errorbars'}
                x = thisplot.x(:);
                y = thisplot.y(:);
                yerr = assign(thisplot, 'yerr', []); yerr = yerr(:);
                [x,y,yerr] = removeinvalidpoints(x,y,yerr);

                % linecolor = assign(thisplot, 'linecolor', [0 0 0]);
                color = assign(thisplot, 'color', [0 0 0]);
                linewidth = assign(thisplot, 'linewidth', 1);
                errorwidth = assign(thisplot, 'errorwidth', []);
                if ~isempty(yerr)
                    if isempty(errorwidth)
                        herr(end+1) = errorbar(x, y, yerr, 'Color', color, 'LineWidth', linewidth,'LineStyle','-');
                        herrsize(end+1) = assign(thisplot, 'errorbar_size', 80);
                    elseif errorwidth == 0                        
                        for j = 1:numel(x)
                            herr(end+1) = plot(x(j)*[1 1], y(j)+yerr(j)*[-1 1], 'Color', color, 'LineWidth', linewidth,'LineStyle','-');
                        end
                        herr(end+1) = plot(x, y, 'Color', color, 'LineWidth', linewidth,'LineStyle','-');
                        herrsize(end+1) = assign(thisplot, 'errorbar_size', 80);
                        
                    end
                end

            case 'line'
                x = thisplot.x(:);
                y = thisplot.y(:);
                yerr = assign(thisplot, 'yerr', []); yerr = yerr(:);
                [x,y,yerr] = removeinvalidpoints(x,y,yerr);
                
                interp = assign(thisplot, 'interp', 0);
                if interp
                    [x,ord] = sort(x);
                    y = y(ord);
                    if ~isempty(yerr); yerr = yerr(ord); end                    
                    xnew = linspace(x(1),x(end),1e3);
                    y = interp1(x,y,xnew,'pchip');
                    if ~isempty(yerr); yerr = interp1(x,yerr,xnew,'spline'); end
                    x = xnew;
                end
                
                if numel(x) == 1
                    x = x + [-0.5 0.5];
                    y = y*[1 1];
                    yerr = yerr*[1 1];
                end
                color = assign(thisplot, 'color', [0 0 0]);
                errColor = assign(thisplot, 'errColor', color);
                linewidth = assign(thisplot, 'linewidth', 2);
                linestyle = assign(thisplot, 'linestyle', '-');
                facealpha = assign(thisplot, 'facealpha', 1);
                if isempty(yerr) || all(isnan(errColor))
                    plot(x,y,'Color',color,'LineWidth',linewidth,'LineStyle',linestyle);
                else
                    if ~isempty(y)
                        idx = isfinite(x) & isfinite(y) & isfinite(yerr);
                        x = x(idx); y = y(idx); yerr = yerr(idx);                        
                        
                        xx = [x(:)', fliplr(x(:)')];
                        yy = [y(:)' + yerr(:)', fliplr(y(:)' - yerr(:)')];
                        
                        if all(isnan(color))
                            h = fill(xx, yy, errColor, 'EdgeColor', 'none');
                            set(h,'FaceAlpha',facealpha);
                            % shadedErrorBar(x,y,yerr,{'Color',errColor,'LineWidth',linewidth,'LineStyle','none'},1,0,1);
                        else
                            h = fill(xx, yy, errColor, 'EdgeColor', color);                 
                            set(h,'FaceAlpha',facealpha);
                            % shadedErrorBar(x,y,yerr,{'Color',errColor,'LineWidth',linewidth,'LineStyle',linestyle},1,0,1);
                        end
                    end
                end

            case 'checkers2d'
                x = thisplot.x(:);
                y = thisplot.y(:);
                z = reshape(thisplot.z(:),length(x),length(y))';            
                xLim = assign(panel,'xLim',get(gca, 'XLim'));
                yLim = assign(panel,'yLim',get(gca, 'YLim'));
                zLim = assign(panel,'zLim',[nanmin(z(:)),nanmax(z(:))]);
                temp = max(min((z-zLim(1))/diff(zLim),1),0);
                z(~isnan(z)) = temp(~isnan(z));
                color = assign(thisplot, 'color', [0.5 0.5 0.5]);
                bkgcolor = assign(thisplot, 'bkgcolor', 0.95*[1 1 1]);            
                linestyle = assign(thisplot, 'linestyle', 'none');            
                % surf(x,y,z,'EdgeColor','none','LineStyle','none');
                dx = diff(xLim)/length(x);
                dy = diff(yLim)/length(y);
                a = 1;
                for ix = 1:length(x)
                    for iy = 1:length(y)
                        a = sqrt(z(ix,iy));
                        if ~isnan(z(ix,iy))
                            patch(0.5*dx*[-1 -1 1 1] + x(ix),0.5*dy*[1 -1 -1 1] + y(iy),0.95*[1 1 1],'FaceColor',bkgcolor,'LineStyle','none'); 
                            %patch(0.5*a*dx*[-1 -1 1 1] + x(ix),0.5*a*dy*[1 -1 -1 1] + y(iy),color,'LineStyle','none'); 
                            patch(0.5*a*dx*[-1 -1 1 1] + x(ix),0.5*a*dy*[1 -1 -1 1] + y(iy),'k','FaceColor',color,'LineStyle',linestyle,'LineWidth',1); 
                        end
                    end
                end

                view([0 90]);
                colormap(gray);

                %linewidth = assign(thisplot, 'linewidth', 1);
                %if ~isempty(yerr)
                %    herr(end+1) = errorbar(x, y, yerr, 'k', 'LineWidth', linewidth,'LineStyle','none');
                %    herrsize(end+1) = assign(thisplot, 'errorbar_size', 80);
                %end                        
                % plot3(x, y, ones(size(x)), markertype, 'MarkerSize', markersize, 'MarkerFaceColor', color, 'MarkerEdgeColor', color, 'LineWidth', linewidth);

            case 'bars2d'
                x = thisplot.x(:);
                y = thisplot.y(:);
                z = reshape(thisplot.z(:),length(x),length(y))';            
                xLim = assign(panel,'xLim',get(gca, 'XLim'));
                yLim = assign(panel,'yLim',get(gca, 'YLim'));
                zLim = assign(panel,'zLim',[nanmin(z(:)),nanmax(z(:))]);
                temp = max(min((z-zLim(1))/diff(zLim),1),0);
                z(~isnan(z)) = temp(~isnan(z));
                color = assign(thisplot, 'color', [0.5 0.5 0.5]);
                bkgcolor = assign(thisplot, 'bkgcolor', 0.95*[1 1 1]);            
                linestyle = assign(thisplot, 'linestyle', 'none');
                dx = diff(xLim)/length(x);
                shift = 0.5*dx*thisplot.shift2d;
                dy = diff(yLim)/length(y);
                for ix = 1:length(x)
                    for iy = 1:length(y)
                        % a = sqrt(z(ix,iy));
                        a = (z(ix,iy));                        
                        if ~isnan(z(ix,iy))
                            patch(0.5*dx*[-1 -1 1 1] + x(ix),0.5*dy*[1 -1 -1 1] + y(iy),0.95*[1 1 1],'FaceColor',bkgcolor,'LineStyle','none'); 
                            patch(0.25*dx*[-1 -1 1 1] + x(ix) + shift,0.5*a*dy*[1 -1 -1 1] + y(iy),'k','FaceColor',color,'LineStyle',linestyle); 
%                            patch(0.25*a*dx*[-1 -1 1 1] + x(ix) + shift,0.5*dy*[1 -1 -1 1] + y(iy),'k','FaceColor',color,'LineStyle',linestyle,'LineWidth',1); 
                        end
                    end
                end

                view([0 90]);
                colormap(gray);
                
        end
    end
end

% Increase tick length
ticklengthmultiplier = assign(panel,'ticklengthmultiplier',1);
set(gca, 'TickLength', ticklengthmultiplier*get(gca,'TickLength'));

xlabelstring = assign(panel,'xlabel',[]);
xlabelinterpreter = assign(panel,'xlabelinterpreter','none');
xlabel(xlabelstring,'FontSize',fontsize,'Interpreter',xlabelinterpreter);

ylabelstring = assign(panel, 'ylabel', []);
ylabelinterpreter = assign(panel,'ylabelinterpreter','none');
ylabel(ylabelstring, 'FontSize', fontsize,'Interpreter',ylabelinterpreter);

titlestring = assign(panel, 'title', []);
titleinterpreter = assign(panel,'titleinterpreter','none');
title(titlestring, 'FontSize', fontsize,'Interpreter',titleinterpreter);

xTick = assign(panel, 'xTick', get(gca, 'xTick'));
yTick = assign(panel, 'yTick', get(gca, 'yTick'));
xTickLabel = assign(panel, 'xTickLabel', []);
if isempty(xTickLabel); for i = 1:length(xTick); xTickLabel{i} = num2str(xTick(i)); end; end
yTickLabel = assign(panel, 'yTickLabel', []);
if isempty(yTickLabel); for i = 1:length(yTick); yTickLabel{i} = num2str(yTick(i)); end; end
set(gca, 'xTick', xTick, 'yTick', yTick, 'xTickLabel', xTickLabel, 'yTickLabel', yTickLabel);

boxflag = assign(panel, 'box', 0);
if boxflag; box on; else box off; end

set(gca, 'FontSize', axesfontsize, 'TickDir', 'out', 'Color', 'none');

xLim = assign(panel,'xLim',get(gca, 'XLim'));
yLim = assign(panel,'yLim',get(gca, 'YLim'));
axis([xLim yLim]);

% Plot zero line
plotzero = assign(panel, 'plotzero', 1);
% if plotzero; plot3(xLim, [0 0], [-1 -1], 'k:', 'LineWidth', 0.5); end
if plotzero; plot(xLim, [0 0], 'k:', 'LineWidth', 0.5); end

% Plot diagonal
plotdiagonal = assign(panel, 'plotdiagonal', 0);
% if plotdiagonal; plot3(xLim, yLim, [-1 -1], 'k:', 'LineWidth', 0.5); end

% Plot diagonal
axissquare = assign(panel, 'axissquare', 0);
if axissquare; axis square; end

% Resize error bars depending on version (need to fix)
d = version('-date');
if str2double(d(end-3:end)) < 2015
    for i = 1:length(herr); errorbar_tick(herr(i), herrsize(i)); end
end

end

%--------------------------------------------------------------------------
%ASSIGN Check if a struct field exists and if nonempty return its value, 
%       otherwise return default
function value = assign(this, field, default)
    if isfield(this, field) && ~isempty(this.(field))
        value = this.(field);
    else
        value = default;
    end
end

%--------------------------------------------------------------------------
function [x,y,yerr] = removeinvalidpoints(x,y,yerr)
%REMOVEINVALIDPOINTS Remove Infs and NaNs
    idx = ~isfinite(x) | ~isfinite(y);
    x(idx) = [];
    y(idx) = [];
    if ~isempty(yerr); yerr(idx) = []; end
end
