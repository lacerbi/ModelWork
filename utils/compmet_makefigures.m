function compmet_makefigures(n)
%COMPATIBILITY_MAKEFIGURES Make figures for compatibility statistic paper.

fontsize = 18;
axesfontsize = 14;

switch n
    case 1
        
        figsize = [383 262 1157 682];
        hg = plotify(1,2,'Margins',[0.15 0.05 0.15 0.05],'Gutter',[0.1 0.05],'Labels',{'a','b'});
        
        col = [0 0 0; 0 0 0.5; 0.1 0.1 0.7; 0.3 0.3 0.9; 0.5 0.5 1];
        
        Nx = 2^10;
        xx = linspace(-5,5,Nx);
        pdf1 = exp(-0.5*xx.^2);
                
        mu = linspace(0,6,1e3);
        mu2 = [0 2 4 6];
        
        C = zeros(1+numel(mu2),numel(mu));
        for i = 1:numel(mu)
            pdf2 = exp(-0.5*(xx-mu(i)).^2);
            C(1,i) = compmet(pdf1,pdf2,'pdf');
            
            for j = 1:numel(mu2)
                pdf3 = exp(-0.5*(xx-mu2(j)).^2);
                C(1+j,i) = compmet(pdf1,pdf2,pdf3,'pdf');
            end
        end
        
        for i = 1:2
            axes(hg(i));
            
            % thresh = [0.05 0.01];
            thresh = [1/21, 1/151];
            
            % log Bayes factor
            if i == 1
                fun = @(x) x;
                yticks = [0,1/21,0.2:0.2:1];
                for j = 1:numel(yticks); ytickslabel{j} = num2str(yticks(j)); end 
                ytickslabel{2} = num2str(thresh(1),'%.3f');
                ybounds = [0 1];
                ystring = '$\gamma_{\theta}$';
            else
                fun = @(x) log((1-x)./x);
                yticks = sort([-4:2:8,fun(thresh(1)),fun(thresh(2))]);
                ytickslabel = {'-4','-2','0','2','log(20)','4','log(150)','6','8'};
                ybounds = [-3 7];
                ystring = '$\log_e B_{10}$';
            end
            
            plot([mu(1), mu(end)], fun(thresh(1)*[1 1]), ':k','LineWidth', 1); hold on;
            plot([mu(1), mu(end)], fun(thresh(2)*[1 1]), ':k','LineWidth', 1); hold on;
            h(1) = plot(mu,fun(C(1,:)),'-','Color',col(1,:), 'LineWidth',3); hold on;
            for j = 1:numel(mu2)
                h(j+1) = plot(mu,fun(C(j+1,:)),'-','Color',col(j+1,:),'LineWidth',2);
            end
                     
            set(gca,'TickDir','out','XTick',0:6,'YTick',yticks,'Yticklabel',ytickslabel,'FontSize',axesfontsize);
            axis([mu(1), mu(end), ybounds]);
            box off;
            xlabel('$\mu_2$','Interpreter','LaTeX','FontSize',fontsize);
            ylabel(ystring,'Interpreter','LaTeX','FontSize',fontsize);
            if i == 1
                hl = legend(h,'$n = 2$', '$n = 3$, $\mu_3 = 0$', '$n = 3$, $\mu_3 = 2$', '$n = 3$, $\mu_3 = 4$', '$n = 3$, $\mu_3 = 6$');
                set(hl,'Box','off','Location','NorthEast','Interpreter','LaTeX');
            end
        end
        
        set(gcf,'Color','w');
    
end

set(gcf,'Position',figsize);