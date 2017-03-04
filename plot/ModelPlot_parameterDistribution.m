function [postpdf,subjpdf,C,P,pxp,xx] = ModelPlot_parameterDistribution(mbag,paramname,bms,modelnames,modelsummary,plottype,color)
%MODELPLOT_PARAMETERDISTRIBUTION Plot distribution of parameter across models

if nargin < 4; modelnames = []; end
if nargin < 5; modelsummary = []; end
if nargin < 6 || isempty(plottype); plottype = 1; end
if nargin < 7; color = []; end

if iscell(mbag); cellflag = 1; else mbag = {mbag}; cellflag = 0; end
if ~iscell(bms); bms = {bms}; end

Nbags = numel(mbag);

Ngrid = 2^8;
% MaxSamples = Inf;
MaxSamples = 1000;

if isempty(modelnames); for m = 1:Nbags; modelnames{m} = []; end; end
if ischar(modelnames{1}); modelnames = {modelnames}; end
if isempty(modelsummary); for m = 1:Nbags; modelsummary{m} = []; end; end

if isempty(color)
    switch Nbags
        case 1; color = [0 0 0];
        case 2; color = [31,120,180; 178,223,138]/255;
        case 3; color = [31,120,180; 178,223,138; 166,206,227]/255;
        case 4; color = [31,120,180; 178,223,138; 166,206,227; 51,160,44]/255;        
    end
end

logflag = 0;

bounds = [Inf,-Inf];

for m = 1:Nbags

    if isempty(modelnames{m}); modellist = bms{m}.modelnames; end
    Nmodels = numel(modellist);
    
    if isfield(mbag{m},'legend')
        legendlabel{m} = mbag{m}.legend;
    else
        legendlabel{m} = [];
    end

    if isempty(modelsummary{m})
        modelsummary{m} = ModelWork_summary(mbag{m});
    end
    ssorder = bms{m}.ssorder;
    exp_r = zeros(1,1,Nmodels);
    g = [];
    pdf{m} = [];

    for i = 1:Nmodels
        mname = modellist{i};
        model_idx = find(strcmp(mname,bms{m}.modelnames));    
        if isempty(model_idx) || numel(model_idx) > 1
            error(['Model ' mname 'should be present and uniquely identified in the BMS structure.']);
        end

        mod_idx = find(strcmp(mname,modelsummary{m}.modelnames),1);    
        mfit = ModelBag_get(mbag{m},modelsummary{m}.dataid,modelsummary{m}.models(mod_idx,:),modelsummary{m}.cnd);
        Nsubjs = numel(mfit);

        if isempty(g); g = zeros(Nsubjs,1,Nmodels); end        
        if isempty(pdf{m}); pdf{m} = zeros(Nsubjs,Ngrid,Nmodels); end        

        param_idx = find(strcmp(paramname,mfit{1}.mp.params),1);

        if isempty(param_idx)
            exp_r(i) = 0;
            g(:,1,i) = 0;        
        else
            LB = mfit{1}.mp.bounds.LB(param_idx);
            UB = mfit{1}.mp.bounds.UB(param_idx);
            if isfield(mfit{1}.mp.bounds,'logflag')
                logflag = mfit{1}.mp.bounds.logflag(param_idx);
            end
            for j = 1:numel(mfit)
                Nsamples = min(MaxSamples,size(mfit{j}.sampling.samples,1));                
                idx = round(linspace(1,size(mfit{j}.sampling.samples,1),Nsamples));
                samples = mfit{j}.sampling.samples(idx,param_idx);
                if 0
                    [b,pdf(j,:,i),xx] = kde(samples,Ngrid,LB,UB);
                else
                    xx = linspace(LB,UB,Ngrid);
                    b = (max(samples)-min(samples))/100;
                    % w = ones(Nsamples,1)/Nsamples;
                    pdf{m}(j,:,i) = sum(exp(-0.5 * bsxfun(@minus, xx, samples).^2/b^2), 1)/ sqrt(2*pi) / (b*Nsamples);
                end

            end

            % Model posterior frequency
            exp_r(i) = bms{m}.exp_r(model_idx);
            g(:,1,i) = bms{m}.g(:,model_idx);        
        end

    end

    % Normalize posteriors
    exp_r = exp_r / sum(exp_r);
    g = bsxfun(@rdivide, g, sum(g,3));

    % Reorder subjects to external order
    g(ssorder,:) = g(1:end,:);

    dx = xx(2)-xx(1);
    subjpdf{m} = sum(bsxfun(@times, ...
        bsxfun(@times, pdf{m}, g), ...
        exp_r),3);
    subjpdf{m} = bsxfun(@rdivide, subjpdf{m}, qtrapz(subjpdf{m},2) * dx);

    postpdf{m} = mean(subjpdf{m},1);
    postcdf = qcumtrapz(postpdf{m})*dx;
    
    bounds(1) = min(bounds(1),xx(find(postcdf >= 0.001,1)));
    bounds(2) = max(bounds(2),xx(find(postcdf >= 0.999,1)));
 
    if bounds(1) > 0 && bounds(1) < 0.05; bounds(1) = 0; end
    
    % Plot individual posteriors
    lightcol = 0.4*color(m,:) + 0.6*[1 1 1];
    switch plottype
        case 1
            for i = 1:Nsubjs
                rescaling = max(postpdf{m}) / max(subjpdf{m}(i,:));
                plot(xx, subjpdf{m}(i,:) * rescaling, '-', 'Color', lightcol, 'LineWidth', 0.5); hold on;
            end
        case 2
            for i = 1:Nsubjs
                cdf = qcumtrapz(subjpdf{m}(i,:))*dx;
                qq(1) = xx(find(cdf >= 0.025,1));
                qq(2) = xx(find(cdf >= 0.25,1));
                qq(3) = xx(find(cdf >= 0.5,1));
                qq(4) = xx(find(cdf >= 0.75,1));
                qq(5) = xx(find(cdf >= 0.975,1));
                plot([qq(1),qq(5)],-[m m]-0.7*(i-(1+Nsubjs)/2)/Nsubjs,'-','Color', lightcol, 'LineWidth', 1); hold on;
                plot([qq(2),qq(4)],-[m m]-0.7*(i-(1+Nsubjs)/2)/Nsubjs,'-','Color', color(m,:), 'LineWidth', 2); hold on;
                
                %rescaling = max(postpdf{m}) / max(subjpdf{m}(i,:));
                %plot(xx, subjpdf{m}(i,:) * rescaling, '-', 'Color', lightcol, 'LineWidth', 0.5); hold on;
            end            
    end
end

if plottype == 1
    % Plot mean posteriors
    for m = 1:Nbags
        h(m) = plot(xx, postpdf{m}, '-', 'Color', color(m,:), 'LineWidth', 3);
    end
end

% Compute compatibility
if Nbags > 1 && nargout > 2
    C = zeros(1,Nsubjs);
    if 0
        for j = 1:Nsubjs
            X = [];
            for m = 1:Nbags; X{m} = subjpdf{m}(j,:); end        
            [C(j),P(j,:)] = compatibility(X{:},'pdf');
        end
    else        
        for m = 1:Nbags; X{m} = subjpdf{m}; end
        [pxp,C,P] = cons('pdf',X{:});        
    end
else
    C = []; P = [];
end

if plottype
    box off;
    set(gca,'TickDir','out');
    set(gcf, 'Color','w');
    paramname(paramname == '_') = '-';
    xlabel(paramname);
    
    if logflag
        xtick = [0.01 0.03 0.1 0.3 1 3 10 30 100];
        for iTick = 1:numel(xtick); xticklabel{iTick} = num2str(xtick(iTick)); end
        set(gca,'XTick',log(xtick),'XTickLabel', xticklabel);
    end
    xlim(bounds);
    
    switch plottype
        case 1
            hl = legend(h,legendlabel{:});
            set(hl,'Box','off','Location','NorthWest');
            ylabel('Average posterior density');
        case 2
            % ytick = -Nbags:1;
            ytick = [];
            set(gca,'YTick',ytick,'YTickLabel',[]);
            ylim([-Nbags-0.75,-0.25]);
            ylabel('Tasks');
    end
    drawnow;
end

% Remove cellification
if ~cellflag
    temp = pdf{1}; pdf = temp;
    temp = postpdf{1}; postpdf = temp;
end

end