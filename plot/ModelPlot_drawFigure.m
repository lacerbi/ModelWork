function [fig,gendata] = ModelPlot_drawFigure(fig,data,mfit,ngen,options)

if nargin<3; mfit = []; end
if nargin<4; ngen = 30; end % Fake datasets per subject
if nargin<5; options = []; end

% figure();
set(gcf, 'Color', 'w');

panelgraph = assign(fig, 'panelgraph', 1:length(fig.panels)); 
intborder = assign(fig, 'intborder', [0.025, 0.1]);
extborder = assign(fig, 'extborder', [0.1, 0.1]);
if numel(fig.panels) > 1
    fig.hg = plotify(panelgraph, 'Gutter', intborder, 'Margins', extborder);
else
    fig.hg = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate fake datasets (this will be revisited)

if ~isempty(mfit)
    if isfield(fig,'project'); project = fig.project; else project = fig.prefix; end
    
    % Fake data generation function and analytics function
    gendataFun = str2func([project '_gendata']);
    analyticsFun = str2func([project '_analytics']); 
    gendata = []; % Fake data struct
    gendata.datamats = [];
    gendata.data = [];
    
    % Passing a previously generated GENDATA object
    if isfield(mfit, 'datamats')
        gendata = mfit;
    % Passing a single MFIT object, generate a number of fake datasets
    elseif isfield(mfit, 'mp')
        gendata.datamats{1} = gendataFun(ngen,mfit);
    % Passing an array of arrays of fake data matrices
    elseif iscell(mfit) && iscell(mfit{1}) && isnumeric(mfit{1}{1})
        gendata.datamats = mfit;        
    % Passing an array of MFIT objects
    elseif isfield(mfit{1}, 'mp') && length(mfit) > 1
        fprintf('Generating subjects'' datasets: ');
        for i = 1:length(mfit)
            fprintf('%d..', i);
            temp = gendataFun(ngen,mfit{i});
            if isnumeric(temp{1})
                gendata.datamats{i} = temp;
                gendata.data{i} = [];
            else
                gendata.data{i} = temp;
            end
        end
        fprintf('\n');
    end
    
    % Generate full datasets
    for i = 1:length(gendata.datamats)
        if isempty(gendata.data) || isempty(gendata.data{i})
            gendata.data{i} = analyticsFun(gendata.datamats{i},options);
        end
    end
    clear mfit temp;
else
    gendata = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare and draw figure

% Loop over panels
for iPanel = 1:length(fig.panels)
    if numel(fig.hg) > 1; axes(fig.hg(iPanel)); end
    thispanel = fig.panels{iPanel}; % Pointer for external functions
    
    % Loop over data plots inside panel
    for jPlot = 1:length(thispanel.plots)
        thisplot = thispanel.plots{jPlot}; % Pointer for external functions
        source = thisplot.source;
        nid = assign(source, 'nid', 0); % Default data source is all subjects
        
        switch lower(source.type)            
            case 'data' % Plot data
                thisplot = PrepareDataPlot(thisplot,data,source,nid);
                % Defaul plot type for data is dots with error bars
                thisplot.type = assign(thisplot, 'type', 'dots'); 

            case 'data2d' % Plot 2D data
                thisplot = PrepareDataPlot2D(thisplot,data,source,nid);
                % Defaul plot type for data is checkerboard matrix
                % thisplot.type = assign(thisplot, 'type', 'checkers2d');
                thisplot.type = assign(thisplot, 'type', 'bars2d');
                thisplot.linestyle = assign(thisplot,'linestyle','none');
                thisplot.shift2d = assign(thisplot,'shift2d',-0.5);                
                
            case {'model','model2d'} % Plot model fit
                if strcmpi(source.type,'model')
                    plotfields = {'x','y','yerr'};
                    prepareDataFun = @PrepareDataPlot; 
                elseif strcmpi(source.type,'model2d')
                    plotfields = {'x','y','z','zerr'};                    
                    prepareDataFun = @PrepareDataPlot2D; 
                end
                
                % source.binfunerr1 = '@(y) std(y)'; % SD of fake data    
                
                source.dataids = assign(source,'dataids',0);
                if source.dataids == 0; source.dataids = 1:numel(gendata.data); end
                
                ii = 1;
                for iGen = source.dataids                  
                    tempplot = prepareDataFun(thisplot,gendata.data{iGen},source,0,1);
                    for iField = 1:length(plotfields)                        
                        temp.(plotfields{iField})(ii,:) = tempplot.(plotfields{iField});
                    end
                    ii = ii + 1;
                end
                                
                thisplot = tempplot;
                for iField = 1:length(plotfields)
                    thisplot.(plotfields{iField}) = nanmean(temp.(plotfields{iField}),1);
                end
                if strcmpi(source.type,'model')
                    if size(temp.yerr, 1) == 1
                        thisplot.yerr = temp.yerr;
                    else
                        thisplot.yerr = stderr(temp.y,[],1);                  
                    end
                    thisplot.type = assign(thisplot, 'type', 'line');
                    thisplot.linestyle = assign(thisplot,'linestyle','none');
                elseif strcmpi(source.type,'model2d')
                    if size(temp.zerr, 1) == 1
                        thisplot.zerr = temp.zerr;
                    else
                        thisplot.zerr = stderr(temp.z,[],1);                    
                    end
                    thisplot.shift2d = assign(thisplot,'shift2d',0.5);
                    % thisplot.type = assign(thisplot,'type','checkers2d');
                    thisplot.type = assign(thisplot,'type','bars2d');
                    thisplot.color = assign(thisplot,'color',[0.2 0.2 0.2]);
                    % thisplot.bkgcolor = assign(thisplot,'bkgcolor',[0.2 0.2 0.2]);
                    % thisplot.linestyle = assign(thisplot,'linestyle','none');
                    thisplot.bkgcolor = assign(thisplot,'bkgcolor','none');
                    % thisplot.linestyle = assign(thisplot,'linestyle','-');
                end
                    
        end
        
        thispanel.plots{jPlot} = thisplot;
    end
    
    fig.panels{iPanel} = thispanel;
    
    ModelPlot_drawPanel(thispanel);
    
end

end

%ASSIGN Check if a struct field exists and if nonempty return its value, 
%       otherwise return default
function value = assign(this, field, default)
    if isfield(this, field) && ~isempty(this.(field))
        value = this.(field);
    else
        value = default;
    end
end

%PREPAREDATAPLOT Prepare a single data plot.
function thisplot = PrepareDataPlot(thisplot,data,source,nid,gendataflag)
    if nargin < 5; gendataflag = 0; end

    % NID 0 means use all subjects' data
    if nid == 0; nid = 1:length(data); end
    
    % Defaut plot method is binning
    method = assign(source, 'method', 'bin');                 

    % Initialize data for binning
    if strcmpi(method, 'bin') || strcmpi(method, 'hist')
        if isfield(source, 'bincenters') && ~isempty(source.bincenters)
            bincenters = source.bincenters;
        else
            x = [];
            for iSubj = 1:length(nid)
                thisdata = data{nid(iSubj)};
                temp = eval(source.x);
                x = [x; temp(:)];
            end
            [~,bincenters] = hist(x);                          
        end
        if numel(bincenters) > 1
            binbounds = [-Inf, bincenters(1:end-1) + 0.5*diff(bincenters), Inf];
        else
            binbounds = [-Inf, Inf];            
        end
        switch lower(method)
            case 'bin'
                binfun = str2func(assign(source, 'binfun', '@(y) nanmean(y)'));
                binfunerr1 = str2func(assign(source, 'binfunerr1', '@(y) stderr(y)'));
            case 'hist'
                binfun = @(y) length(y(:));
                binfunerr1 = @(y) 0;            
        end

        yy = []; yyerr1 = [];
    end
    xall = []; yall = [];

    % Collect data from all subjects
    for iSubj = 1:length(nid)
        thisdata = data{nid(iSubj)};
        
        try            
            yold = eval(source.y); % There must be a Y
            if isfield(source,'x') && ~isempty(source.x)
                xold = eval(source.x);
            else
                xold = yold;
            end
            x = xold; y = yold;
            if isfield(source,'xfun') && ~isempty(source.xfun)
                fun = str2func(source.xfun);
                x = fun(xold,yold);
            end
            if isfield(source,'yfun') && ~isempty(source.yfun)
                fun = str2func(source.yfun);
                y = fun(xold,yold);
            end
            
            % Add some uniform jitter to the x-axis
            xjitter = assign(thisplot, 'xjitter', 0);
            x = x + (rand(size(x))-0.5)*xjitter;
            % Add some uniform jitter to the y-axis
            yjitter = assign(thisplot, 'yjitter', 0);
            y = y + (rand(size(y))-0.5)*yjitter;
            
            switch lower(method)
                case {'bin','hist'} % Bin data or build histogram
                    [~,pos] = histc(x, binbounds);
                    for k = 1:length(bincenters)
                        % binfun(y(pos == k))
                        yy(iSubj,k) = binfun(y(pos == k));
                        if length(nid) == 1
                            yyerr1(k) = binfunerr1(y(pos == k));
                        end
                    end
                    % Normalize histograms
                    if strcmpi(method,'hist')
                        yy = bsxfun(@rdivide, yy, sum(yy,2));
                    end

                case 'pool'
                    xall = [xall; x(:)];
                    yall = [yall; y(:)];
            end
            
        catch
            display(['Cannot plot data for dataset #' num2str(nid(iSubj)) '.']);
            
            switch lower(method)
                case {'bin','hist'} % Bin data or build histogram
                    yy(iSubj,1:numel(bincenters)) = NaN;
                case 'pool'
            end
            
        end

    end
    
    % Store data
    switch lower(method)
        case {'bin','hist'}
            binshift = assign(thisplot, 'binshift', 0);
            thisplot.x = bincenters + binshift;
            thisplot.y = nanmean(yy,1);
            
            if length(nid) == 1
                if isfield(source, 'yerr') && ~isempty(source.yerr)
                    thisdata = data{1};
                    thisplot.yerr = eval(source.yerr);
                elseif ~isempty(yyerr1)
                    thisplot.yerr = yyerr1;
                else
                    thisplot.yerr = [];                        
                end
            else
                if gendataflag
                    thisplot.yerr = nanstd(yy,[],1);                 
                else
                    thisplot.yerr = stderr(yy,[],1);
                end
            end
        case 'pool'
            thisplot.x = xall;
            thisplot.y = yall;
            thisplot.yerr = [];                                
    end
        
end


%PREPAREDATAPLOT2D Prepare a single 2D data plot.
function thisplot = PrepareDataPlot2D(thisplot,data,source,nid,gendataflag)
    if nargin < 5; gendataflag = 0; end

    % NID 0 means use all subjects' data
    if nid == 0; nid = 1:length(data); end
    
    % Defaut plot method is binning
    method = assign(source, 'method', 'bin');                 

    % Initialize data for binning
    if strcmpi(method, 'bin') || strcmpi(method, 'hist')
        datafield = {'x','y'};
        for ibin = 1:2
            if isfield(source, [datafield{ibin} 'bincenters']) && ~isempty(source.([datafield{ibin} 'bincenters']))
                bincenters{ibin} = source.([datafield{ibin} 'bincenters']);
            else
                x = [];
                for iSubj = 1:length(nid)
                    thisdata = data{nid(iSubj)};
                    temp = eval(source.(datafield{ibin}));
                    x = [x; temp(:)];
                end
                [~,xbincenter{ibin}] = hist(x);
            end
            binbounds{ibin} = [-Inf, bincenters{ibin}(1:end-1) + 0.5*diff(bincenters{ibin}), Inf];
        end
        
        switch lower(method)
            case 'bin'
                binfun = str2func(assign(source, 'binfun', '@(y) nanmean(y)'));
                binfunerr1 = str2func(assign(source, 'binfunerr1', '@(y) stderr(y)'));
            case 'hist'
                binfun = @(y) length(y(:));
                binfunerr1 = @(y) 0;            
        end

        nxBins = length(binbounds{1})-1;
        nyBins = length(binbounds{2})-1;
        nBinsTot = nxBins*nyBins;
        zz = zeros(length(nid),nBinsTot); zzerr1 = zeros(1, nBinsTot);
    end
    xall = []; yall = []; zall = [];

    % Collect data from all subjects
    for iSubj = 1:length(nid)
        thisdata = data{nid(iSubj)};
        xold = eval(source.x);
        yold = eval(source.y);
        zold = eval(source.z);
        x = xold; y = yold; z = zold;
        if isfield(source,'xfun') && ~isempty(source.xfun)
            fun = str2func(source.xfun);
            x = fun(xold,yold,zold);
        else
            x = xold;
        end
        if isfield(source,'yfun') && ~isempty(source.yfun)
            fun = str2func(source.yfun);
            y = fun(xold,yold,zold);
        else
            y = yold;
        end
        if isfield(source,'zfun') && ~isempty(source.zfun)
            fun = str2func(source.zfun);
            z = fun(xold,yold,zold);
        else
            z = zold;
        end
        % Add some uniform jitter to the axes
        xjitter = assign(thisplot, 'xjitter', 0);
        x = x + (rand(size(x))-0.5)*xjitter;
        yjitter = assign(thisplot, 'yjitter', 0);
        y = y + (rand(size(y))-0.5)*yjitter;
        zjitter = assign(thisplot, 'zjitter', 0);
        z = z + (rand(size(z))-0.5)*zjitter;

        switch lower(method)
            case {'bin','hist'} % Bin data or build histogram
                [~,posx] = histc(x, binbounds{1});
                [~,posy] = histc(y, binbounds{2});
                for ix = 1:nxBins
                    for iy = 1:nyBins
                        index = (ix-1)*nxBins + iy;
                        zz(iSubj,index) = binfun(z(posx == ix & posy == iy));
                        if length(nid) == 1
                            zzerr1(index) = binfunerr1(z(posx == ix & posy == iy));
                        end
                    end
                end
                % Normalize histograms
                if strcmpi(method,'hist')
                    zz = bsxfun(@rdivide, zz, sum(zz,2));
                end
                
            case 'pool'
                error('Unknown method ''pool'' for 2D data.')
                xall = [xall; x(:)];
                yall = [yall; y(:)];
                zall = [zall; z(:)];
        end

    end    
    
    % Store data
    switch lower(method)
        case {'bin','hist'}
            xbinshift = assign(thisplot, 'xbinshift', 0);
            ybinshift = assign(thisplot, 'ybinshift', 0);
            thisplot.x = bincenters{1} + xbinshift;
            thisplot.y = bincenters{2} + ybinshift;
            thisplot.z = nanmean(zz,1);
                        
            if length(nid) == 1
                if isfield(source, 'zerr') && ~isempty(source.zerr)
                    thisdata = data{1};
                    thisplot.zerr = eval(source.zerr);
                elseif ~isempty(zzerr1)
                    thisplot.zerr = zzerr1;
                else
                    thisplot.zerr = [];                        
                end
            else
                if gendataflag
                    thisplot.zerr = nanstd(zz,[],1);                  
                else
                    thisplot.zerr = stderr(zz,[],1);
                end
            end
        case 'pool'
            thisplot.x = xall;
            thisplot.y = yall;
            thisplot.yerr = [];                                
    end

end
