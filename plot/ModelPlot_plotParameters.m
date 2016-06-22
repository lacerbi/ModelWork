function ModelPlot_plotParameters(project,varargin)
%MODELPLOT_PLOTPARAMETERS Plot model parameters.

mfit = [];
for i = 1:numel(varargin)
    % Check that input is a model fit structure
    if isfield(varargin{i},'X') && isfield(varargin{i},'mp') && isfield(varargin{i},'metrics')
        mfit{end+1} = varargin{i};
    elseif iscell(varargin{1}) && ...
            isfield(varargin{1}{1},'X') && isfield(varargin{1}{1},'mp') && isfield(varargin{1}{1},'metrics')
        mfit = varargin{1};
    elseif ~isempty(varargin{i})
        error('Input is not a model fit structure.');
    end
end

% Remove empty fits
mfitnew = [];
for i = 1:numel(mfit); if ~isempty(mfit{i}); mfitnew{end+1} = mfit{i}; end; end
mfit = mfitnew; clear mfitnew;

% Loop over models, get all parameters
mfit = ModelWork_loadFields(project,mfit);

params = [];
for i = 1:numel(mfit)
    newp = mfit{i}.mp.params;
    for j = 1:numel(newp)
        if ~any(strcmp(newp{j},params)); params{end+1} = newp{j}; end
    end
end

% Define subplot grid
gridsize = [1 2; 1 3; 2 2; 2 3; 2 3; 2 4; 2 4; 3 3; 3 4; 3 4; 3 4; ...
    4 4; 4 4; 4 4; 4 4; 4 5; 4 5; 4 5; 4 5];
colors = [  ...
    141 211 199; ...          
    251 128 114; ...            
    128 177 211; ...
    253 180 98; ...
    160 120 100; ...
    70 70 233; ...
    252 205 229; ...
    165 211 195; ...
    159 212 105; ...
    188 128 189; ...
    212 148 169; ...      
    0 0 0 ...
    ]/255;


nparams = numel(params);

for j = 1:nparams
    subplot(gridsize(nparams,1),gridsize(nparams,2),j);
    for i = 1:numel(mfit)
        idx = find(strcmp(params{j},mfit{i}.mp.params),1);
        if isempty(idx); continue; end
        samples = mfit{i}.sampling.samples(:,idx);
        bounds = [mfit{i}.mp.bounds.LB(idx),mfit{i}.mp.bounds.UB(idx)];
        
        % Plot histogram of posterior over parameter
        nbins = 40;
        xx = linspace(min(samples),max(samples),nbins);
        dx = xx(2)-xx(1);
        pdf = histc(samples,xx);
        xx = [xx(1)-sqrt(eps), xx, xx(end)+sqrt(eps)];
        pdf = [0; pdf/sum(pdf)/dx; 0];

        xlim(bounds);
        plot(xx,pdf,'-','Color',colors(i,:),'LineWidth',1);
        hold on;        
    end
    
    textstring = params{j};
    textstring(textstring == '_') = '-'; 
    xlabel(textstring);
    ylabel('Probability density');
    
    set(gca,'TickDir','out');
    box off;
    
end

set(gcf,'Color','w');

% Write legend
subplot(gridsize(nparams,1),gridsize(nparams,2),prod(gridsize(nparams,:)));
for i = 1:numel(mfit)
    plot([0 0], [0 0], '-', 'Color',colors(i,:),'LineWidth',1); hold on;
    names{i} = VestBMS_getModelName(mfit{i}.model);
end

hl = legend(names{:});
legend('boxoff');
set(hl,'Location','NorthWest');
box off;
axis off;




