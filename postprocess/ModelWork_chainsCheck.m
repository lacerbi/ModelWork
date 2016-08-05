function ModelWork_chainsCheck(mbag)

if isfield(mbag,'bag') && isfield(mbag,'project')
    mfit = mbag.bag;
elseif ~iscell(mbag)
    mfit = {mbag};
else
    mfit = mbag;
end

table = [1 1; 1 2; 1 3; 2 2; 2 3; ...
        2 3; 2 4; 2 4; 2 5; 2 5; ...
        3 4; 3 4; 3 5; 3 5; 3 5; ...
        4 5; 4 5; 4 5; 4 5; 4 5];

screensize = get(groot,'Screensize');

if isempty(mfit{1}.mp)
    temp = load('VestBMS_data.mat');
    mfit = ModelWork_loadFields('VestBMS',mfit,temp.data);
end

for i = 1:numel(mfit)
    figure;
    
    nchains = mfit{i}.sampling.nchains;
    nparams = size(mfit{i}.maptheta,2);
    nrows = table(nparams+1,1);
    ncols = table(nparams+1,2);
        
    for iParam = 1:nparams
        subplot(nrows,ncols,iParam);
        sumstats.R = mfit{i}.sampling.sumstats.R(iParam);
        sumstats.neff = mfit{i}.sampling.sumstats.neff(iParam);
        bounds = [mfit{i}.mp.bounds.LB(iParam),mfit{i}.mp.bounds.UB(iParam)];
        plotChains(nchains,mfit{i}.sampling.samples(:,iParam),mfit{i}.mp.params{iParam},bounds,sumstats);
    end
    
    % Add log likelihoods chain
    subplot(nrows,ncols,iParam+1);
    plotChains(nchains,mfit{i}.sampling.logliks,'LL',[],[]);
    
    set(gcf,'Color','w');
    set(gcf,'Position',screensize);
end

end

%--------------------------------------------------------------------------
function plotChains(nchains,samples,xstring,bounds,sumstats)

plotcol = [0 0 0];

nSamplesTot = size(samples,1);
n = floor(linspace(0,nSamplesTot,nchains+1));
if isempty(bounds)
    xmin = prctile(samples,0.1);
    xmax = prctile(samples,99.9);    
else
    xmin = bounds(1);
    xmax = bounds(2);
end

for i = 1:nchains
    X{i} = samples(n(i)+1:n(i+1));    
    [~,y{i},xmesh{i}] = kde(X{i},2^12,xmin,xmax);
    plot(xmesh{i},y{i},'Color',plotcol,'LineWidth',1); hold on;    
end

ybounds = ylim;
plotcols = [0.7 0.3 0.7; 0.3 0.9 0.3; 0.3 0.3 0.9];

for i = 1:nchains
    idx = round(linspace(1,size(X{1},1),503));
    yy = linspace(ybounds(1),ybounds(2),numel(idx));
    plot(X{i}(idx),yy,'Color',plotcols(i,:),'LineWidth',0.25); hold on;    
end

for i = 1:nchains
    plot(xmesh{i},y{i},'Color',plotcol,'LineWidth',1); hold on;    
end
if isempty(bounds)
    xbounds = xlim;
else
    xbounds = [xmin,xmax];
end

xstring(xstring == '_') = '-';
xlabel(xstring);

axis([xbounds,ybounds]);

if isempty(sumstats)
    [sumstats.R,sumstats.neff] = psrf(X{:});
end

title(['R = ' num2str(sumstats.R,'%.2f') ', neff = ' num2str(round(sumstats.neff))]);

set(gca,'YTick',[]);
set(gca,'TickLength',2*get(gca,'TickLength'),'TickDir','out');
box off;


end