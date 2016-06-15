function [exitflag,minstats] = ModelWork_optimizationCheck(mbag,thresh)
% MODELWORK_OPTIMIZATIONCHECK check optimization results.
%
% MODELSUMMARY = MODELWORK_SUMMARY(MBAG) merges batches with identifiers in 
% list ID returning the merged model bag MBAG.
%

if nargin < 2 || isempty(thresh); thresh = 0.1; end

figure;

compute_stats = 1;
if isfield(mbag,'ferr_median') && isfield(mbag,'minval')
    minstats = mbag;
    compute_stats = 0;
elseif isfield(mbag,'bag') && isfield(mbag,'table')
    mfit = mbag.bag;
elseif ~iscell(mbag)
    mfit = {mbag};
else
    mfit = mbag;
end

clear mbag;

%% Compute bootstrapped minima

if compute_stats; n = numel(mfit); else n = numel(minstats); end
nrows = ceil(sqrt(n+1));
ncols = ceil((n+1)/nrows);

%nrows = 8;
%ncols = 10;
exitflag = zeros(1,n);
g = zeros(nrows,ncols);
for i = 1:n
    icol = mod((i-1),ncols)+1;
    irow = 1 + floor((i-1)/ncols);
    subplot(nrows,ncols,i);
    if compute_stats
        [exitflag(i),minstats(i)] = optimizationCheck(mfit{i},thresh);
        drawnow;
    else
        [exitflag(i),minstats(i)] = optimizationCheck(minstats(i),thresh);
    end
    if icol ~= 1; set(gca,'YtickLabel',[]); end
    g(irow,icol) = 1;
end

icol = 1; irow = ceil(nrows/2);
subplot(nrows,ncols,(irow-1)*ncols + icol);
ylabel('Estimated error on the minimum');
if ~g(irow,icol); axis off; end

icol = ceil(ncols/2); irow = nrows;
subplot(nrows,ncols,(irow-1)*ncols + icol);
if g(irow,icol)
    xlabel('Bootstrapped restarts');
else
    axis([0 1 0 1]);
    text(0.5,0.5,'Bootstrapped restarts','Units','Normalized','HorizontalAlignment','Center');
    axis off;
end

subplot(nrows,ncols,nrows*ncols);
plot([0 0],[0 0],'k:','LineWidth',1); hold on;
plot([0 0],[0 0],'k-','LineWidth',1);    
plot([0 0],[0 0],'k--','LineWidth',1);
plot([0 0],[0 0], '-', 'LineWidth',2,'Color',0.6*[1 1 1]);

hl = legend('Min 95% UCB','Min median','Error threshold','Reached threshold');
set(hl,'Box','off');
axis([0 1 0 1]);
axis off;

%% Plot summary statitics of bootstrapped minima

figure;
subplot(2,1,1);

t = [];
for i = 1:n
    % Threshold t
    idx = find(minstats(i).ferr_ucb < thresh,1);
    if isempty(idx)
        t = [t, Inf]; 
    else
        t = [t, minstats(i).nrestarts(idx)/minstats(i).nrestarts(end)]; 
    end
end

dt = 0.05;
t = round(t/dt)*dt;

tvals = sort(unique(t(isfinite(t))));

bincounts = histc(t,[0,tvals(1:end-1) + 0.5*diff(tvals),tvals(end)+1]);
bincounts = [bincounts(1:end-1), sum(isinf(t))];

x = [tvals, 1.1];

bar(x,bincounts,'FaceColor',0.7*[1 1 1],'EdgeColor','none');
for i = 1:numel(x)
    text(x(i),bincounts(i)+0.05*max(bincounts),[num2str(bincounts(i))],'HorizontalAlignment','Center');
end

xticks = x;
for i = 1:numel(xticks)-1; xticklabels{i} = num2str(xticks(i)); end 
xticklabels{numel(xticks)} = 'Inf';

set(gca,'Xtick',xticks,'XTickLabel',xticklabels,'TickDir','out');
box off;
set(gcf,'Color','w');
xlabel('Fraction restarts before hitting threshold with 95% confidence');
ylabel('# optimization problems');


subplot(2,1,2);

med = []; ucb = [];
for i = 1:n
    med = [med,minstats(i).ferr_median(end)];
    ucb = [ucb,minstats(i).ferr_ucb(end)];
end
TolFun = 0.01;
med = log(max(med,TolFun));
ucb = log(max(ucb,TolFun));

plot(1:n,ucb,'k:','LineWidth',1); hold on;
plot(1:n,med,'k-','LineWidth',1);
plot([1 n],log(thresh)*[1,1],'k--','LineWidth',1)

xlabel('Optimization problems');
ylabel('Estimated error');
set(gca,'TickDir','out');
box off;

ymin = TolFun;
ymax = (thresh^2/TolFun);
xticks = [1,n];
yticks = [ymin,thresh,ymax];
for i = 1:numel(yticks); yticklabel{i} = num2str(yticks(i)); end
set(gca,'Xtick',xticks,'Ytick',log(yticks),'YtickLabel',yticklabel);

hl = legend('Min 95% UCB','Min median','Error threshold');
set(hl,'Box','off');

end

%--------------------------------------------------------------------------
function [exitflag,minstats] = optimizationCheck(mfit,thresh)

TolFun = 0.01;  % Maximum precision
logflag = 1;    % Plot in log coordinates
nboot = 1e4;    % Bootstrap samples

if isfield(mfit,'minval') && isfield(mfit,'ferr_median')
    minstats = mfit;
    n = numel(minstats.ferrs);
else
    % Read optimization results
    fvals = [];
    data = mfit.optimization.output;
    for i = 1:numel(data); fvals = [fvals; data{i}.fval(:)]; end

    n = numel(fvals);

    % Prepare min stats table
    frac = 0.75;

    minstats.minval = min(fvals);
    minstats.ferrs = fvals - minstats.minval;
    minstats.nrestarts = unique(sort(round(n*frac.^(0:8))));

    for i = 1:numel(minstats.nrestarts)
        idx = randi(n,[minstats.nrestarts(i),nboot]);
        m = min(minstats.ferrs(idx),[],1);
        minstats.ferr_median(i) = median(m);
        minstats.ferr_ucb(i) = prctile(m,95);
    end

end

clear mfit;

ferr_median = max(minstats.ferr_median,TolFun);
ferr_ucb = max(minstats.ferr_ucb,TolFun);

% xticks = [10,30,100,1e3,1e4];
ymin = TolFun;
ymax = (thresh^2/TolFun);
xticks = minstats.nrestarts;
for i = [1,ceil(numel(xticks)/2),numel(xticks)]; xticklabel{i} = num2str(xticks(i)); end
yticks = [ymin,thresh,ymax];
for i = 1:numel(yticks); yticklabel{i} = num2str(yticks(i)); end
threshi = thresh*10;
xthresh = minstats.nrestarts(find(minstats.ferr_ucb < thresh, 1));
if isempty(xthresh); xthresh = Inf; end

xx = minstats.nrestarts;

if logflag
    ymin = log(ymin);
    ymax = log(ymax);
    xticks = log(xticks);
    yticks = log(yticks);
    xx = log(xx);
    ferr_median = log(ferr_median);
    ferr_ucb = log(ferr_ucb);
    thresh = log(thresh);
    threshi = log(threshi);
    xthresh = log(xthresh);    
end

% Plot graph
plot(xx,ferr_ucb,'k:','LineWidth',1); hold on;
plot(xx,ferr_median,'k-','LineWidth',1);    
plot([xx(1),xx(end)], thresh*[1 1],'k--','LineWidth',1);
box off;
set(gca,'TickDir','out','XTick',xticks,'XTickLabel',xticklabel,'YTick',yticks,'YTickLabel',yticklabel);
set(gca,'TickLength',6*get(gca,'TickLength'));
axis([xx(1),xx(end), ymin ymax])

if isfinite(xthresh); plot(xthresh*[1 1], [ymin ymax], '-', 'LineWidth',2,'Color',0.6*[1 1 1]); end

if ferr_ucb(end) > threshi || ferr_median(end) > thresh
    exitflag = 2;
    set(gca,'Color',[1 0.7 0.7]);
elseif ferr_ucb(end) > thresh  || xthresh == xx(end)
    exitflag = 1;
    set(gca,'Color',[1 1 0.9]);
else
    exitflag = 0;
end

set(gcf,'Color','w');

end