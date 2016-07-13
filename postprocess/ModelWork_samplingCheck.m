function exitflag = ModelWork_samplingCheck(mbag,method,type,nreplicas)
% MODELWORK_SAMPLINGCHECK check sampling results.
%

if nargin < 2 || isempty(method); method = 'plot'; end
if nargin < 3; type = []; end
if nargin < 4 || isempty(nreplicas); nreplicas = 3; end

switch lower(method(1))
    case 'p'
        info.plotflag = 1;
        info.writeflag = 0;
    case 'w'
        info.plotflag = 0; 
        info.writeflag = 1;
    otherwise
        error('Unknown METHOD. METHOD can be ''plot'' or ''write''.');
end

if info.plotflag; figure; end

compute_stats = 1;
if isfield(mbag,'bag') && isfield(mbag,'table')
    mfit = mbag.bag;
elseif ~iscell(mbag)
    mfit = {mbag};
else
    mfit = mbag;
end

n = numel(mfit);
info.project = mbag.project;
info.nrows = ceil(sqrt(n+1));
info.ncols = ceil((n+1)/info.nrows);
info.n = n;

exitflag = zeros(1,n);
g = zeros(info.nrows,info.ncols);
for i = 1:n
    icol = mod((i-1),info.ncols)+1;
    irow = 1 + floor((i-1)/info.ncols);
    if info.plotflag; tight_subplot(info.nrows,info.ncols,irow,icol); end
    [exitflag(i)] = samplingCheck(mfit{i},info);
    if info.plotflag
        drawnow;
        if icol ~= 1; set(gca,'YtickLabel',[]); end
    end
    g(irow,icol) = 1;
end

if any(exitflag)
    fprintf('%d datasets with major issues and %d datasets with minor issues.\n', numel(find(exitflag == 2)), numel(find(exitflag == 1)));
else
    fprintf('All datasets and model present no issues.\n');
end

if info.writeflag
    if isempty(type); error('TYPE not specified.'); end
    iProc = 1;
    for i = find(exitflag)
        model = packuint(mfit{i}.model);
        dataid = packuint(mfit{i}.dataid);
        cnd = mfit{i}.cnd;
        for rep = 1:nreplicas
            jobfilename = [num2str(iProc) '.job'];
            fout = fopen(jobfilename,'w');
            fprintf(fout,'%d %s %s %s %s\n',type,model,dataid,numarray2str(cnd),numarray2str(rep));
            fclose(fout);
            iProc = iProc + 1;
        end
    end    
end

end

%--------------------------------------------------------------------------
function [exitflag,minstats] = samplingCheck(mfit,info)

minstats = [];

thresh = [100,300];

% Read diagnostics
R = mfit.sampling.sumstats.R;
neff = mfit.sampling.sumstats.neff;

if max(R) > 1.1 || min(neff) < thresh(1)
    exitflag = 2;
    bkgcol = [1 0.7 0.7];
    % set(gca,'Color',[1 0.7 0.7]);
    plotcol = [1 0 0];
elseif max(R) > 1.05  || min(neff) < thresh(2)
    exitflag = 1;
    bkgcol = [1 1 0.9];
    % set(gca,'Color',[1 1 0.9]);        
    plotcol = [0.7 0.7 0];
else
    exitflag = 0;
    bkgcol = [1 1 1];
    plotcol = [0 0 0];
end

if ~info.plotflag; return; end
    
% Take log likelihood chains
logliks = sum(mfit.sampling.logliks,2);
nchains = mfit.sampling.nchains;
nSamplesTot = size(logliks,1);
n = floor(linspace(0,nSamplesTot,nchains+1));
xmin = prctile(logliks,0.25);
xmax = prctile(logliks,99.75);

for i = 1:nchains
    ll{i} = logliks(n(i)+1:n(i+1),:);    
    [~,y,xmesh] = kde(ll{i},2^12,xmin,xmax);
    plot(xmesh,y,'Color',plotcol,'LineWidth',1); hold on;    
end


box off;
set(gca,'TickDir','out');
xticks = [];
yticks = [];
xticklabel = [];
yticklabel = [];
set(gca,'XTick',xticks,'XTickLabel',xticklabel,'YTick',yticks,'YTickLabel',yticklabel);
set(gca,'TickLength',6*get(gca,'TickLength'));
set(gca,'ButtonDownFcn',@(src,evt) axesCallback(src,evt,mfit,info));
% axis([xx(1),xx(end), ymin ymax])
% axis off;
set(gca,'XColor','none','YColor','none');

t1 = []; t2 = [];
if max(R) > 1.05; t1 = ['R = ' num2str(max(R),'%.2f')]; end
if min(neff) < thresh(2); t2 = ['neff = ' num2str(round(min(neff)),'%g')]; end

if ~isempty(t1); text(0,0.8,t1,'Units','Normalized'); end
if ~isempty(t2); text(0,0.5,t2,'Units','Normalized'); end
set(gcf,'Color','w');

end

%--------------------------------------------------------------------------
function axesCallback(src,evt,mfit,info)

n = info.n;

for i = 1:n
    icol = mod((i-1),info.ncols)+1;
    irow = 1 + floor((i-1)/info.ncols);
    tight_subplot(info.nrows,info.ncols,irow,icol);
    set(gca,'Color',[1 1 1]);
end

icol = mod(n,info.ncols)+1;
irow = 1 + floor(n/info.ncols);
tight_subplot(info.nrows,info.ncols,irow,icol);
set(src,'Color',[1 1 0.85]);

cla;

getModelName = str2func([info.project '_getModelName']);

s1 = ['Dataset: ' numarray2str(mfit.dataid,[],',','','')];
% s2 = ['Model: ' numarray2str(mfit.model,[],[],'','')];
s2 = ['Model: ' getModelName(mfit.model)];
s3 = ['Cnd: ' numarray2str(mfit.cnd,[],[],'','')];
text(0.1,0.6,s1,'Units','Normalized');
text(0.1,0.35,s2,'Units','Normalized');
text(0.1,0.1,s3,'Units','Normalized');

axis off;

end
