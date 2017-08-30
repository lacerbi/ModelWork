function gendata = ModelWork_gendata(project,N,varargin)
%MODELWORK_GENDATA Generate a number of fake datasets.
%   GENDATA = MODELWORK_GENDATA(PROJECT,N,MFIT) generates N fake datasets 
%   from individual model struct MFIT. GENDATA is a cell array of generated 
%   data matrices.
%
%   GENDATA = MODELWORK_GENDATA(PROJECT,N,X,MP) generates fake datasets from 
%   data matrix X and model parameter structure MP.
%
%  By Luigi Acerbi <luigi.acerbi@gmail.com>

if nargin < 5
    mfit = ModelWork_loadFields(project,varargin{1});
    X = mfit.X;
    mp = mfit.mp;
    infostruct = mfit.infostruct;
else
    X = varargin{1};
    mp = varargin{2};
    infostruct = varargin{3};
end

if ischar(varargin{end}) && strncmpi(varargin{end},'r',1)
    regenerate = 1;
else
    regenerate = 0;
end

setupModelFun = str2func([project '_setupModel']);
gendataFun = str2func([project '_gendata']);

gendata = [];
if ~isempty(mfit.sampling) && ~isempty(mfit.sampling.samples)
    Nsamples = size(mfit.sampling.samples,1);
    idx = round(linspace(1,Nsamples,N));
    trueparams = mfit.sampling.samples(idx,:);
    for i = 1:N
        mp = setupModelFun(mp, trueparams(i,:), mfit.model, infostruct);
        XX = gendataFun(1,mp,X,infostruct,regenerate);
        gendata{i} = squeeze(XX(:, :, 1));
    end
else
    trueparams = repmat(mfit.maptheta,[N 1]);
    mp = setupModelFun(mp, mfit.maptheta, mfit.model, infostruct);    
    XX = gendataFun(N,mp,X,infostruct,regenerate);    
    for i = 1:N; gendata{i} = squeeze(XX(:, :, i)); end    
end

end