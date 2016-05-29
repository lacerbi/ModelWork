function [data,datafile] = ModelWork_loaddatafile(project,data,type)
%MODELWORK_LOADDATAFILE Load main data file for a given project.

batchrunInitFun = str2func([project '_batchrunInit']);

% Load data structure if not provided
if isempty(data)
    temp = batchrunInitFun(data,type);
    if isfield(temp,'datafile')
        data = temp.datafile;
    end
end

% Standard data file name if not provided
if isempty(data); data = [project '_data.mat']; end

% Try loading file
if ischar(data)
    try
        [datapath,dataname,dataext] = fileparts(data);
        if isempty(dataext); dataext = '.mat'; end
        data = load(fullfile(datapath,[dataname,dataext]));
        datafile = [dataname,dataext];
    catch
        error(['Could not load data file ''' data ''' from disk.']);
    end
else
    datafile = [];
end

end