function ModelWork_table(modelsummary,modellist,metric,notes,printnames)
% MODELWORK_TABLE LaTeX table of model comparison scores for given metric.

if nargin < 4 || isempty(notes); notes = 'Table notes.'; end
if nargin < 5; printnames = []; end

rankflag = 1;

% Convert from model names to model vectors
if iscell(modellist) && ischar(modellist{1})
    temp = [];
    for k = 1:numel(modellist)
        idx = find(strcmp(modellist{k},modelsummary.modelnames));        
        if isempty(idx) || numel(idx) > 1
            error(['Cannot find a unique model in MODELSUMMARY matching model #' num2str(k) '.']);
        end
        temp(k,:) = modelsummary.models(idx,:);
    end
    modellist = temp;
end

% Collect model comparison metric
tab = [];
getnamesflag = isempty(printnames);
for k = 1:size(modellist,1)
    % Find model in full model list
    model = modellist(k,:);
    idx = find(all(bsxfun(@eq, modellist(k,:), modelsummary.models),2));
    if isempty(idx) || numel(idx) > 1
        error(['Cannot find a unique model in MODELSUMMARY matching model #' num2str(k) '.']);
    end    
    tab(k,:) = modelsummary.(metric)(:,idx)';
    
    % Printed model names (might differ from stored model names)
    if getnamesflag
        if isfield(modelsummary, 'modelnames') && ~isempty(modelsummary.modelnames)
            printnames{k} = modelsummary.modelnames{idx};
        else
            printnames{k} = ['M' num2str(idx)];
        end
    end
end

[~,bestidx] = max(mean(tab,2));
tab = bsxfun(@minus,tab,tab(bestidx,:));

% Order model from best to worst
if rankflag
   meantab = mean(tab,2);
   [~,ord] = sort(meantab,'descend');   
   tab = tab(ord,:);
   printnames = printnames(ord);
   bestidx = 1;
end    
    
% Print LaTeX table

fprintf('\n{\n\\begin{center}\n');
% fprintf('\\begin{adjustwidth}{-2.25in}{0in} %% Comment out/remove adjustwidth environment if table fits in text column.\n');
fprintf('\\centering\n');
fprintf('%% \\caption{{\\bf Table caption}}\n');
fprintf('{ %% \\scriptsize %% Adjust table font size\n');
fprintf('\\begin{tabular}{l|');
for i = 1:size(tab,2); fprintf('c'); end
fprintf('|c}\n');


fprintf('{\\bf Model} ')
for i = 1:size(tab,2); fprintf(' & {\\bf S%d}', i); end
fprintf(' & {\\bf Mean} $\\bm{\\pm}$ {\\bf SE}\\\\\n\\hline\n');

for k = 1:size(tab,1)    
    if k == bestidx
        fprintf(' {\\bf %s} ',printnames{k});
        for i = 1:size(tab,2); fprintf(' & $\\bm{%.1f}$ ', tab(k,i)); end
        fprintf(' & $\\bm{%.1f \\pm %.1f}$ ', mean(tab(k,:)), stderr(tab(k,:)));        
    else
        fprintf(' %s ',printnames{k});
        for i = 1:size(tab,2); fprintf(' & $%.1f$ ', tab(k,i)); end
        fprintf(' & $%.1f \\pm %.1f$ ', mean(tab(k,:)), stderr(tab(k,:)));
    end
    fprintf('\\\\\n');
end

fprintf('\\end{tabular}\n}\n')
fprintf('\\begin{flushleft} %s \\end{flushleft}\n', notes);
% fprintf('\\end{adjustwidth} %% Comment out/remove adjustwidth environment if table fits in text column.\n');
fprintf('\\end{center}\n}\n\n');


% \hline
% \caption{{\bf Observer models by task.}}
% \textbf{Observer model} & \textbf{Parameters} & \textbf{\#} \\
% \hline
% \begin{flushleft} Table notes. \end{flushleft}



end