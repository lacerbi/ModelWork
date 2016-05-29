function Z = mylinkage(Y, method, pdistArg)
%LINKAGE Create hierarchical cluster tree.
%   Z = LINKAGE(X), where X is a matrix with two or more rows, creates a
%   matrix Z defining a tree of hierarchical clusters of the rows of X.
%   Clusters are based on the single linkage algorithm using Euclidean
%   distances between the rows of X. Rows of X correspond to observations
%   and columns to variables.
%
%   Z = LINKAGE(X,METHOD) creates a hierarchical cluster tree using the
%   the specified algorithm. The available methods are:
%
%      'single'    --- nearest distance (default)
%      'complete'  --- furthest distance
%      'average'   --- unweighted average distance (UPGMA) (also known as
%                      group average)
%      'weighted'  --- weighted average distance (WPGMA)
%      'centroid'  --- unweighted center of mass distance (UPGMC)
%      'median'    --- weighted center of mass distance (WPGMC)
%      'ward'      --- inner squared distance (min variance algorithm)
%
%   Z = LINKAGE(X,METHOD,METRIC) performs clustering based on the distance
%   metric METRIC between the rows of X. METRIC can be any of the distance
%   measures accepted by the PDIST function. The default is 'euclidean'.
%   For more information on PDIST and available distances, type HELP PDIST.
%   The centroid, median, and Ward's methods are intended only for the
%   Euclidean distance metric.
%
%   Z = LINKAGE(X, METHOD, PDIST_INPUTS) enables you to pass extra input
%   arguments to PDIST. PDIST_INPUTS should be a cell array containing
%   input arguments to be passed to PDIST.
%
%   Z = LINKAGE(Y) and Z = LINKAGE(Y,METHOD) are alternative syntaxes that
%   accept a vector representation Y of a distance matrix. Y may be a
%   distance matrix as computed by PDIST, or a more general dissimilarity
%   matrix conforming to the output format of PDIST.
%
%   The output matrix Z contains cluster information. Z has size m-1 by 3,
%   where m is the number of observations in the original data. Column 1
%   and 2 of Z contain cluster indices linked in pairs to form a binary
%   tree. The leaf nodes are numbered from 1 to m. They are the singleton
%   clusters from which all higher clusters are built. Each newly-formed
%   cluster, corresponding to Z(i,:), is assigned the index m+i, where m is
%   the total number of initial leaves. Z(i,1:2) contains the indices of
%   the two component clusters which form cluster m+i. There are m-1 higher
%   clusters which correspond to the interior nodes of the output
%   clustering tree. Z(i,3) contains the corresponding linkage distances
%   between the two clusters which are merged in Z(i,:), e.g. if there are
%   total of 30 initial nodes, and at step 12, cluster 5 and cluster 7 are
%   combined and their distance at this time is 1.5, then row 12 of Z will
%   be (5,7,1.5). The newly formed cluster will have an index 12+30=42. If
%   cluster 42 shows up in a latter row, that means this newly formed
%   cluster is being combined again into some bigger cluster.
%
%   The centroid and median methods can produce a cluster tree that is not
%   monotonic. This occurs when the distance from the union of two
%   clusters, r and s, to a third cluster is less than the distance between
%   r and s. In such a case, in a dendrogram drawn with the default
%   orientation, the path from a leaf to the root node takes some downward
%   steps. You may want to use another method when that happens.
%
%   You can provide the output Z to other functions including DENDROGRAM to
%   display the tree, CLUSTER to assign points to clusters, INCONSISTENT to
%   compute inconsistent measures, and COPHENET to compute the cophenetic
%   correlation coefficient.
%
%   Example: Compute four clusters of the Fisher iris data using Ward
%            linkage and ignoring species information, and see how the
%            cluster assignments correspond to the three species.
%
%       load fisheriris
%       Z = linkage(meas,'ward','euclidean');
%       c = cluster(Z,'maxclust',4);
%       crosstab(c,species)
%       dendrogram(Z)
%
%   See also PDIST, INCONSISTENT, COPHENET, DENDROGRAM, CLUSTER,
%   CLUSTERDATA, KMEANS, SILHOUETTE.

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.5 $

% Check for size and type of input
[k, n] = size(Y);
m = ceil(sqrt(2*n)); % m = (1+sqrt(1+8*n))/2, but works for large n
if k>1  % data matrix input
    if nargin<2
        method = 'single';
    end
    if nargin<3
        pdistArg = 'euclidean';
    end
    nargs = 3;
else % distance matrix input or bad input
    nargs = nargin;
end

if nargs==3 % should be data input
    if k == 1 && m*(m-1)/2 == n
        warning('stats:linkage:CallingPDIST',...
                ['You have used the syntax to call PDIST from within LINKAGE, but the first ',...
                 'input argument to LINKAGE appears to be a distance matrix already.']);
    end
    if k < 2
        error('stats:linkage:TooFewDistances',...
              'You have to have at least two observations to do a linkage.');
    end
    callPdist = true;
    if ischar(pdistArg)
        pdistArg = {pdistArg};
    elseif ~iscell(pdistArg)
        error('stats:linkage:BadPdistArgs',...
              'Third input must be a string or a cell array.');
    end
else % should be distance input
    callPdist = false;
    if n < 1
        error('stats:linkage:TooFewDistances',...
              'You have to have at least one distance to do a linkage.');
    end
    if k ~= 1 || m*(m-1)/2 ~= n
        error('stats:linkage:BadSize',...
              'The first input does not appear to be a distance matrix because its size is not compatible with the output of the PDIST function. A data matrix input must have more than one row.');
    end
end

% Selects appropriate method
methods = {'single',   'nearest'; ...
           'complete', 'farthest'; ...
           'average',  'upgma'; ...
           'weighted', 'wpgma'; ...
           'centroid', 'upgmc'; ...
           'median',   'wpgmc'; ...
           'ward''s',  'incremental'};
if nargs == 1 % set default switch to be 'si'
    s = 1;
else
    s = find(strncmpi(method,methods,length(method)));
    if isempty(s)
        error('stats:linkage:BadMethod','Unknown method name: %s.',method);
    elseif length(s)>1
        error('stats:linkage:BadMethod','Ambiguous method name: %s.',method);
    else
        if s>size(methods,1), s = s - size(methods,1); end
    end
end
methodStr = methods{s};
method = methodStr(1:2);


% The recursive distance updates for these three methods only make sense
% when the distance matrix contains Euclidean distances (which will be
% squared) or the distance metric is Euclidean
if  ~isempty(strmatch(method,['ce';'me';'wa']))
    if ~callPdist
        if (any(~isfinite(Y)) || ~iseuclidean(Y))
            warning('stats:linkage:NotEuclideanMatrix',...
                '%s linkage specified with non-Euclidean dissimilarity matrix.',methodStr);
        end
    else
        nonEuc = false;
        if (~isempty(pdistArg))
            if (~ischar (pdistArg{1}))
                nonEuc = true;
            else
                distMethods = {'euclidean'; 'minkowski';'mahalanobis'; };
                %pdistArg{1} represents the distance metric
                i = strmatch(lower(pdistArg{1}), distMethods);
                if length(i) > 1
                    error('stats:linkage:BadDistance',...
                        'Ambiguous ''DISTANCE'' argument:  %s.', pdistArg{1});
                elseif (isempty(i) || i == 3 || ...
                  (i == 2 && length(pdistArg) ~= 1 && isscalar(pdistArg{2}) && pdistArg{2} ~= 2) )
                    nonEuc = true;
                end
            end
            
        end
        if (nonEuc)
            warning('stats:linkage:NotEuclideanMethod',...
                '%s linkage specified with non-Euclidean distance metric.',methodStr);
        end
    end
end

if exist('linkagemex','file')==3
    % call mex file
    if callPdist
        Z = linkagemex(Y,method,pdistArg);
    else
        Z = linkagemex(Y,method);
    end
else
    %warning('stats:linkage:NoMexFilePresent',...
    %        '''mex'' file for linkage is not available, running ''m'' version.');
    if callPdist
        Y = mypdist(Y,pdistArg{:});
    end
    % optional old linkage function (use if mex file is not present)
    Z = mylinkageold(Y,method);
end

% Check if the tree is monotonic and warn if not.  Z is built so that the rows
% are in non-decreasing height order, thus we can look at the heights in order
% to determine monotonicity, rather than having to explicitly compare each parent
% with its children.
zdiff = diff(Z(:,3));
if any(zdiff<0)
    % With distances that are computed recursively (average, weighted, median,
    % centroid, ward's), errors can accumulate.  Two nodes that are really
    % at the same height in the tree may have had their heights calculated in
    % different ways, making them differ by +/- small amounts.  Make sure that
    % doesn't produce false non-monotonicity warnings.
    negLocs = find(zdiff<0);
    if any(abs(zdiff(negLocs)) > eps(Z(negLocs,3))) % eps(the larger of the two values)
        warning('stats:linkage:NonMonotonicTree',...
                'Non-monotonic cluster tree -- the %s linkage is probably not appropriate.',methodStr);
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%% OLD LINKAGE FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = mylinkageold(Y, method)
%LINKAGEOLD Create hierarchical cluster tree using only MATLAB code.

n = size(Y,2);
m = ceil(sqrt(2*n)); % (1+sqrt(1+8*n))/2, but works for large n
if isa(Y,'single')
   Z = zeros(m-1,3,'single'); % allocate the output matrix.
else
   Z = zeros(m-1,3); % allocate the output matrix.
end

% during updating clusters, cluster index is constantly changing, R is
% a index vector mapping the original index to the current (row, column)
% index in Y.  N denotes how many points are contained in each cluster.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;

% Square the distances so updates are easier.  The cluster heights will be
% square-rooted back to the original scale after everything is done.
if ~isempty(strmatch(method,['ce';'me';'wa']))
   Y = Y .* Y;
end

for s = 1:(n-1)
   if strcmp(method,'av')
      p = (m-1):-1:2;
      I = zeros(m*(m-1)/2,1);
      I(cumsum([1 p])) = 1;
      I = cumsum(I);
      J = ones(m*(m-1)/2,1);
      J(cumsum(p)+1) = 2-p;
      J(1)=2;
      J = cumsum(J);
      W = N(R(I)).*N(R(J));
      [v, k] = min(Y./W);
   else
      [v, k] = min(Y);
   end

   i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
   j = k - (i-1)*(m-i/2)+i;

   Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A

   % Update Y. In order to vectorize the computation, we need to compute
   % all the indices corresponding to cluster i and j in Y, denoted by I
   % and J.
   I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables
   U = [I1 I2 I3];
   I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
   J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

   switch method
   case 'si' % single linkage
      Y(I) = min(Y(I),Y(J));
   case 'co' % complete linkage
      Y(I) = max(Y(I),Y(J));
   case 'av' % average linkage
      Y(I) = Y(I) + Y(J);
   case 'we' % weighted average linkage
      Y(I) = (Y(I) + Y(J))/2;
   case 'ce' % centroid linkage
      K = N(R(i))+N(R(j));
      Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v)./K)./K;
   case 'me' % median linkage
      Y(I) = (Y(I) + Y(J))/2 - v /4;
   case 'wa' % Ward's linkage
      Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
	  N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
   end
   J = [J i*(m-(i+1)/2)-m+j];
   Y(J) = []; % no need for the cluster information about j.

   % update m, N, R
   m = m-1;
   N(n+s) = N(R(i)) + N(R(j));
   R(i) = n+s;
   R(j:(n-1))=R((j+1):n);
end

if ~isempty(strmatch(method,['ce';'me';'wa']))
   Z(:,3) = sqrt(Z(:,3));
end

Z(:,[1 2])=sort(Z(:,[1 2]),2);
end

%--------------------------------------------------------------------------
function Y = mypdist(X,dist,varargin)
%PDIST Pairwise distance between observations.
%   D = PDIST(X) returns a vector D containing the Euclidean distances
%   between each pair of observations in the M-by-N data matrix X. Rows of
%   X correspond to observations, columns correspond to variables. D is a
%   1-by-(M*(M-1)/2) row vector, corresponding to the M*(M-1)/2 pairs of
%   observations in X.
%
%   D = PDIST(X, DISTANCE) computes D using DISTANCE.  Choices are:
%
%       'euclidean'   - Euclidean distance (default)
%       'seuclidean'  - Standardized Euclidean distance. Each coordinate
%                       difference between rows in X is scaled by dividing
%                       by the corresponding element of the standard
%                       deviation S=NANSTD(X). To specify another value for
%                       S, use D=PDIST(X,'seuclidean',S).
%       'cityblock'   - City Block distance
%       'minkowski'   - Minkowski distance. The default exponent is 2. To
%                       specify a different exponent, use
%                       D = PDIST(X,'minkowski',P), where the exponent P is
%                       a scalar positive value.
%       'chebychev'   - Chebychev distance (maximum coordinate difference)
%       'mahalanobis' - Mahalanobis distance, using the sample covariance
%                       of X as computed by NANCOV. To compute the distance
%                       with a different covariance, use
%                       D =  PDIST(X,'mahalanobis',C), where the matrix C
%                       is symmetric and positive definite.
%       'cosine'      - One minus the cosine of the included angle
%                       between observations (treated as vectors)
%       'correlation' - One minus the sample linear correlation between
%                       observations (treated as sequences of values).
%       'spearman'    - One minus the sample Spearman's rank correlation
%                       between observations (treated as sequences of values).
%       'hamming'     - Hamming distance, percentage of coordinates
%                       that differ
%       'jaccard'     - One minus the Jaccard coefficient, the
%                       percentage of nonzero coordinates that differ
%       function      - A distance function specified using @, for
%                       example @DISTFUN.
%
%   A distance function must be of the form
%
%         function D2 = DISTFUN(XI, XJ),
%
%   taking as arguments a 1-by-N vector XI containing a single row of X, an
%   M2-by-N matrix XJ containing multiple rows of X, and returning an
%   M2-by-1 vector of distances D2, whose Jth element is the distance
%   between the observations XI and XJ(J,:).
%
%   The output D is arranged in the order of ((2,1),(3,1),..., (M,1),
%   (3,2),...(M,2),.....(M,M-1)), i.e. the lower left triangle of the full
%   M-by-M distance matrix in column order.  To get the distance between
%   the Ith and Jth observations (I < J), either use the formula
%   D((I-1)*(M-I/2)+J-I), or use the helper function Z = SQUAREFORM(D),
%   which returns an M-by-M square symmetric matrix, with the (I,J) entry
%   equal to distance between observation I and observation J.
%
%   Example:
%      % Compute the ordinary Euclidean distance
%      X = randn(100, 5);                 % some random points
%      D = pdist(X, 'euclidean');         % euclidean distance
%
%      % Compute the Euclidean distance with each coordinate difference
%      % scaled by the standard deviation
%      Dstd = pdist(X,'seuclidean');
%
%      % Use a function handle to compute a distance that weights each
%      % coordinate contribution differently
%      Wgts = [.1 .3 .3 .2 .1];           % coordinate weights
%      weuc = @(XI,XJ,W)(sqrt(bsxfun(@minus,XI,XJ).^2 * W'));
%      Dwgt = pdist(X, @(Xi,Xj) weuc(Xi,Xj,Wgts));
%
%   See also SQUAREFORM, LINKAGE, SILHOUETTE, PDIST2.

%   An example of distance for data with missing elements:
%
%      X = randn(100, 5);     % some random points
%      X(unidrnd(prod(size(X)),1,20)) = NaN; % scatter in some NaNs
%      D = pdist(X, @naneucdist);
%
%      function d = naneucdist(XI, XJ) % euclidean distance, ignoring NaNs
%      [m,p] = size(XJ);
%      sqdx = bsxfun(@minus,XI,XJ).^2;
%      pstar = sum(~isnan(sqdx),2); % correction for missing coords
%      pstar(pstar == 0) = NaN;
%      d = sqrt(nansum(sqdx,2) .* p ./ pstar);
%
%
%   For a large number of observations, it is sometimes faster to compute
%   the distances by looping over coordinates of the data (though the code
%   is more complicated):
%
%      function d = nanhamdist(XI, XJ) % hamming distance, ignoring NaNs
%      [m,p] = size(XJ);
%      nesum = zeros(m,1);
%      pstar = zeros(m,1);
%      for q = 1:p
%          notnan = ~(isnan((XI(q)) | isnan(XJ(:,q)));
%          nesum = nesum + (XI(q) ~= XJ(:,q)) & notnan;
%          pstar = pstar + notnan;
%      end
%      nesum(any() | nans((i+1):n)) = NaN;
%      d(k:(k+n-i-1)) = nesum ./ pstar;

%   Copyright 1993-2009 The MathWorks, Inc.
%   $Revision: 1.1.8.2 $ $Date: 2010/05/10 17:59:09 $

if nargin < 2
    dist = 'euc';
else
    if ischar(dist)
        methods = {'euclidean'; 'seuclidean'; 'cityblock'; 'chebychev'; ...
            'mahalanobis'; 'minkowski'; 'cosine'; 'correlation'; ...
            'spearman'; 'hamming'; 'jaccard'};
        i = strmatch(lower(dist), methods);
        if length(i) > 1
            error('stats:pdist:BadDistance',...
                'Ambiguous ''DISTANCE'' argument:  %s.', dist);
        elseif isempty(i)
            % Assume an unrecognized string is a user-supplied distance
            % function name, change it to a handle.
            distfun = str2func(dist);
            distargs = varargin;
            dist = 'usr';
        else
            dist = lower(methods{i}(1:3));
        end
    elseif isa(dist, 'function_handle') ||  isa(dist, 'inline')
        distfun = dist;
        distargs = varargin;
        dist = 'usr';
    else
        error('stats:pdist:BadDistance',...
            'The ''DISTANCE'' argument must be a string or a function.');
    end
end

% Integer/logical/char/anything data may be handled by a caller-defined
% distance function, otherwise it is converted to double.  Complex floating
% point data must also be handled by a caller-defined distance function.
if ~strcmp(dist,'usr')
    if ~isfloat(X)
        warning('stats:pdist:DataConversion', ...
            'Converting %s data to double.',class(X));
        X = double(X);
    elseif any(imag(X(:)))
        error('stats:pdist:InvalidData', ...
            'PDIST does not accept complex data for built-in distances.');
    end
end

[n,p] = size(X);

% Degenerate case, just return an empty of the proper size.
if n < 2
    if ~strcmp(dist,'usr')
        Y = zeros(1,0,class(X)); % X was single/double, or cast to double
    elseif isa(X,'single')
           Y = zeros(1,0,'single');
    else
           Y = zeros(1,0);
    end
    return;
end

additionalArg = [];
switch dist
    case 'seu' % Standardized Euclidean weights by coordinate variance
        if nargin < 3
            additionalArg =  nanvar(X);
            if any(additionalArg == 0)
                warning('stats:pdist:ConstantColumns',...
                    ['Some columns have zero standard deviation. ',...
                    'You may want to use other inverse weights or other distance metrics. ']);
            end
            additionalArg = 1./ additionalArg;
        else
            additionalArg = varargin{1};
            if ~(isvector(additionalArg) && length(additionalArg) == p...
                    && all(additionalArg >= 0))
                error('stats:pdist:InvalidWeights',...
                    ['The inverse weights for the standardized Euclidean metric must be a vector of ', ...
                    'non-negative values, with length equal to the number of columns in X.']);
            end
            if any(additionalArg == 0)
                  warning('stats:pdist:ZeroInverseWeights',...
                    ['Some columns of the inverse weight are zeros. ',...
                    'You may want to use other inverse weights or another distance metric. ']);
            end
            %We will apply the inverse weight on each coordinate in the sum
            %of squares.
            additionalArg = 1./ (additionalArg .^2);
        end
        
    case 'mah' % Mahalanobis
        if nargin < 3
            additionalArg = nancov(X);
           [T,flag] = chol(additionalArg);
        else
            additionalArg = varargin{1};
            if ~isequal(size(additionalArg),[p,p])
                error('stats:pdist:InvalidCov',...
                    ['The covariance matrix for the Mahalanobis metric must be a ',...
                    'square matrix with the same number of columns as X.']);
            end
            %use cholcov because we also need to check whether the matrix is symmetric
            [T,flag] = cholcov(additionalArg,0);
        end
            if flag ~= 0
                error('stats:pdist:InvalidCov',...
                    ['The covariance for matrix the Mahalanobis metric must be symmetric ',...
                    'and positive definite.']);
            end
        if ~issparse(X) 
             additionalArg = T \ eye(p);  %inv(T)
        end   
       
    case 'min' % Minkowski distance needs a third argument
        if nargin < 3  % use default value for exponent
            additionalArg = 2;
        elseif ( isscalar(varargin{1}) && varargin{1} > 0)
            additionalArg = varargin{1}; % get exponent from input args
        else
            error('stats:pdist:InvalidExponent',...
                'The exponent for the Minkowski metric must be a positive scalar.');
        end
    case 'cos' % Cosine
        [X,flag] = mynormalizeX(X);
        if flag
            warning('stats:pdist:zeroPoints',...
                ['Some points have small relative magnitudes, making them ', ...
                'effectively zero. ',...
                'Cosine metric may not be appropriate for these points.']);
        end
    case 'cor' % Correlation
        X = bsxfun(@minus,X,mean(X,2));
        [X, flag] = mynormalizeX(X);
        if flag
            warning('stats:pdist:constantPoints',...
                ['Some points have small relative standard deviations, making ', ...
                'them effectively constant. ',...
                'Correlation metric may not be appropriate for these points.']);
        end
    case 'spe'  %Spearman
        X = tiedrank(X')'; % treat rows as a series
        X = X - (p+1)/2; % subtract off the (constant) mean
        [X,flag] = mynormalizeX(X);
        if flag
            warning('stats:pdist:constantPoints',...
                ['Some points have too many ties, making ', ...
                'them effectively constant. ',...
                'Rank correlation metric may not be appropriate for these points.']);
        end
end

% Note that if there is any code for case 'che','euc' or 'cit' in the
% above switch dist block, the some code need to be repeated in the
% corresponding block below.
if strcmp(dist,'min') % Minkowski distance
    if isinf(additionalArg) %the exponent is inf
        dist = 'che';
        additionalArg = [];
    elseif additionalArg == 2 %the exponent is 2
        dist = 'euc';
        additionalArg = [];
    elseif additionalArg == 1 %the exponent is 1
        dist = 'cit';
        additionalArg = [];
    end
end

% Call a mex file to compute distances for the standard distance measures
% and full real double or single data.
if ~strcmp(dist,'usr') && (isfloat(X) && ~issparse(X)) ...
        && exist('pdistmex','file')==3 % ~usr => ~complex
    additionalArg = cast(additionalArg,class(X));
    Y = pdistmex(X',dist,additionalArg);
    
    % This M equivalent assumes real single or double.  It is currently only
    % called for sparse inputs, but it may also be useful as a template for
    % customization.
elseif ~strcmp(dist,'usr') && isfloat(X) % ~usr => ~complex
    if strmatch(dist, {'ham' 'jac' 'che'})
        nans = any(isnan(X),2);
    end
    outClass = class(X);
    Y = zeros(1,n*(n-1)./2, outClass);
    k = 1;
    for i = 1:n-1
        switch dist
            case 'euc'    % Euclidean
                dsq = zeros(n-i,1,outClass);
                for q = 1:p
                    dsq = dsq + (X(i,q) - X((i+1):n,q)).^2;
                end
                Y(k:(k+n-i-1)) = sqrt(dsq);
                
            case 'seu'    % Standardized Euclidean
                wgts = additionalArg;
                dsq = zeros(n-i,1,outClass);
                for q = 1:p
                    dsq = dsq + wgts(q) .* (X(i,q) - X((i+1):n,q)).^2;
                end
                Y(k:(k+n-i-1)) = sqrt(dsq);
                
            case 'cit'    % City Block
                d = zeros(n-i,1,outClass);
                for q = 1:p
                    d = d + abs(X(i,q) - X((i+1):n,q));
                end
                Y(k:(k+n-i-1)) = d;
                
            case 'mah'    % Mahalanobis
     
                 del = bsxfun(@minus, X(i,:), X((i+1):n,:));
                 dsq = sum((del/T) .^ 2, 2);
                 Y(k:(k+n-i-1)) = sqrt(dsq);
            case 'min'    % Minkowski
                expon = additionalArg;
                dpow = zeros(n-i,1,outClass);
                for q = 1:p
                    dpow = dpow + abs(X(i,q) - X((i+1):n,q)).^expon;
                end
                Y(k:(k+n-i-1)) = dpow .^ (1./expon);
                
            case {'cos' 'cor' 'spe'}   % Cosine, Correlation, Rank Correlation
                % This assumes that data have been appropriately preprocessed
                d = zeros(n-i,1,outClass);
                for q = 1:p
                    d = d + (X(i,q).*X((i+1):n,q));
                end
                d(d>1) = 1; % protect against round-off, don't overwrite NaNs
                Y(k:(k+n-i-1)) = 1 - d;
                
            case 'ham'    % Hamming
                nesum = zeros(n-i,1,outClass);
                for q = 1:p
                    nesum = nesum + (X(i,q) ~= X((i+1):n,q));
                end
                nesum(nans(i) | nans((i+1):n)) = NaN;
                Y(k:(k+n-i-1)) = nesum ./ p;
                
            case 'jac'    % Jaccard
                nzsum = zeros(n-i,1,outClass);
                nesum = zeros(n-i,1,outClass);
                for q = 1:p
                    nz = (X(i,q) ~= 0 | X((i+1):n,q) ~= 0);
                    ne = (X(i,q) ~= X((i+1):n,q));
                    nzsum = nzsum + nz;
                    nesum = nesum + (nz & ne);
                end
                nesum(nans(i) | nans((i+1):n)) = NaN;
                Y(k:(k+n-i-1)) = nesum ./ nzsum;
                
            case 'che'    % Chebychev
                dmax = zeros(n-i,1,outClass);
                for q = 1:p
                    dmax = max(dmax, abs(X(i,q) - X((i+1):n,q)));
                end
                dmax(nans(i) | nans((i+1):n)) = NaN;
                Y(k:(k+n-i-1)) = dmax;
                
        end
        k = k + (n-i);
    end
    
    % Compute distances for a caller-defined distance function.
else % if strcmp(dist,'usr')
    try
        Y = feval(distfun,X(1,:),X(2,:),distargs{:})';
    catch ME
        if strcmp('MATLAB:UndefinedFunction', ME.identifier) ...
                && ~isempty(strfind(ME.message, func2str(distfun)))
            error('stats:pdist:DistanceFunctionNotFound',...
                'The distance function ''%s'' was not found.', func2str(distfun));
        end
        % Otherwise, let the catch block below generate the error message
        Y = [];
    end
    
    % Make the return have whichever numeric type the distance function
    % returns, or logical.
    if islogical(Y)
        Y = false(1,n*(n-1)./2);
    else % isnumeric
        Y = zeros(1,n*(n-1)./2, class(Y));
    end
    
    k = 1;
    for i = 1:n-1
        try
            Y(k:(k+n-i-1)) = feval(distfun,X(i,:),X((i+1):n,:),distargs{:})';
        catch ME
            if isa(distfun, 'inline')
                throw(addCause(MException('stats:pdist:DistanceFunctionError',...
                    'Error evaluating inline distance function.'),...
                    ME));
            else
                throw(addCause(MException('stats:pdist:DistanceFunctionError',...
                    'Error evaluating distance function ''%s''.',...
                    func2str(distfun)),...
                    ME));
            end
        end
        k = k + (n-i);
    end
end
end

%---------------------------------------------
% Normalize the rows in X to have unit norm.
function [X,flag] = mynormalizeX(X)
% Rescale each row by largest element to prevent over/underflow, and
% compute the 2-norm.
Xmax = max(abs(X),[],2);
X2 = bsxfun(@rdivide,X,Xmax);
Xnorm = sqrt(sum(X2.^2, 2));

% The norm will be NaN for rows that are all zeros, fix that for the test
% below.
Xnorm(Xmax==0) = 0;

% The norm will be NaN for rows of X that have any +/-Inf. Those should be
% Inf, but leave them as is so those rows will not affect the test below.
% The points can't be normalized, so any distances from them will be NaN
% anyway.

% Find points that are effectively zero relative to the point with largest norm.
flag = any(Xnorm <= eps(max(Xnorm)));

% Back out the rescaling, and normalize rows of X to have unit 2-norm.
% Rows can't be normalized become all NaN.
Xnorm = Xnorm .* Xmax;
X = bsxfun(@rdivide,X,Xnorm);
end
