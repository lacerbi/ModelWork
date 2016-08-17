function q = minqoverlap(X,Y)
%MINQOVERLAP Minimum quantile interval overlap metric.
%   Q = MINQOVERLAP(X,Y) returns the minimum quantile interval, symmetric
%   around the median, at which two empirical distributions X and Y begin 
%   to overlap.

X = sort(X);
Y = sort(Y);
nx = numel(X);
ny = numel(Y);
N = max(sqrt([nx,ny]));

qq = linspace(0,1,N);

o = overlap(qq);
idx = find(o,1);
if isempty(idx); q = Inf; return; end
if idx == 1; q = 0; return; end

qq = linspace(qq(idx-1),qq(idx),N);
o = overlap(qq);
idx = find(o,1);

q = qq(idx);

return;


function o = overlap(q)
    xx = X(round(1 + (nx-1)*(0.5 + q(:)/2*[-1 1])));
    yy = Y(round(1 + (ny-1)*(0.5 + q(:)/2*[-1 1])));
    
    o = ( xx(:,1) <= yy(:,2) & xx(:,2) >= yy(:,1) ) | ...
        ( yy(:,1) <= xx(:,2) & yy(:,2) >= xx(:,1) );
end

end

