function [perPair, raw, nEff] = SEnergy(objs, s, refRange)
%SENERGY  Riesz s-energy diversity of a point set (lower = more uniform).
%
%   [perPair, raw, nEff] = SEnergy(objs, s, refRange)
%
%   INPUTS
%     objs     : (N x M) point set (a P^(t) or a Z^(t)).
%     s        : Riesz exponent alpha. Default M-1. (s<=0 -> 1.)
%     refRange : (2 x M) [lo; hi] used to normalize the space BEFORE computing
%                distances, so values are comparable across generations and
%                problems. Pass the reference-front range [min(PF);max(PF)].
%                If omitted/empty, the set's own range is used (NOT recommended
%                for {Z^(t)} drift analysis; kept only as a fallback).
%
%   OUTPUTS
%     perPair : E / nchoosek(nEff,2)   -- energy normalized by number of pairs.
%               Comparable across sets of different cardinality (adaptive
%               methods change |Z|). USE THIS for analysis.
%     raw     : E = sum_{i~=j} 1/||x_i - x_j||^s  (classic total energy).
%     nEff    : number of points AFTER de-duplication (diagnostic). Lets you
%               see how many near-coincident points were collapsed.
%
%   NUMERICS
%     The classic energy explodes when two points nearly coincide
%     (1/0^s). We (1) normalize the space by refRange, then (2) collapse
%     points closer than a relative tolerance TOL into one before computing
%     energy. This removes the spurious 1e10+ spikes seen in early
%     generations without inventing a magic distance floor.
%
%   E_s(X) = sum over unordered pairs, counted once, times 2 for the
%   symmetric double sum. Hardin & Saff (2004); minimizing E yields
%   asymptotically uniform distributions.

TOL = 1e-9;   % relative dedup tolerance in the normalized [0,1] space

X = objs;
if isempty(X)
    perPair = nan; raw = nan; nEff = 0; return;
end
[N, M] = size(X);
if nargin < 2 || isempty(s); s = M - 1; end
if s <= 0; s = 1; end

% ----- normalize space -----
if nargin >= 3 && ~isempty(refRange)
    lo = refRange(1,:);
    hi = refRange(2,:);
else
    lo = min(X,[],1);
    hi = max(X,[],1);
end
span = hi - lo;
span(span <= 0) = 1;                 % guard degenerate dimensions
Xn = (X - repmat(lo,N,1)) ./ repmat(span,N,1);

% ----- de-duplicate near-coincident points -----
Xn = dedup(Xn, TOL);
nEff = size(Xn,1);
if nEff < 2
    perPair = nan; raw = nan; return;  % diversity undefined for <2 points
end

% ----- energy -----
D   = pdist(Xn);                     % unique pairwise distances
raw = 2 * sum(1 ./ (D .^ s));        % factor 2 for the symmetric sum
perPair = raw / (nEff * (nEff - 1)); % == raw / (2*nchoosek(nEff,2))
end

% =========================================================================
function Y = dedup(X, tol)
%DEDUP  Collapse points within Euclidean distance < tol into a single point.
%   Greedy: keep a point, drop all later points within tol of it.
n = size(X,1);
keep = true(n,1);
for i = 1:n
    if ~keep(i); continue; end
    d = sqrt(sum((X(i+1:end,:) - repmat(X(i,:), n-i, 1)).^2, 2));
    close = false(n,1);
    close(i+1:end) = d < tol;
    keep(close) = false;
end
Y = X(keep,:);
end