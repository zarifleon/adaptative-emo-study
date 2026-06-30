function score = SolowPolasky(objs, theta)
%SolowPolasky  Solow-Polasky diversity measure (higher = more diverse).
%   score = SolowPolasky(objs, theta)
%     objs  : (N x M) point set
%     theta : scaling parameter of the correlation kernel. Default: 10.
%
%   Builds the (N x N) matrix  C_ij = exp(-theta * ||x_i - x_j||)  and returns
%   the Solow-Polasky measure  e^T C^{-1} e  (e = ones vector). Ranges in
%   [1, N]; 1 = all points identical, N = maximally dissimilar. Used as a
%   spread/diversity indicator for the PFA and for {Z^(t)} (Sec. 4, axis iii).
%
%   Reference: Solow & Polasky, "Measuring biological diversity",
%   Environmental and Ecological Statistics, 1994.
%
%   The kernel matrix can be ill-conditioned when points are very close; a
%   tiny ridge is added to C for numerical stability of the inverse.

X = objs;
N = size(X,1);
if N < 1; score = nan; return; end
if N == 1; score = 1; return; end
if nargin < 2 || isempty(theta); theta = 10; end

Dist = squareform(pdist(X));          % N x N
C    = exp(-theta .* Dist);
C    = C + 1e-10 * eye(N);            % ridge for invertibility
e    = ones(N,1);
score = e' * (C \ e);
end