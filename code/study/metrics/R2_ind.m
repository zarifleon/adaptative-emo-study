function score = R2_ind(objs, W, z)
%R2_ind  Unary R2 indicator (utility-based, weighted Tchebycheff).
%   score = R2_ind(objs, W, z)
%     objs : (N x M) approximation set (minimization)
%     W    : (L x M) set of weight vectors (direction set lambda).
%            Default: UniformPoint(100, M) if omitted/empty.
%     z    : (1 x M) ideal/reference point. Default: min over objs per dim.
%
%   R2(A) = (1/|W|) * sum_{w in W} min_{a in A} g_tch(a | w, z)
%   with augmented-free weighted Tchebycheff  g = max_j w_j * |a_j - z_j|.
%   Smaller is better. Weakly Pareto compliant -> valid for monotonicity.
%
%   NOTE: weights with a zero component are handled by replacing 0 with a
%   tiny epsilon (standard practice) so the Tchebycheff product is defined.
%
%   Reference: Brockhoff, Wagner, Trautmann, "On the properties of the R2
%   indicator", GECCO 2012.

A = objs;
if isempty(A); score = nan; return; end
M = size(A,2);

if nargin < 2 || isempty(W)
    W = UniformPoint(100, M);
end
if nargin < 3 || isempty(z)
    z = min(A, [], 1);
end

W(W==0) = 1e-6;                 % avoid zero weights
L = size(W,1);
Na = size(A,1);

util = zeros(1,L);
for k = 1 : L
    w  = W(k,:);
    % weighted Tchebycheff of every a w.r.t. (w,z): max_j w_j*|a_j - z_j|
    g  = max( repmat(w, Na, 1) .* abs(A - repmat(z, Na, 1)), [], 2 );
    util(k) = min(g);
end
score = mean(util);
end