function score = IGD_mat(objs, optimum)
%IGD_mat  Inverted Generational Distance, matrix interface, optimized.
%   Identical to PlatEMO metric/IGD.m. NOT Pareto compliant; reported for
%   comparability with prior work only (do NOT use for monotonicity).
%
%   IGD = mean over reference points r of  min over a in A of ||a - r||.
%   Loops over the N approximation points (small) instead of Nr reference
%   points (large); accumulates the per-reference running minimum.

    A = objs; O = optimum;
    if isempty(A) || size(A,2) ~= size(O,2)
        score = nan; return;
    end
    Nr   = size(O,1);
    N    = size(A,1);

    best = inf(Nr,1);
    for j = 1 : N
        diff = repmat(A(j,:),Nr,1) - O;     % Nr x M
        dj   = sqrt(sum(diff.^2, 2));       % Nr x 1
        best = min(best, dj);
    end
    score = mean(best);
end