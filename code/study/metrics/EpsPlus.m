function score = EpsPlus(objs, optimum)
%EpsPlus  Unary additive epsilon indicator I_{eps+}, matrix interface.
%   score = EpsPlus(objs, optimum)
%     objs    : (N x M) approximation set (minimization)
%     optimum : (Nr x M) reference set (sampled true PF)
%
%   I_{eps+}(A,R) = max_{r in R} min_{a in A} max_j ( a_j - r_j ).
%   Smallest additive eps such that every reference point is weakly
%   dominated by some shifted a. Smaller is better; weakly Pareto compliant
%   -> valid for monotonicity tracking (Sec. 4, axis ii).
%
%   OPTIMIZATION: loop over the N approximation points (small) rather than
%   the Nr reference points (large), keeping a running per-reference minimum
%   of  max_j (a_j - r_j). Same result, faster when N << Nr.
%
%   Reference: Zitzler et al., IEEE TEVC 2003.

    A = objs; O = optimum;
    if isempty(A) || size(A,2) ~= size(O,2)
        score = nan; return;
    end
    Nr = size(O,1);
    N  = size(A,1);

    bestPerRef = inf(Nr,1);                 % min over a of max_j(a_j - r_j)
    for j = 1 : N
        dj = max(repmat(A(j,:),Nr,1) - O, [], 2);   % Nr x 1
        bestPerRef = min(bestPerRef, dj);
    end
    score = max(bestPerRef);
end