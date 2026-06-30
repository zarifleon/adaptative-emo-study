function score = DeltaP_mat(objs, optimum)
%DeltaP_mat  Averaged Hausdorff distance Delta_p (p=1), matrix interface.
%   Identical to PlatEMO metric/DeltaP.m: max( GD-like, IGD-like ).
%
%   Delta_p = max( mean_a min_r ||a-r||,  mean_r min_a ||a-r|| ).
%   Both directions come from ONE distance matrix. pdist2 is MATLAB's
%   compiled pairwise-distance routine (fast); we compute it once and reduce
%   along both axes, instead of two separate loops.

    A = objs; O = optimum;
    if isempty(A) || size(A,2) ~= size(O,2)
        score = nan; return;
    end
    D = pdist2(O, A);                       % Nr x N, computed once
    igd_like = mean(min(D, [], 2));         % mean over references of nearest a
    gd_like  = mean(min(D, [], 1));         % mean over a of nearest reference
    score    = max(gd_like, igd_like);
end