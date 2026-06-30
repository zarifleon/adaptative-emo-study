function score = HV_ref(PopObj,optimum,refFactor,refPoint)
%HV_ref  Hypervolume with an EXPOSED reference point (factor OR absolute).
%
%   score = HV_ref(PopObj,optimum,refFactor)            % factor mode
%   score = HV_ref(PopObj,optimum,[],refPoint)          % absolute-point mode
%
%   Same computation as PlatEMO's metric/HV.m, but the reference point is
%   configurable so the study can recompute HV per-problem
%   choice WITHOUT re-running any experiment.
%
%   INPUTS
%     PopObj    : (N x M) approximation set objective vectors.
%     optimum   : (K x M) sampled true Pareto front (for nadir estimation).
%     refFactor : scalar multiplier on the normalized span (PlatEMO default
%                 1.1). Used when refPoint is empty/omitted.
%     refPoint  : (1 x M) ABSOLUTE reference point in ORIGINAL objective
%                 units. If provided, it overrides refFactor: the objective
%                 space is normalized so this point maps to ones(1,M).
%
%   MODES
%     factor mode   : normalize by (fmax-fmin)*refFactor, fmax=max(optimum),
%                     fmin=min(min(PopObj),0). RefPoint = ones (PlatEMO).
%     absolute mode : normalize by (refPoint - fmin). Points beyond refPoint
%                     (>1 after normalization) are discarded, exactly as the
%                     factor mode discards points beyond the ref point.
%

    if nargin < 3; refFactor = 1.1; end
    if nargin < 4; refPoint  = []; end

    if size(PopObj,2) ~= size(optimum,2)
        score = nan; return;
    end
    [N,M] = size(PopObj);
    fmin  = min(min(PopObj,[],1),zeros(1,M));

    if ~isempty(refPoint)
        % ----- absolute reference point mode -----
        rp   = refPoint(:)';
        span = rp - fmin;
        span(span<=0) = 1;
        PopObj = (PopObj - repmat(fmin,N,1)) ./ repmat(span,N,1);
    else
        % ----- factor mode (PlatEMO-compatible) -----
        fmax   = max(optimum,[],1);
        PopObj = (PopObj-repmat(fmin,N,1))./repmat((fmax-fmin)*refFactor,N,1);
    end

    PopObj(any(PopObj>1,2),:) = [];
    RefPoint = ones(1,M);
    if isempty(PopObj)
        score = 0;
    elseif M < 4
        pl = sortrows(PopObj);
        S  = {1,pl};
        for k = 1 : M-1
            S_ = {};
            for i = 1 : size(S,1)
                Stemp = Slice(cell2mat(S(i,2)),k,RefPoint);
                for j = 1 : size(Stemp,1)
                    temp(1) = {cell2mat(Stemp(j,1))*cell2mat(S(i,1))};
                    temp(2) = Stemp(j,2);
                    S_      = Add(temp,S_);
                end
            end
            S = S_;
        end
        score = 0;
        for i = 1 : size(S,1)
            p     = Head(cell2mat(S(i,2)));
            score = score + cell2mat(S(i,1))*abs(p(M)-RefPoint(M));
        end
    else
        SampleNum = 1e6;
        MaxValue  = RefPoint;
        MinValue  = min(PopObj,[],1);
        Samples   = unifrnd(repmat(MinValue,SampleNum,1),repmat(MaxValue,SampleNum,1));
        for i = 1 : size(PopObj,1)
            domi = true(size(Samples,1),1);
            m    = 1;
            while m <= M && any(domi)
                domi = domi & PopObj(i,m) <= Samples(:,m);
                m    = m + 1;
            end
            Samples(domi,:) = [];
        end
        score = prod(MaxValue-MinValue)*(1-size(Samples,1)/SampleNum);
    end
end

function S = Slice(pl,k,RefPoint)
    p  = Head(pl);  pl = Tail(pl);  ql = [];  S = {};
    while ~isempty(pl)
        ql  = Insert(p,k+1,ql);
        p_  = Head(pl);
        cell_(1,1) = {abs(p(k)-p_(k))};
        cell_(1,2) = {ql};
        S   = Add(cell_,S);
        p   = p_;  pl = Tail(pl);
    end
    ql = Insert(p,k+1,ql);
    cell_(1,1) = {abs(p(k)-RefPoint(k))};
    cell_(1,2) = {ql};
    S  = Add(cell_,S);
end

function ql = Insert(p,k,pl)
    flag1 = 0; flag2 = 0; ql = [];  hp = Head(pl);
    while ~isempty(pl) && hp(k) < p(k)
        ql = [ql;hp];  pl = Tail(pl);  hp = Head(pl);
    end
    ql = [ql;p];  m = length(p);
    while ~isempty(pl)
        q = Head(pl);
        for i = k : m
            if p(i) < q(i); flag1 = 1; else; if p(i) > q(i); flag2 = 1; end; end
        end
        if ~(flag1 == 1 && flag2 == 0); ql = [ql;Head(pl)]; end
        pl = Tail(pl);
    end
end

function p = Head(pl)
    if isempty(pl); p = []; else; p = pl(1,:); end
end

function ql = Tail(pl)
    if size(pl,1) < 2; ql = []; else; ql = pl(2:end,:); end
end

function S_ = Add(cell_,S)
    n = size(S,1);  m = 0;
    for k = 1 : n
        if isequal(cell_(1,2),S(k,2))
            S(k,1) = {cell2mat(S(k,1))+cell2mat(cell_(1,1))};
            m = 1;  break;
        end
    end
    if m == 0; S(n+1,:) = cell_(1,:); end
    S_ = S;
end