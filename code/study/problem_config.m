function cfg = problem_config(probName, M)
%PROBLEM_CONFIG  Per-problem overrides for metric computation.
%
%   cfg = problem_config(probName, M)
%   Compute_metrics calls
%   this once per run and uses the returned overrides. While a problem has no
%   entry, it falls back to the defaults below (the current "naive" setting).
%
%   Because metrics are post-processed, changing a value here and recomputing
%   does NOT require re-running any algorithm: delete that problem's CSVs in
%   csv_logs/ and re-run compute_metrics (skip-existing regenerates only them).
%
%   RETURNED STRUCT
%     cfg.hv_ref   : scalar factor on the normalized span (default 1.1).
%                    Used when cfg.hv_point is empty.
%     cfg.hv_point : (1 x M) ABSOLUTE reference point in original objective
%                    units. If set (non-empty), it OVERRIDES cfg.hv_ref.
%                    Provide per (problem, M) because the nadir scale changes
%                    with M.
%
%   HOW TO FILL 
%     - FACTOR (e.g. 1.05 for a given problem):
%         set cfg.hv_ref = 1.05 in that problem's case.
%     - ABSOLUTE POINT (e.g. [1.2 1.2 1.2] for M=3):
%         set cfg.hv_point = [1.2 1.2 1.2]; (match length to M).
%       For M-dependent points, branch on M inside the case.

% ---------- defaults (current naive setting) ----------
cfg = struct();
cfg.hv_ref   = 1.1;
cfg.hv_point = [];

% ---------- per-problem overrides ----------
% Add a case per problem reference point.
% Examples (commented; replace with real values):
switch upper(probName)

    % case 'DTLZ2'
    %     cfg.hv_ref = 1.1;                 % factor mode

    % case 'ZCAT20'
    %     if M == 3
    %         cfg.hv_point = [1.2 1.2 1.2]; % absolute mode, M=3
    %     elseif M == 10
    %         cfg.hv_point = repmat(1.2,1,10);
    %     end

    % case 'WFG1'
    %     cfg.hv_point = ...;

    otherwise
        % keep defaults
end
end