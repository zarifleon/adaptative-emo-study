function TrajLog(action, varargin)
%TrajLog - Trajectory logger for the adaptive-reference-set study.
%
%   Captures {Z^(t)} (reference set) and {P^(t)} (population objectives)
%   across generations without modifying PlatEMO's core. Driven by one call
%   per generation from inside each algorithm's main loop, plus a begin and
%   an onCleanup-guarded end.
%
%   ACTIONS
%     TrajLog('begin', Algorithm, Problem, tag, seed)
%         Call once, right after the population is initialized.
%         seed : the rng seed used for this run (stored in the log).
%
%     TrajLog('rec', Algorithm, Problem, Population, Z)
%         Call once per generation, after Z is (re)computed.
%         Z : THIS algorithm's reference set (Z / RefPoint / W / V ...),
%             passed as a raw numeric matrix (N x M).
%
%     TrajLog('end', Algorithm, Problem)
%         Call via  onCleanup(@() TrajLog('end',Algorithm,Problem))
%         so it fires no matter how main() exits (normal, NotTerminated
%         exception, Ctrl-C, or error).
%
%   COMPLETION / ATOMICITY
%     During the run, records are buffered in memory. On 'end', the buffer
%     is written to  <tag>.mat.tmp  and renamed to  <tag>.mat  ONLY if the
%     run actually exhausted its evaluation budget (Problem.FE >= maxFE).
%     An interrupted run (Ctrl-C, crash) leaves no <tag>.mat, so the runner
%     re-runs it. A 'completed' flag is also stored inside the struct.
%
%   CADENCE
%     Every DELTA_T-th generation is stored (env TRAJLOG_DELTAT, default 1
%     = every generation). The first generation is always stored.
%
%   OUTPUT  code/study/traj_logs/<tag>.mat  with struct  traj :
%     .tag .alg .prob .M .D .N .seed        metadata
%     .completed (logical)                  true iff budget was exhausted
%     .FE   (1xK)   function evals at each stored step
%     .gen  (1xK)   generation index at each stored step
%     .numGen       total generations run (= last gen index)
%     .Zcard (1xK)  |Z^(t)| at each stored step
%     .Z    (1xK cell)  each cell (|Z| x M) reference set
%     .Pobj (1xK cell)  each cell (|P| x M) population objectives
%     .runtime      wall-clock seconds (Algorithm.metric.runtime)

    persistent S

    switch lower(action)
        % -----------------------------------------------------------------
        case 'begin'
            Algorithm = varargin{1};
            Problem   = varargin{2};
            tag       = varargin{3};
            if numel(varargin) >= 4 && ~isempty(varargin{4})
                seed = varargin{4};
            else
                sv = getenv('TRAJLOG_SEED');
                if ~isempty(sv); seed = str2double(sv); else; seed = NaN; end
            end

            S = struct();
            S.tag     = tag;
            S.alg     = class(Algorithm);
            S.prob    = class(Problem);
            S.M       = Problem.M;
            S.D       = Problem.D;
            S.N       = Problem.N;        % effective N after PlatEMO adjustment
            S.seed    = seed;
            S.FE      = [];
            S.gen     = [];
            S.Zcard   = [];
            S.Z       = {};
            S.Pobj    = {};
            S.callcnt = 0;
            S.gencnt  = 0;

        % -----------------------------------------------------------------
        case 'rec'
            if isempty(S); return; end   % defensive: begin missing
            Problem    = varargin{2};
            Population = varargin{3};
            Z          = varargin{4};

            DELTA_T   = getDeltaT();
            S.callcnt = S.callcnt + 1;
            S.gencnt  = S.gencnt + 1;

            if mod(S.callcnt-1, DELTA_T) == 0
                S.gen(end+1)   = S.gencnt;
                S.FE(end+1)    = Problem.FE;
                Zmat           = double(Z);
                S.Z{end+1}     = Zmat;
                S.Zcard(end+1) = size(Zmat,1);
                S.Pobj{end+1}  = Population.objs;
            end

        % -----------------------------------------------------------------
        case 'end'
            if isempty(S); return; end
            Algorithm = varargin{1};
            Problem   = varargin{2};

            completed = Problem.FE >= Problem.maxFE;
            S.completed = completed;
            S.numGen    = S.gencnt;
            if isstruct(Algorithm.metric) && isfield(Algorithm.metric,'runtime')
                S.runtime = Algorithm.metric.runtime;
            else
                S.runtime = NaN;
            end

            traj = rmfield(S, {'callcnt','gencnt'});

            outdir = fullfile(fileparts(mfilename('fullpath')), 'traj_logs');
            if ~exist(outdir,'dir'); mkdir(outdir); end
            base    = sanitize(S.tag);
            tmpfile = fullfile(outdir, [base '.mat.tmp']);
            outfile = fullfile(outdir, [base '.mat']);

            S = [];   % reset persistent BEFORE I/O so a failed save can't leave stale state

            if ~completed
                % interrupted: do not produce a definitive .mat
                fprintf('[TrajLog] run INCOMPLETE (FE<%g); no .mat written for %s\n', ...
                        Problem.maxFE, base);
                return;
            end

            save(tmpfile, 'traj', '-v7');
            movefile(tmpfile, outfile, 'f');   % atomic-ish promotion
            fprintf('[TrajLog] wrote %s  (%d steps, %d gens, completed)\n', ...
                    outfile, numel(traj.FE), traj.numGen);

        otherwise
            error('TrajLog:badaction','Unknown action "%s".', action);
    end
end

% =========================================================================
function dt = getDeltaT()
    v = getenv('TRAJLOG_DELTAT');
    if ~isempty(v)
        dt = str2double(v);
        if isnan(dt) || dt < 1; dt = 1; end
    else
        dt = 1;     % default: every generation
    end
end

% =========================================================================
function s = sanitize(s)
    if isempty(s); s = 'UNTAGGED'; end
    s = regexprep(s, '[^\w\-\.]', '_');
end