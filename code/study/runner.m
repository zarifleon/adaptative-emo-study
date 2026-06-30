function runner()
%
%   Orchestrates:  algorithms x problems x M x seeds x mode
%   - logs full trajectory {Z^(t)},{P^(t)} per run via TrajLog
%   - fixes rng(seed) before each run for reproducibility
%   - skips runs whose .mat already exists (resumable: Ctrl-C and relaunch)
%   - records the effective N, the seed, and a 'completed' flag in each log
%
%   Run from the PlatEMO folder (where platemo.m lives), or adjust paths.
%
%   maxFE is fixed and identical across algorithms (fair budget). Population
%   N is requested per M; PlatEMO may adjust it and the EFFECTIVE N is what
%   TrajLog stores, so cross-algorithm N differences are detectable in data.

%% ----- paths -----
here = fileparts(mfilename('fullpath'));      % code/study
addpath(here);                                % so TrajLog/HV_ref resolve
% PlatEMO is assumed to be on the path already (run from its folder).

%% ----- experiment definition -----
% Problem families.
dtlzs = {@DTLZ2};
zcats = {@ZCAT1, @ZCAT20};
wfgs  = {};                 % e.g. {@WFG4} when you want to add them
mafs  = {};                 % e.g. {@MaF2}
imops = {};
sdtlzs= {};                 % SDTLZ1/2 = scaled, useful for the scale axis

final_problems = [dtlzs, zcats, wfgs, mafs, imops, sdtlzs];


%algorithms = {@NSGAIII, @ANSGAIII, @ARMOEA, @RVEA, @MOEAD, @MOEADUR, @AdaW };
algorithms = {@ARMOEA};

%Ms     = [10, 2, 3, 5];
Ms     = [10];
nSeeds = 1;            % seeds will be baseSeed+0 .. baseSeed+nSeeds-1
baseSeed = 42;

% Population size requested per M (PlatEMO adjusts to a valid simplex N).
% Keyed lookup: Npref(M).
Npref = containers.Map('KeyType','double','ValueType','double');
Npref(2)  = 100;
Npref(3)  = 105;
Npref(5)  = 126;
Npref(10) = 275;

maxFE = 100000;

% Normalization mode. Internal homogenization not implemented yet
% we run 'native' but keep the dimension in the tag
% so logs need not be renamed later.
modes = {'native'};

% Cadence: every generation.
setenv('TRAJLOG_DELTAT','1');

%% ----- output dir for logs (informational) -----
logdir = fullfile(here,'traj_logs');
if ~exist(logdir,'dir'); mkdir(logdir); end

%% ----- main sweep -----
nTotal = numel(algorithms)*numel(final_problems)*numel(Ms)*nSeeds*numel(modes);
fprintf('=== Smoke test: %d planned runs (maxFE=%d) ===\n', nTotal, maxFE);
done = 0; skipped = 0; failed = 0;

for ia = 1:numel(algorithms)
    algF = algorithms{ia};
    algN = func2str(algF); algN = strrep(algN,'@','');
    for ip = 1:numel(final_problems)
        probF = final_problems{ip};
        probN = func2str(probF); probN = strrep(probN,'@','');
        for im = 1:numel(Ms)
            M = Ms(im);
            if isKey(Npref,M); Nreq = Npref(M); else; Nreq = 100; end
            for is = 0:(nSeeds-1)
                seed = baseSeed + is;
                for ix = 1:numel(modes)
                    mode = modes{ix};

                    tag = sprintf('%s_%s_M%d_seed%d_%s', ...
                        algN, probN, M, seed, mode);
                    outfile = fullfile(logdir, [tag '.mat']);

                    % --- resume: skip if already completed ---
                    if exist(outfile,'file')
                        skipped = skipped + 1;
                        fprintf('skip  %s\n', tag);
                        continue;
                    end

                    % --- expose tag + seed to TrajLog via env ---
                    setenv('TRAJLOG_TAG', tag);
                    setenv('TRAJLOG_SEED', num2str(seed));

                    % --- reproducible run ---
                    rng(seed);
                    fprintf('run   %s ... ', tag);
                    t0 = tic;
                    try
                        platemo('algorithm', algF, ...
                            'problem',   probF, ...
                            'N',         Nreq, ...
                            'M',         M, ...
                            'maxFE',     maxFE, ...
                            'save',      0);
                        done = done + 1;
                        fprintf('ok (%.1fs)\n', toc(t0));
                    catch ME
                        % NotTerminated throws by design at budget end;
                        % that is NOT a failure. Distinguish real errors.
                        if strcmp(ME.identifier,'PlatEMO:Terminate') || ...
                                contains(ME.message,'PlatEMO') || ...
                                isempty(ME.identifier)
                            done = done + 1;
                            fprintf('ok* (%.1fs)\n', toc(t0));
                        else
                            failed = failed + 1;
                            fprintf('FAIL: %s\n', ME.message);
                        end
                    end
                end
            end
        end
    end
end

fprintf('=== done: %d run, %d skipped, %d failed (of %d) ===\n', ...
    done, skipped, failed, nTotal);
end