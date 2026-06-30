function compute_metrics(varargin)
%COMPUTE_METRICS  Post-process trajectory logs into per-run tidy metrics CSVs.
%
%   compute_metrics()                 % process all logs, default options
%   compute_metrics('hv_ref',1.1,'refN',10000,'s_exp',[],'r2N',100)
%   compute_metrics('combine',true)   % also merge all per-run CSVs into one
%
%   For every  traj_logs/<tag>.mat  this writes  csv_logs/<tag>.csv  (one
%   tidy file per run). PER-FILE OUTPUT gives:
%     - resumability: an existing <tag>.csv is skipped, so an interrupted
%       pass resumes where it stopped;
%     - parallelism: several MATLAB instances can run this over disjoint
%       subsets without colliding (same pattern as the runner);
%     - smaller memory: rows are flushed per run, not held for all runs.
%
%   This is the MATLAB->Python handoff: ALL numeric values are computed here
%   (comparability); Python only reads the CSVs to aggregate/plot.
%
%   Metric panel:
%     on P^(t): HV, IGDp, DeltaP, IGD, EpsPlus, R2, SEnergy(+raw,+nEff), SolowPolasky
%     on Z^(t): IGDp, EpsPlus, R2, SEnergy(+raw,+nEff), SolowPolasky

    %% ----- options -----
    p = inputParser;
    addParameter(p,'hv_ref',1.1);
    addParameter(p,'refN',10000);
    addParameter(p,'s_exp',[]);
    addParameter(p,'r2N',100);
    addParameter(p,'sp_theta',10);
    addParameter(p,'logdir','');
    addParameter(p,'csvdir','');
    addParameter(p,'combine',false);
    addParameter(p,'overwrite',false);
    parse(p,varargin{:});
    opt = p.Results;

    here = fileparts(mfilename('fullpath'));
    studyDir = fileparts(here);
    addpath(here); addpath(studyDir);

    if isempty(opt.logdir); opt.logdir = fullfile(studyDir,'traj_logs'); end
    if isempty(opt.csvdir); opt.csvdir = fullfile(studyDir,'csv_logs');  end
    if ~exist(opt.csvdir,'dir'); mkdir(opt.csvdir); end

    files = dir(fullfile(opt.logdir,'*.mat'));
    if isempty(files)
        fprintf('No .mat logs in %s\n', opt.logdir); return;
    end
    fprintf('Found %d logs. Output dir: %s\n', numel(files), opt.csvdir);

    pfCache = containers.Map('KeyType','char','ValueType','any');
    nDone = 0; nSkip = 0; nIncomplete = 0;

    for fi = 1:numel(files)
        matName = files(fi).name;
        [~,base] = fileparts(matName);
        outcsv   = fullfile(opt.csvdir, [base '.csv']);

        if exist(outcsv,'file') && ~opt.overwrite
            nSkip = nSkip + 1;
            continue;
        end

        S = load(fullfile(opt.logdir, matName));
        if ~isfield(S,'traj'); continue; end
        t = S.traj;
        if isfield(t,'completed') && ~t.completed
            nIncomplete = nIncomplete + 1;
            fprintf('  skip incomplete: %s\n', matName); continue;
        end

        M    = double(t.M);
        prob = char(t.prob);
        alg  = char(t.alg);
        seed = double(t.seed);
        mode = parseMode(char(t.tag));

        key = sprintf('%s_M%d', prob, M);
        if isKey(pfCache,key); PF = pfCache(key);
        else; PF = getReferenceFront(prob, M, opt.refN); pfCache(key) = PF; end

        Wr2      = UniformPoint(opt.r2N, M);
        refRange = [min(PF,[],1); max(PF,[],1)];

        % per-problem HV reference override (defaults to factor)
        pcfg = problem_config(prob, M);

        rows = {};
        K = numel(t.FE);
        for k = 1:K
            P = t.Pobj{k}; Z = t.Z{k};
            gen = t.gen(k); FE = t.FE(k);

            add('P','HV',         HV_ref(P,PF,pcfg.hv_ref,pcfg.hv_point));
            add('P','IGDp',       IGDp_mat(P,PF));
            add('P','DeltaP',     DeltaP_mat(P,PF));
            add('P','IGD',        IGD_mat(P,PF));
            add('P','EpsPlus',    EpsPlus(P,PF));
            add('P','R2',         R2_ind(P,Wr2,[]));
            [pp,raw,ne] = SEnergy(P,opt.s_exp,refRange);
            add('P','SEnergy',pp); add('P','SEnergy_raw',raw); add('P','SEnergy_nEff',ne);
            add('P','SolowPolasky', SolowPolasky(P,opt.sp_theta));

            add('Z','IGDp',       IGDp_mat(Z,PF));
            add('Z','EpsPlus',    EpsPlus(Z,PF));
            add('Z','R2',         R2_ind(Z,Wr2,[]));
            [zp,zraw,zne] = SEnergy(Z,opt.s_exp,refRange);
            add('Z','SEnergy',zp); add('Z','SEnergy_raw',zraw); add('Z','SEnergy_nEff',zne);
            add('Z','SolowPolasky', SolowPolasky(Z,opt.sp_theta));
        end

        T = cell2table(rows, 'VariableNames', ...
            {'alg','prob','M','seed','mode','gen','FE','target','metric','value'});
        tmp = [outcsv '.tmp'];
        writetable(T, tmp);
        movefile(tmp, outcsv, 'f');
        nDone = nDone + 1;
        fprintf('  [%d/%d] %s : %d steps\n', fi, numel(files), base, K);
    end

    fprintf('=== metrics: %d written, %d skipped, %d incomplete ===\n', ...
            nDone, nSkip, nIncomplete);

    if opt.combine
        combineCSVs(opt.csvdir, fullfile(studyDir,'metrics_all.csv'));
    end

    function add(target,metric,value)
        rows(end+1,:) = {alg,prob,M,seed,mode,gen,FE,target,metric,value}; %#ok<AGROW>
    end
end

function combineCSVs(csvdir, outfile)
    files = dir(fullfile(csvdir,'*.csv'));
    if isempty(files); fprintf('No per-run CSVs to combine.\n'); return; end
    allT = [];
    for i = 1:numel(files)
        Ti = readtable(fullfile(csvdir, files(i).name));
        allT = [allT; Ti]; %#ok<AGROW>
    end
    writetable(allT, outfile);
    fprintf('Combined %d CSVs -> %s (%d rows)\n', numel(files), outfile, height(allT));
end

function PF = getReferenceFront(probName, M, refN)
    try
        probFcn = str2func(probName);
        Prob = probFcn('M', M);
        PF   = Prob.GetOptimum(refN);
    catch
        try
            Prob = probFcn('M', M); PF = Prob.optimum;
        catch err
            warning('Could not get PF for %s M=%d: %s', probName, M, err.message);
            PF = nan(1, M);
        end
    end
end

function mode = parseMode(tag)
    parts = strsplit(tag,'_');
    if isempty(parts); mode = 'native'; else; mode = parts{end}; end
end