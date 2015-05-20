function [prm, x, x_true, E, stats, EStore] = main(prm, x_true0, E0, stats0)
%% function [prm, x, x_true, E, stats, EStore] = main(prm, x_true0, E0, stats0)
%
% Runs the system using parameters contained in a structure variable.
%
% @param prm - system parametes (structure, see get_prmstruct.m)
% @param x_true0 - initial true field (n x 1, optional)
% @param E0 - initial ensemble (n x m, optional)
% @param stats0 - initial statistics (structure, optional)
% @param prm - updated prm structure
% @return x - final analysis (n x 1)
% @return x_true - final true field (n x 1)
% @return E - final ensemble (n x m)
% @return stats - final statistics (structure) [note: not used for real_world cases].
% @return EStore - total ensemble history
%
% Usage:
%
% [prm, x, x_true, E, stats] = main(prm)
%   - conduct a run with a freshly generated ensemble
% [prm, x, x_true, E, stats] = main(prm, x_true0)
%   - conduct a run with a freshly generated ensemble and specified true field
% [prm, x, x_true, E, stats] = main(prm, x_true0, E0)
%   - conduct a run with specified true field and ensemble
% [prm, x, x_true, E, stats] = main(prm, x_true0, E0, stats0)
%   - conduct a run with specified true field and ensemble; and append the
%     existing statics rather than start it from scratch
%
% File:           main.m
%
% Created:        23/03/2009
%
% Last modified:  14/02/2013
%
% Author:         Tobias Siegfried (TS, hydrosolutions ltd.), based on the
%                 original EnKF-Matlab code by Pavel Sakov (PS, NERSC).           
%
% Purpose:        Runs the system using parameters contained in a structure
%                 variable. % (TS) The state here is an augmented state-parameter vector with
%                 x = [Q; S; G; SN; ETact; alpha1; alpha2; alpha3; Smax; T0; Tmelt]; (11-D)
%                 if prm.customprm.snow == 1. if prm.customprm.snow == 0, then
%                 x = [Q; S; G; ETact; alpha1; alpha2; alpha3; Smax]; (8-D)
%
% Description:    Main file of the data assimilation system.
%
% Revisions:      21.12.2012 TS: starting work on assimilating an arbitrary
%                 number of catchments simultaneously.
%                 13.12.2012 TS: accounting for the snow / no snow option.
%                 30.11.2012 TS: starting implementation of variable (in time)
%                 observations - DONE: 10.12.2012
%                 06.11.2012 TS: adaptation to iMoMo
%                 29.10.09 PS:
%                   - Swapped integration and observation+assimilation in the
%                     main cycle. Replaced (step - 1) by (step) in tests on
%                     whether observation or assimilation should be conducted.
%                     Introduced reporting of accumulated time spent on
%                     integration and assimilation.
%                 1.11.2009 PS:
%                   - added measurement of total time for assimilation and
%                     propagation
%                 25/08/2010 PS:
%                   - This file, formerly known as main_a.m, now replaced the
%                     previous (synchronous) version of main.m
% 
% Copyright (C) 2014 hydrosolutions ltd., Lindenbachstrasse 11, CH-8006
% Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.
%% PREPARATIONS
% number of catchments
prm.customprm.nC = size(E0,3);
% storage for ensembles
EStore = zeros(prm.n * prm.customprm.nC, prm.m, prm.nstep);
m = prm.m; % m = number of ensembles
% check the resume state later!
resumed = 0;
% path
% prm = setpath(prm); % deleted, 29/04/2013

rand('state', prm.seed); % initialise to seed
randn('state', prm.seed); % initialise to seed

if ~exist('E0', 'var') || isempty(E0)
    % START
    [x_true, E, t, prm] = generate(prm);  %%function generate(...) does not exist, E0 anyway always exists
    if exist('x_true0', 'var') && ~isempty(x_true0)
        x_true = x_true0;
    end
else
    % RESUME
    resumed = 1;
    if ~prm.realWorld, x_true = x_true0; end
    E = E0(:, 1 : m, :); % E: n x m - matrix of an ensemble of model states
    if exist('stats0', 'var') && ~isempty(stats0)
        rand('state', stats0.rand_finalstate);
        randn('state', stats0.randn_finalstate);
    end
end

if exist('stats0', 'var') && ~isempty(stats0)
    t = stats0.t(end);
    step0 = stats0.step(end) + 1;
else
    stats0 = [];
    step0 = 1;
    if exist('E0', 'var')
        t = 0;
    end
end
step1 = step0 + prm.nstep - 1;
prm.current.step = step0;

% here, for the multi-catchment case, simply reshape E.
rsE = squeeze(E(:,:,1));
for i = 1 : prm.customprm.nC - 1
    rsE = [rsE;squeeze(E(:,:,i+1))]; %n.C times states, n. of ensembles -> "usual" shape of E
end
E = abs(rsE); % abs() - start with consistent meaningful values

if ~prm.realWorld, x_true = x_true(:); end
%% ENSEMBLE ANOMALIES
x = mean(E,2); % x: mean over ensemble members - E state-vec matrix!
A = E - repmat(x, 1, m); % ensemble anomalies.

if resumed && ~isempty(stats0) && ~prm.realWorld
    stats = stats0;
    % if this is a resumed run, stats.lasttime will have the value of 1
    % perhaps, which may need to be corrected
    stats.lasttime = floor(step1 / prm.assim_step) == floor(step0 / prm.assim_step);
elseif ~prm.realWorld
    % calculate first forecast and analysis stats from the initial state
    stats = calc_stats(prm, x_true, x, A, step0, step1, t, stats0); % Updates Filter statistics, stats0 empty because first forecast.
    stats = calc_stats(prm, x_true, x, A, step0, step1, t, stats);
end

y = []; HE = []; pos = [];
%% OBSERVATION - CALCULATE H
% ToDo: Dynamic updating once live iMoMo sensed data comes in!
if strcmp(prm.obs_spacing, 'regular') || strcmp(prm.obs_spacing, 'variable')
    [prm,pos_now] = get_pos(prm);
    H_now = calc_h(prm, pos_now);
end
%% MAIN CYCLE PREPARATION

if prm.verbose > 1
    fprintf('  main cycle:');
end

prm.trueX = []; prm.obsX = [];

prm.crit1=NaN(1,step1); % for truncating.
prm.crit2=NaN(1,step1); % ditto
prm.crit3=NaN(1,step1); % ditto

if prm.customprm.snow
    prm.obs_variance = [prm.customprm.obsVarQ prm.customprm.obsVarS ...
        prm.customprm.obsVarG prm.customprm.obsVarSN prm.customprm.obsVarET]';
else
    prm.obs_variance = [prm.customprm.obsVarQ prm.customprm.obsVarS ...
        prm.customprm.obsVarG prm.customprm.obsVarET]';
    prm.obs_variance = [prm.obs_variance; zeros(prm.nprm,1)];
    prm.obs_variance = repmat(prm.obs_variance,1,prm.customprm.nC);
    prm.obs_variance = prm.obs_variance(:);
end
%% MAIN CYCLE
for step = step0 : step1

    prm.current.step = step;
    
    if prm.current.step == 1
        prm.storage.SG = zeros(prm.customprm.nStorage * prm.customprm.nC, prm.nstep, m);
        % here, we really need to figure out a way to initialize S0 / G0 properly. 
        if prm.customprm.snow % this needs to be updated for snow!
            prm.storage.SG(:,step,:) = repmat([prm.customprm.S0; ...
                prm.customprm.G0; prm.customprm.SN0],[prm.customprm.nC 1 m]);
        else 
            S0 = load(prm.path.S0G0,'S0');
            prm.customprm.S0 = S0.S0';
            G0 = load(prm.path.S0G0,'G0')';
            prm.customprm.G0 = G0.G0';
            S0G0 = [S0.S0';G0.G0']; S0G0 = S0G0(:);
            prm.storage.SG(:,step,:) = repmat(S0G0,[1 1 m]);
        end
        % Time
        ye = num2str(prm.customprm.simStartY);
        mo = sprintf('%02d\n',prm.customprm.simStartM);
        da = sprintf('%02d\n',prm.customprm.simStartD);
        prm.current.serialDate = DateStr2Num(strcat(ye,'-',mo,'-',da),29);
        prm.customprm.serialDateStart = prm.current.serialDate;
    end
    
    prm.current.SerialDateVec = datevec(prm.current.serialDate)';
    prm.current.doy = datenum2doy(prm.current.SerialDateVec);
    %% STATE PROPAGATION FOR ENSEMBLE
    [E,t] = model_step(prm, E, t, 0); % Performs one model step.
    
    if ~prm.customprm.snow
        prm.storage.SG(1:2:end-1,step,:) = E(2:8:end,:); % S
        prm.storage.SG(2:2:end,step,:) = E(3:8:end,:); % G
    else % to check still
        prm.storage.SG(1:3:end-2,step,:) = E(2:5:(prm.n-prm.p)*prm.customprm.nC-3,:); % S
        prm.storage.SG(2:3:end-1,step,:) = E(3:5:(prm.n-prm.p)*prm.customprm.nC-2,:); % G
        prm.storage.SG(3:3:end,step,:) = E(4:5:(prm.n-prm.p)*prm.customprm.nC-1,:); % G
    end
    
    if ~isreal(E) || any(isinf(E(:))) || any(isnan(E(:))) || any(E(:)<0)
        disp('ERROR in ENSEMBLE! ? CHECK!'); keyboard;
    end
    
    %% MEAN FIELD CHECK
    r = max(abs(x));
    if isnan(r) || r > 1e10 % fudge factor here - check and hardwire and parameters (TS)
        if prm.verbose
            fprintf('\n    warning: main(): mean field: invalid magnitude at t = %f, step = %d\n', t, step);
            fprintf('    max(||x||) = %.3g\n', r);
            fprintf('    bailing out...\n');
        end
        break;
    end
    
    %% TRUE FIELD (in lieu of observations for testing)
    if ~prm.realWorld
        [x_true, t] = model_step(prm, x_true, t, 1);
        prm.trueX = [prm.trueX x_true];
        
        %% TRUE FIELD CHECK
        r = max(abs(x_true));
        if isnan(r) || r > 1e10 % fudge factor here - check and hardwire and parameters (TS)
            if prm.verbose
                fprintf('\n    warning: main(): true field: invalid magnitude at t = %f, step = %d\n', t, step);
                fprintf('    max(||x_true||) = %.3g\n', r);
                fprintf('    bailing out...\n');
            end
            break;
        end
    end
    %% OBSERVATIONS every prm.obs_step step
    if ~mod(step, prm.obs_step) && step < prm.assim_stop % prm.obs_step by default set to prm.assim_step
        
        %
        if ~mod(step, prm.assim_step)
            fprintf('O%d - ', step);
        end
        
        % in case of time varying obs locations -- calculate H
        if ~strcmp(prm.obs_spacing, 'regular')
            [prm,pos_now] = get_pos(prm);
            H_now = calc_h(prm, pos_now);
        end
        
        if ~prm.realWorld
            if ~isempty(stats.randn_obsstate)
                randn('state', stats.randn_obsstate);
            end
        end
        
        if ~prm.realWorld
            y_now = model_getobs(prm, x_true, H_now); % y_now: actual vector of observation
        else
            y_now = zeros(length(pos_now),1);
            for o = 1 : length(pos_now)
                tt = prm.current.serialDate - prm.customprm.serialDateStart + 1;
                if prm.customprm.snow % not yet adapated to multi-catchment assimilation.
                    if pos_now(o) == 1
                        y_now(o,1) = prm.customprm.Q(tt);
                    elseif pos_now(o) == 2
                        y_now(o,1) = prm.customprm.SM(tt);
                    elseif pos_now(o) == 3
                        y_now(o,1) = prm.customprm.G(tt);
                    elseif pos_now(o) == 4
                        y_now(o,1) = prm.customprm.SN(tt);
                    elseif pos_now(o) == 5
                        y_now(o,1) = prm.customprm.ET(tt);
                    end
                else % no snow - not elegant but working and not generizable.
                    obs = [prm.customprm.Q(tt,:); prm.customprm.SM(tt,prm.customprm.cID); ...
                        prm.customprm.G(tt,:); prm.customprm.ET(tt,:)];
                    obs = [obs;zeros(prm.n-prm.nprm,prm.customprm.nC)]; obs = obs(:);
                    y_now = obs(pos_now);
                end
            end
        end
        
        stats.randn_obsstate = randn('state'); % save the state of randn()
        y = [y; y_now]; % collects and stores actual vectors of observations
        prm.obsX= [prm.obsX y];
        
        % Default procedure for getting ensemble observations.
        HE_now = model_getHE(prm, H_now, E); % calculated states of obsveration points
        HE = [HE; HE_now];
        pos = [pos; pos_now];
        
    end
    
    %% (iii) ASSIMILATE every prm.assim_step model steps;
    if mod(step, prm.assim_step) == 0 && step < prm. assim_stop
        
        if ~mod(step, 1), fprintf(' A%d - ', step); end
        
        % calculate ensemble anomalies if necessary
        x = mean(E,2);
        A = E - repmat(x, 1, m);
        
        % calculate forecast stats
        %stats = calc_stats(prm, x_true, x, A, step, step1, t, stats); % note: step1 is last time whereas step is current step
        
        if ~isempty(y) % there are observations, so assimilate
            Hx = mean(HE,2); % mean of ensemble observations
            dy = y - Hx; % actual vec. of observations
            HA = HE - repmat(Hx, 1, m); % ensemble observation anomalies
            [dx, A] = assimilateEnKF(prm, A, HA, pos, dy, stats);
            
            x = x + dx;
            E = A + repmat(x, 1, m);
            
            %% CHECK whether states and parameters are physically meaningful
            % I. alpha1, alpha2 and d (groundwater parameter)
            idP = repmat((prm.p+1:prm.p+3)',prm.customprm.nC,1);
            addID = repmat(1:prm.customprm.nC,3,1); addID = addID(:) - 1;
            idP = idP + addID * prm.n;
            EE = E(idP,:);
            EE = truncateE(EE, prm,0.001,.999, step); % hardwired bounds here
            E(idP,:) = EE;
            
            if prm.customprm.snow % Smax and Tmelt
                EE = E([prm.p+4 prm.p+prm.nprm],:);
                EE = truncateE(EE, prm,0,realmax, step);
                E([prm.p+4 prm.p+prm.nprm],:) = EE;
            else % Smax only
                idP = repmat((prm.nprm+prm.p)',prm.customprm.nC,1);
                addID = repmat(1:prm.customprm.nC,1,1); addID = addID(:) - 1;
                idP = idP + addID * prm.n;
                EE = E(idP,:);
                EE = truncateE(EE, prm,0.001,realmax, step);
                E(idP,:) = EE;
            end
            % II. CHECK: any unphysical states
            idP = repmat((1:4)',prm.customprm.nC,1);
            addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
            idP = idP + addID * prm.n;
            tempE = E(idP,:);
            tempE(tempE<0)=0;
            E(idP,:) = tempE;
            % III. CHECK: Smax >= S(t)
            if step > 1
                E = SmaxCheck(E,prm,step);
            end
            
            % update storage so as to reflect ensemble corrections
            if ~prm.customprm.snow
                prm.storage.SG(1:2:end-1,step,:) = E(2:8:end,:); % S
                prm.storage.SG(2:2:end,step,:) = E(3:8:end,:);   % G
            else % to check still
                prm.storage.SG(1:3:end-2,step,:) = E(2:5:(prm.n-prm.p)*prm.customprm.nC-3,:); % S
                prm.storage.SG(2:3:end-1,step,:) = E(3:5:(prm.n-prm.p)*prm.customprm.nC-2,:); % G
                prm.storage.SG(3:3:end,step,:) = E(4:5:(prm.n-prm.p)*prm.customprm.nC-1,:);   % SN
            end
            
            % redo A for stats calc and check on uncrecov. divergence.
            x = mean(E,2);
            A = E - repmat(x, 1, m);
            
        end % if ~isempty(y) % there are obs
        
        % calculate analysis stats
        if ~prm.realWorld
            stats = calc_stats(prm, x_true, x, A, step, step1, t, stats);
        end
        %% check for crash or "unrecoverable" divergence
        if ~prm.realWorld
            if isnan(stats.rmse_a(1, end)) || (prm.rmse_max > 0 && stats.rmse_a(1, end) > prm.rmse_max)
                if prm.verbose
                    fprintf('\n    warning: main(): rmse = %.3g > rmse_max = %.3g at t = %f, step = %d\n', stats.rmse_a(1, end), prm.rmse_max, t, step);
                    fprintf('    bailing out...\n');
                end
                break;
            end
        end
        
        y = [];
        HE = [];
        pos = [];
        
        if prm.verbose > 1
            fprintf('.');
        end
    end % assimilate
    
    %% Display output
    if ~mod(step,10)
        disp(['end step: ' num2str(step)]);
    end
    
    EStore(:,:,step) = E;
    
    if ~(step==step1) % don't update at the last step!
        prm.current.serialDate = prm.current.serialDate + prm.customprm.dt;
    end
end % step = step0 : step1

prm.crit1 = prm.crit1; prm.crit2 = prm.crit2; prm.crit3 = prm.crit3;

%% print output
if prm.verbose && ~prm.realWorld
    fprintf('\n');
    i2 = length(stats.step);
    i1 = max(1, i2 - prm.plot_diag_np + 1);
    fprintf('  a fingerprint:\n', i2 - i1 + 1);
    fprintf('    <analysis rmse, step=%d-%d> = %.10g\n', stats.step(i1), stats.step(i2), mean(stats.rmse_a(1, i1 : i2)));
    fprintf('    <analysis spread, step=%d-%d> = %.10g\n', stats.step(i1), stats.step(i2), mean(stats.spread_a(1, i1 : i2)));
    fprintf('    <last analysis rmse> = %.10g\n', stats.rmse_a(end));
    fprintf('    <last analysis spread> = %.10g\n', stats.spread_a(end));
end
%% final state estimate
x = mean(E,2);

stats.rand_finalstate = rand('state');
stats.randn_finalstate = randn('state');

if prm.realWorld, x_true = []; end % just to make sure that the algo does not complain.

% write back last S and G for later reuse as S0 and G0
prm.customprm.S0 = mean(squeeze(prm.storage.SG(1:2:end-1,step,:)),2)'; S0 = prm.customprm.S0;
prm.customprm.G0 = mean(squeeze(prm.storage.SG(2:2:end,step,:)),2)'; G0 = prm.customprm.G0;
SG = squeeze(prm.storage.SG(:,end,:));
save(prm.path.S0G0,'S0');
save(prm.path.S0G0,'G0','-append');
save(prm.path.S0G0,'SG','-append');

fileID = fopen([prm.path.S0G0(1:end-8) 'simTime.txt'],'w');
fprintf(fileID,'%e %e %e %e %e %e',datevec(prm.current.serialDate));
fclose(fileID);

return

%% NESTED FUNCTIONS
    function [stats] = calc_stats(prm, x_true, x, A, step, step1, t, stats)
        
        [n, m, o] = size(A);
        np = prm.nprm; % number of parameters
        mult = prm.mult; % ? = 1
        ns = n - np * prm.customprm.nC; % number of states
        nv = ns / mult;
        elem = prm.rank_elem;
        
        E = A + repmat(x, 1, m); % restore E
        
        if ~exist('stats', 'var') | isempty(stats)
            stats = struct(                  ...
                'wdir',               {pwd}, ...
                'prm',                {[]},  ...
                'step',               {[]},  ...
                't',                  {[]},  ...
                'rmse_f',             {[]},  ...
                'spread_f',           {[]},  ...
                'rmse_a',             {[]},  ...
                'spread_a',           {[]},  ...
                'rmse_p_f',           {[]},  ...
                'spread_p_f',         {[]},  ...
                'rmse_p_a',           {[]},  ...
                'spread_p_a',         {[]},  ...
                'skewness_f',         {[]},  ...
                'skewness_a',         {[]},  ...
                'kurtosis_f',         {[]},  ...
                'kurtosis_a',         {[]},  ...
                'prm_t',              {[]},  ...
                'prm_f',              {[]},  ...
                'prm_a',              {[]},  ...
                'rank_true_f',        {[]},  ...
                'rank_mean_f',        {[]},  ...
                'rank_true_a',        {[]},  ...
                'rank_mean_a',        {[]},  ...
                'corr_f',             {[]},  ...
                'corr_a',             {[]},  ...
                'var95_f',            {[]},  ...
                'var95_a',            {[]},  ...
                'bestrmse',           {[]},  ...
                'svdspectrum_start',  {[]},  ...
                'svdspectrum_end',    {[]},  ...
                'sum_diff_f',         {[]},  ...
                'sum_diff_sq_f',      {[]},  ...
                'sum_diff_a',         {[]},  ...
                'sum_diff_sq_a',      {[]},  ...
                'lasttime',           {[]},  ...
                'rand_finalstate',    {[]},  ...
                'randn_finalstate',   {[]},  ...
                'randn_obsstate',     {[]},  ...
                'restart',            {[]}   ...
                );
        end
        
        forecast = size(stats.rmse_f, 2) == size(stats.rmse_a, 2);
        analysis = ~forecast;
        
        rms = zeros(mult, 1);
        spread = zeros(mult, 1);
        corr = zeros(mult, 1);
        skew = zeros(mult, 1);
        kurt = zeros(mult, 1);
        var95 = zeros(mult, 1);
        if np > 0                       % np: number of parameters (probably, TS)
            rmse_p = zeros(np, 1);      % length(nC)
            spread_p = zeros(np, 1);    % length(nC)
        end
        
        % note below:   x is ensemble mean
        %               A ensemble anomalies
        
        for mm = 1 : mult
            i1 = (mm - 1) * nv + 1; % first state (TS)
            i2 = mm * nv;           % last state entry (TS)
            % This is how it used to be. Wrong! (Subtracts the sample mean.)
            % rms(mm) = std(x(i1 : i2) - x_true(i1 : i2), 1);
            rms(mm) = rmse(x(i1 : i2) - x_true(i1 : i2));
            spread(mm) = mean(std(A(i1 : i2, :), 1));
            tmp = corrcoef([x(i1 : i2) x_true(i1 : i2)]);
            corr(mm) = tmp(1, 2);
            sumabs = sum(abs(A(i1 : i2, :)'));
            good = find(sumabs > 0) + i1 - 1;
            skew(mm) = sum(abs(skewness(A(good, :)'))) / length(good);
            kurt(mm) = sum(kurtosis(A(good, :)')) / length(good);
            energy = sum(A .* A);
            if isnan(sum(energy))
                if prm.fatalnan
                    error(sprintf('\n  error: A contains NaNs'));
                else
                    var95(mm) = NaN;
                end
            else
                energy = sort(energy);
                energy = energy(end : -1 : 1);
                var95(mm) = min(find(cumsum(energy) >= 0.95 * sum(energy)));
            end
        end
        for par = 1 : np
            rmse_p(par) = rmse(x(ns + par) - x_true(ns + par));
            spread_p(par) = mean(std(A(ns + par, :), 1));
        end
        
        if elem ~= 0
            rank_true = 0;
            rank_mean = 0;
            for i = 1 : m
                if E(elem, i) < x_true(elem)
                    rank_true = rank_true + 1;
                end
                if A(elem, i) < 0
                    rank_mean = rank_mean + 1;
                end
            end
        end
        
        if analysis && prm.calc_bestrmse && isempty(stats.bestrmse)
            if prm.verbose
                fprintf('  calculating bestrmse...');
            end
            if strcmp(prm.model, 'L40') | ns <= m
                bestrmse = 0;
            elseif strcmp(prm.model, 'S')
                bestrmse = calc_bestrmse(x_true(1 : ns), E(1 : ns, :));
            else
                bestrmse = calc_bestrmse(x_true(1 : ns), E(1 : ns, :));
            end
            if bestrmse < 1.0e-10
                bestrmse = 0;
            end
            if prm.verbose
                fprintf('%g\n', bestrmse);
            end
            stats.bestrmse = [bestrmse];
        end
        
        if forecast && isempty(stats.svdspectrum_start)
            [U, S, V] = svd(A(1 : nv, :), 0);
            s = diag(S);
            stats.svdspectrum_start = s / s(1);
        end
        
        %
        % if last assimilation
        %
        step1 = floor(step1 / prm.assim_step) *  prm.assim_step;
        stats.lasttime = analysis & floor((step1 - 1) / prm.assim_step) ...
            == floor((step - 1) / prm.assim_step);
        
        if stats.lasttime
            [U, S, V] = svd(A(1 : nv, :), 0);
            s = diag(S);
            stats.svdspectrum_end = s / s(1);
        end
        
        %
        % store the results
        %
        if analysis
            stats.step = [stats.step step];
            stats.t = [stats.t t];
        end
        
        if forecast
            stats.rmse_f = [stats.rmse_f rms];
            stats.spread_f = [stats.spread_f spread];
            stats.skewness_f = [stats.skewness_f skew];
            stats.kurtosis_f = [stats.kurtosis_f kurt];
            if np > 0
                stats.rmse_p_f = [stats.rmse_p_f rmse_p];
                stats.spread_p_f = [stats.spread_p_f spread_p];
                stats.prm_t = [stats.prm_t x_true(ns + 1 : ns + np)];
                stats.prm_f = [stats.prm_f x(ns + 1 : ns + np)];
            end
            if elem ~= 0
                stats.rank_true_f = [stats.rank_true_f rank_true];
                stats.rank_mean_f = [stats.rank_mean_f rank_mean];
            end
            stats.corr_f = [stats.corr_f corr];
            stats.var95_f = [stats.var95_f var95];
            if isempty(stats.sum_diff_f)
                stats.sum_diff_f = x - x_true;
                stats.sum_diff_sq_f = (x - x_true) .^ 2;
            else
                stats.sum_diff_f = stats.sum_diff_f + x - x_true;
                stats.sum_diff_sq_f = stats.sum_diff_sq_f + (x - x_true) .^ 2;
            end
        else
            stats.rmse_a = [stats.rmse_a rms];
            stats.spread_a = [stats.spread_a spread];
            stats.skewness_a = [stats.skewness_a skew];
            stats.kurtosis_a = [stats.kurtosis_a kurt];
            if np > 0
                stats.rmse_p_a = [stats.rmse_p_a rmse_p];
                stats.spread_p_a = [stats.spread_p_a spread_p];
                stats.prm_a = [stats.prm_a x(ns + 1 : ns + np)];
            end
            if elem ~= 0
                stats.rank_true_a = [stats.rank_true_a rank_true];
                stats.rank_mean_a = [stats.rank_mean_a rank_mean];
            end
            stats.corr_a = [stats.corr_a corr];
            stats.var95_a = [stats.var95_a var95];
            if isempty(stats.sum_diff_a)
                stats.sum_diff_a = x - x_true;
                stats.sum_diff_sq_a = (x - x_true) .^ 2;
            else
                stats.sum_diff_a = stats.sum_diff_a + x - x_true;
                stats.sum_diff_sq_a = stats.sum_diff_sq_a + (x - x_true) .^ 2;
            end
        end
        
        return
    end
    function [prm,pos] = get_pos(prm)
        
        n = prm.n * prm.customprm.nC; % state-par dimension
        np = prm.nprm * prm.customprm.nC; % number of parameters
        ns = (n - np); % number of states
        nv = ns / prm.mult; % prm.mult: number of model variables... since mult=1, nv probably also same meaning
        p = prm.p * prm.customprm.nC; % no. of observations
        
        if prm.nx <= 1 % 1D case - (TS) nx might be the number of spatial dimensions?
            pos = zeros(p, 1);
        else % 2D case
            pos = zeros(p, 2);
            nx = prm.nx;
            ny = nv / prm.nx;
        end
        
        if strcmp(prm.obs_spacing, 'regular')
            if prm.nx <= 1 % 1D case
                for i = 1 : p
                    pos(i) = ceil(nv / p * (i - 0.5));
                end
            else % 2D case
                nn = (nx - 1) * (ny - 1);
                dn = nn / p;
                for i = 1 : p
                    pp = dn * (i - 0.5);
                    pos(i, 1) = mod(pp, nx - 1) + 1;
                    pos(i, 2) = pp / (nx - 1) + 1;
                end
            end
            prm.current.obs = pos;
        elseif strcmp(prm.obs_spacing, 'variable')
            % TS: this is the iMoMo situation where variables can become
            % available or not, depending on the crowd-sourced activity and
            % the avail. of remotely-sensed data. +> temp: cange
            % accordingly for q, G and SN...
            % 10.12.2012, now with snow observtation (including all zero)
            %tt1=prm.current.serialDate - prm.customprm.serialDateStart + 1;
            if prm.customprm.snow % +> not yet adopted to multi-catchment assimilation
                pos = [NaN;prm.customprm.SM(prm.current.step);NaN;...
                    prm.customprm.SN(prm.current.step); prm.customprm.ET(prm.current.step)];
            else
                % has to be dynamically updated given current data
                % availability. +> TODO and adapt to dynamic situation
                % where iMoMo database needs to be queried and new data
                % assimilated in.
                pos = repmat([NaN ; NaN ; NaN ; 1],prm.customprm.nC,1);
            end
            
            idPos = repmat((1:4)',prm.customprm.nC,1);
            addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
            idPos = idPos + addID * prm.n;
            
            pos = pos .* idPos; pos = pos(pos==pos);
            
            prm.current.obs = pos;
            
        else
            error(sprintf('\n  EnKF: error: get_pos(): unknown obs_spacing \"%s\"', prm.obs_spacing));
        end
        
        if strcmp(prm.observe, 'B')
            if prm.mult < 2
                error(sprintf('\n  error: getpos(): prm.observe = B is only allowed for models with more than one variable'));
            else
                pos = pos + nv;
            end
        end
        if strcmp(prm.observe, 'A&B')
            if prm.mult < 2
                error(sprintf('\n  error: getpos(): prm.observe = A&B is only allowed for models with more than one variable'));
            else
                pos = [pos; pos + nv];
            end
        end
        
        if prm.p_inverse
            pos = pos(end : -1 : 1);
        end
        
        return
    end
    function [H] = calc_h(prm, pos)
        
        %keyboard
        
        n = prm.n * prm.customprm.nC;
        np = prm.nprm * prm.customprm.nC;
        mult = prm.mult;
        nv = (n - np) / mult;
        p = size(pos, 1);
        nx = prm.nx;
        
        H = spalloc(p, n, p * 4);
        %
        % '* 4' - to handle bilinear 2D interpolation
        % (should be expanded if a need to handle 2D cases observing A&&B occurs)
        
        if prm.mult == 1 || strcmp(prm.observe, 'A') || strcmp(prm.observe, 'B') || strcmp(prm.observe, 'A&&B')
            for i = 1 : p
                pp = pos(i, :);
                if (nx <= 1) % 1D case
                    if (pp == floor(pp))
                        H(i, pp) = 1;
                    else
                        frac = pp - floor(pp);
                        pp = floor(pp);
                        offset = floor((pp - 1) / nv) * nv; % to handle pp > nv
                        pp1 = mod(pp, nv) + 1 + offset;
                        H(i, pp) = 1 - frac;
                        H(i, pp1) = frac;
                    end
                else % 2D case
                    px = floor(pp(1));
                    py = floor(pp(2));
                    fx = pp(1) - px;
                    fy = pp(2) - py;
                    py = py - 1; % for convenience
                    if fx ~= 0 && fy ~= 0
                        H(i, px + py * nx) = (1 - fx) * (1 - fy);
                        H(i, px + 1 + py * nx) = fx * (1 - fy);
                        H(i, px + (py + 1) * nx) = ( 1 - fx) * fy;
                        H(i, px + 1 + (py + 1) * nx) = fx * fy;
                    elseif fx == 0 && fy == 0
                        H(i, px + py * nx) = 1;
                    elseif fx == 0
                        H(i, px + py * nx) = (1 - fy);
                        H(i, px + (py + 1) * nx) = fy;
                    elseif fy == 0
                        H(i, px + py * nx) = 1 - fx;
                        H(i, px + 1 + py * nx) = fx;
                    end
                end
            end
            
        elseif strcmp(prm.observe, 'A+B')
            
            if (nx > 1) % 1D case
                error(sprintf('\n  EnKF: error: calc_h(): prm.observe = A+B is not implemented for 2D geometry'));
            end
            
            for i = 1 : p
                pp = pos(i);
                if (pp == floor(pp))
                    H(i, pp) = 0.5;
                    H(i, pp + nv) = 0.5;
                else
                    frac = pp - floor(pp);
                    pp = floor(pp);
                    offset = floor((pp - 1) / nv) * nv;  % to handle pp > nv
                    pp1 = mod(pp, nv) + 1 + offset;
                    H(i, pp) = (1 - frac) / 2.0;
                    H(i, pp + nv) = (1 - frac) / 2.0;
                    H(i, pp1) = frac / 2.0;
                    H(i, pp1 + nv) = frac / 2.0;
                end
            end
            
        else
            error(sprintf('\n  EnKF: error: calc_h(): no code to observe \"%s\"', prm.observe));
        end
        
        return;
    end
    function [HE] = model_getHE(prm, H, E)
        
        HE = H * E;
        
        return
    end
    function [dx, A] = assimilateEnKF(prm, A, HA, pos, dy, stats)
        % function [dx, A] = assimilate(prm, A, HA, pos, dy, stats)
        % Copyright (C) 2009 Pavel Sakov (see his assimilate.m)
        
        m = prm.m; % ensemble size
        n = prm.n * prm.customprm.nC; % state dimension
        r = prm.obs_variance(prm.current.obs) * prm.rfactor1;
        rfactor = prm.rfactor2;
        p = size(HA, 1);
        np = prm.nprm * prm.customprm.nC; % number of parameters
        ns = (n - np); % number of states
        if p < 1
            dx = zeros(n, 1);
            return
        end
        
        s = dy ./ sqrt(r * (m - 1)); % s: scaled innovation vector!
        S = HA ./ repmat(sqrt(r * (m - 1)),1,m); %(TS) - old version see orig PS enkf
        
        
        D = randn(p, m) / (rfactor * sqrt(m - 1));
        d = mean(D,2);
        % D is scaled matrix of random pertubations.
        D = sqrt(m / (m - 1)) * (D - repmat(d, 1, m)); % p x m matrix
        
        if prm.loc_len == 0 || isnan(prm.loc_len) % no localisation
            %% COVARIANCE LOCALIZATION - NONE
            if m <= p
                G =  (speye(m) + S' * S) \ S';
            else
                G = S' / (speye(p) + S * S');
            end
            dx = A * G * s;
            if rfactor ~= 1
                S = S / sqrt(rfactor);
                if m <= p
                    G =  (speye(m) + S' * S) \ S';
                else
                    G = S' / (speye(p) + S * S');
                end
            end
            A = A * (speye(m) + G * (D - S));
            
        elseif strcmp(prm.loc_method, 'LA') % local analysis
            %% COVARIANCE LOCALIZATION - Local Analysis
            % Note that we only have to pass A in the arguments because of the local
            % analysis... In practice, with a large-scale system... this can perhaps
            % be avoided by reading, updating, and then writing back slabs from
            % NetCDF data files.
            dx = zeros(n, 1);
            for i = 1 : n
                [localobs, coeffs] = find_localobs(prm, i, pos);
                ploc = length(localobs);
                if ploc == 0
                    continue
                end
                Sloc = S(localobs, :) .* repmat(coeffs, 1, m);
                if m <= ploc
                    Gloc = (speye(m) + Sloc' * Sloc) \ Sloc';
                else
                    Gloc = Sloc' / (speye(ploc) + Sloc * Sloc');
                end
                dx(i) = A(i, :) * Gloc * (s(localobs) .* coeffs);
                if rfactor ~= 1
                    Sloc = Sloc / sqrt(rfactor);
                    if m <= ploc
                        Gloc = (speye(m) + Sloc' * Sloc) \ Sloc';
                    else
                        Gloc = Sloc' / (speye(ploc) + Sloc * Sloc');
                    end
                end
                
                Dloc = D(localobs, :) .* repmat(coeffs, 1, m);
                A(i, :) = A(i, :) + A(i, :) * Gloc * (Dloc - Sloc);
            end
        elseif strcmp(prm.loc_method, 'CF') % covariance filtering
            %% COVARIANCE LOCALIZATION - Covariance Filtering
            K = calc_k(prm, A, HA, r, pos);
            dx = K * dy;
            
            if rfactor ~= 1
                K = calc_k(prm, A, HA, r * rfactor, pos);
            end
            D = randn(p, m) * sqrt(r * rfactor);
            % Subtract the ensemble mean from D to ensure that update of the
            % anomalies does not perturb the ensemble mean. This reduces the
            % variance of each sample by a factor of 1 - 1/m (I think).
            d = mean(D')';
            D = sqrt(m / (m - 1)) * (D - repmat(d, 1, m));
            A = A + K * (D - HA);
            
        else
            error(sprintf('\n  EnKF: error: assimilate(): prm.loc_method = \"%s\" is not suported', prm.loc_method));
        end
        
        idStates = repmat((1:4)',prm.customprm.nC,1); % only no snow
        addID = repmat(1:prm.customprm.nC,prm.n-prm.nprm,1); addID = addID(:) - 1;
        idStates = idStates + addID * prm.n;
        
        idP = repmat((5:8)',prm.customprm.nC,1);
        addID = repmat(1:prm.customprm.nC,prm.n-prm.nprm,1); addID = addID(:) - 1;
        idP = idP + addID * prm.n;
        
        % inflate "normal" (observed) elements
        if prm.inflation ~= 1
            A(idStates, :) = A(idStates, :) * prm.inflation;
        end
        % inflate unobserved elements ("parameters")
        if prm.inflation_prm ~= 1
            A(idP, :) = A(idP, :) * prm.inflation_prm;
        end
        % randomly rotate the ensemble every prm.rotate steps
        if prm.rotate && mod(stats.step(end), prm.rotate) == 0 && prm.rotate_ampl ~= 0
            if prm.rotate_ampl == 1
                A(idStates, :) = A(idStates, :) * genU(m);
            else
                A(idStates, :) = A(idStates, :) * genUeps(m, prm.rotate_ampl);
            end
        end
        return
        
    end % assimilateEnKF
    function [x, t] = model_step(prm, x, t, true_field)
        
        idPStates = repmat((1:4)',prm.customprm.nC,1);
        addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
        idPStates = idPStates + addID * prm.n;
        
        idPar = repmat((5:8)',prm.customprm.nC,1);
        addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
        idPar = idPar + addID * prm.n;
        
        [x(idPStates,:)] = budykoStepFSVec(prm,x(idPStates,:),x(idPar,:),true_field);
        
        t = t + prm.customprm.dt;
        
        return;
        
    end
    function [y] = model_getobs(prm, x, H)
        p = size(H, 1);
        variance = prm.obs_variance(prm.current.obs);
        obsEr = randn(p, 1) .* sqrt(variance);
        y = H * x + obsEr;
        while any(y < 0) % negative observations not permitted!
            y = H * x + randn(p, 1) .* sqrt(variance);
        end
        return
    end
    function [EE] = rescaleE(AA,EE,lb,ub,prm,xold,dx)
        % truncates parameter ensemble members to physically meaningful
        % values.
        lb_EE = EE < lb;
        ub_EE = EE > ub;
        if any(lb_EE(:)) || any(ub_EE(:))
            [x_lb,y_lb] = ind2sub(size(EE),find(lb_EE));
            [x_ub,y_ub] = ind2sub(size(EE),find(ub_EE));
            x_lb_F = x_lb + prm.n-prm.nprm; AAlb = AA(lb_EE); AAlb = AAlb(:);
            x_ub_F = x_ub + prm.n-prm.nprm; AAub = AA(ub_EE); AAub = AAub(:);
            EE(lb_EE) = lb + abs((AAlb + xOld(x_lb_F) - lb) ./ dx(x_lb_F));
            EE(ub_EE) = ub - abs((AAub + xOld(x_ub_F) - ub) ./ dx(x_ub_F));
        end
        return
    end
    function EE = truncateE(EE, prm,lb,ub, step)
        E1 = EE(:) < lb; E2 = EE(:) > ub;
        if any(E1)
            prm.crit1(step) = nansum(E1);
            EE(E1) =  lb;
        end
        if any(E2)
            prm.crit2(step) = nansum(E2);
            EE(E2) =  ub;
        end
    end
    function E = SmaxCheck(E,prm,step) % sets Smax accordingly if problems arise
        
        idPS = (1 : 2 : prm.customprm.nC * 2)';
        
        idPSmax = repmat((prm.nprm+prm.p)',prm.customprm.nC,1);
        addID = repmat(1:prm.customprm.nC,1,1); addID = addID(:) - 1;
        idPSmax = idPSmax + addID * prm.n;
        
        idPSE = repmat(2,prm.customprm.nC,1);
        addID = repmat(1:prm.customprm.nC,1,1); addID = addID(:) - 1;
        idPSE = idPSE + addID * prm.n;
        
        mSCheck = E(idPSE,:);
        
        if prm.customprm.snow % snow
            i = mSCheck > E(end-2,:); % note, Smax is now stored in the end-2 slot!
            prm.crit3(step) = nansum(i);
            if any(i)
                E(end-2,i) = mSCheck(i);
            end
        else % no snow
            Smax = E(idPSmax,:);
            id2corr = mSCheck>Smax;
            mSCheck(id2corr) = Smax(id2corr);
            E(idPSE,:) = mSCheck;
            prm.crit3(step) = nansum(i);
        end
    end
    function [K] = calc_k(prm, A, HA, r, pos)
        % Pavel Sakov, 2008
        [n, m] = size(A);
        p = size(HA, 1);
        
        np = prm.nprm;
        ns = n - np;
        mult = prm.mult;
        nv = ns / mult;
        nx = prm.nx;
        loclen = prm.loc_len;
        
        HPHT = HA * HA' / (m - 1);
        PHT = A * HA' / (m - 1);
        R = speye(size(HA, 1)) * r;
        
        if loclen == 0 || isnan(loclen) % no localisation
            
            K = PHT / (HPHT + R);
            
        else % localisation
            
            % masking HPH^T
            %
            if p > 1
                S = zeros(p, p);
                if nx <= 1 % 1D case
                    pos = mod(pos - 1, nv) + 1; % to handle prm.observe = 'A&B' and 'B'
                    for o = 1 : p
                        xo = pos(o);
                        dist = abs(pos(o : end) - xo);
                        if prm.periodic
                            dist = min(dist, nv - dist);
                        end
                        coeffs = calc_loccoeffs(loclen, prm.loc_function, dist);
                        S(o, o : end) = coeffs;
                        S(o, 1 : o - 1) = S(1 : o - 1, o);
                    end
                else  % 2D case
                    for o = 1 : p
                        xo = pos(o, 1);
                        yo = pos(o, 2);
                        dist = hypot(pos(o : end, 1) - xo, pos(o : end, 2) - yo);
                        coeffs = calc_loccoeffs(loclen, prm.loc_function, dist);
                        S(o, o : end) = coeffs;
                        S(o, 1 : o - 1) = S(1 : o - 1, o);
                    end
                end
                HPHT = S .* HPHT;
            end % if p > 1
            
            % now masking PH^T
            %
            if nx <= 1 % 1D case
                grid = 1 : nv;
                if ~prm.periodic
                    dist = grid - 1;
                else
                    dist = min(grid - 1, nv + 1 - grid);
                end
                coeffs = calc_loccoeffs(loclen, prm.loc_function, dist);
                
                v = zeros(nv, 1);
                for o = 1 : p
                    xo = round(pos(o));
                    offset = 0;
                    for k = 1 : mult
                        v(xo + offset : nv + offset) = coeffs(1 : nv + 1 - xo);
                        v(1 + offset : xo - 1 + offset) = coeffs(nv + 2 - xo : nv);
                        offset = offset + nv;
                    end
                    PHT(:, o) = PHT(:, o) .* v;
                end
            else % 2D case
                ny = nv / nx;
                if floor(ny) ~= nx
                    error(sprintf('\n  EnKF: error: calc_k(): only square domains are currently handled in 2D cases\n'));
                end
                
                I = zeros(nx, ny);
                for i = 1 : nx
                    I(i, :) = i;
                end
                J = zeros(nx, ny);
                for j = 1 : ny
                    J(:, j) = j;
                end
                
                for o = 1 : p
                    xo = pos(o, 1);
                    yo = pos(o, 2);
                    D = hypot(I - xo, J - yo);
                    dist = reshape(D, nv, 1);
                    coeffs = calc_loccoeffs(loclen, prm.loc_function, dist);
                    PHT(:, o) = PHT(:, o) .* coeffs;
                end
            end
            
            K = PHT / (HPHT + R);
        end
        return
    end
end
