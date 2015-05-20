function [status, x_true2, E2, stats] = ienkf_cycle(prm, x_true1, E1, cycle)

    nitermax = prm.nitermax; %maximum interation
    cycle_nstep = prm.cycle_nstep; %time steps per assimilation cycle
    obs_nstep = prm.obs_nstep; %time step when observation is available and assimilation is peformed
    infl = prm.inflation; %inflation parameter
    propagate = prm.propagate; %model propagation
    observe = prm.observe; %function to produce observation error from x_true and stdv
    r = sqrt(prm.obs_variance); %standard deviation of observations

    method = prm.method;
    if strcmp(method, 'IEnKF') % Not used
    elseif strcmp(method, 'IEKF')
        if ~isfield(prm, 'C') 
            C = 10^4;
        else
            C = prm.C;
        end
    elseif strcmp(method, 'EnRML')
    elseif strcmp(method, 'EnKF')
        nitermax = 1; %only one iteration with EnKF
    else
        error(sprintf('\n  error: ienkf_step(): input: prm.obs_method = \"%s\" not known', method));
    end
    
    [n, m] = size(E1); % n=number of states; m=number of ensembles
   
    calc_stats = nargout == 4; %output calculation stats?
    
    if calc_stats %initialize stats
        stats.niter = 0;
        stats.rmse_f = NaN;
        stats.rmse_a = NaN;
        stats.spread_f = NaN;
        stats.spread_a = NaN;
        stats.kurt_f = NaN;
        stats.kurt_a = NaN;
        stats.rmse = zeros(1, nitermax);
        stats.inc = zeros(1, nitermax);
        stats.dy = zeros(1, nitermax);
        if prm.calc_norms
            stats.normA1_f = NaN;
            stats.normA2_f = NaN;
            stats.normT = NaN;
        end
    end
    
    xf = mean(E1')'; % Ensemble Mean
    Af = E1 - repmat(xf, 1, m); % Deviation from Ensemble Mean
    
    y = [];
    pos = [];
    x1 = xf; %Ensemble Mean
    A1 = Af; %Deviation from Ensemble mean
    x_true = x_true1; %true states
    
    if  strcmp(method, 'IEnKF') | strcmp(method, 'EnRML') %not used
        T = speye(m);
    end
    if nitermax > 1 %not used
        AA = pinv_ps(Af' * Af);
    end
    
    iter = 1;
    while 1

        % scale the ensemble anomalies %not used
        %
        if strcmp(method, 'IEnKF') | strcmp(method, 'EnRML')
            A1 = Af * T;
        elseif strcmp(method, 'IEKF')
            A1 = Af / C;
        end
        
        % propagate the true field (once), collect observations
        %
        if iter == 1
            for step = 1 : cycle_nstep
                [x_true] = propagate(x_true); %propagate model from true field
                if mod(step, obs_nstep) == 0 %if assimilation step
                    y_now = observe(x_true, r); %produce observation from true field and stdv
                    y = [y; y_now];
                end
            end
            x_true2 = x_true; %new true field for next assimilation step
        end
            
        % propagate the ensemble, calculate ensemble observations
        %
        E2 = repmat(x1, 1, m) + A1; %E2 = E1, Ensemble Mean + Deviations
        HE = [];
        for step = 1 : cycle_nstep %each time step
            for e = 1 : m %each ensemble
                E2(:, e) = propagate(E2(:, e)); %propagate Model from previos model state
            end
            if mod(step, obs_nstep) == 0 %if assimilation step
                for e = 1 : m %each ensemble
                    HE_now(:, e) = observe(E2(:, e), 0); %HE_now = E2 = model states of each ensemble at assimilation step
                end
                HE = [HE; HE_now];
            end
        end
        % calculate ensemble observation anomalies
        %
        Hx = mean(HE')'; % Ensemble Mean
        HA = HE - repmat(Hx, 1, m); % Deviation of Ensembles from Ensemble Mean
                            
        dy = y - Hx; %Deviation of Ensemble mean from observation
        
        p = size(HE, 1); % number of states

        % rescale the ensemble anomalies or use regression
        % not used
        if strcmp(method, 'IEnKF') & iter > 1
            if prm.epsilon == 0 
                HA = HA * inv(T);
            else
                HA = HA * U * diag(1 ./ lsqrt) * U';
            end
        elseif strcmp(method, 'IEKF')
            HA = HA * C;
        elseif strcmp(method, 'EnRML')
            HM = HA * pinv_ps(A1); % ("Gl" is used for "HM"
                                   % in Gu and Oliver 2007)

            HA = HM * Af; % "would-be" propagated ensemble anomalies if the
                          % forecast anomalies were propagated with improved
                          % gradients (M and H)
        end
        %HA = HA * infl;

        % calculate standardised innovation and standardised ensemble anomalies
        %
        s = dy / (r * sqrt(m - 1)); %Deviation of Ensemble mean from observation
        S = HA / (r * sqrt(m - 1)); %Deviation of Ensembles from Ensemble Mean
        
        if prm.epsilon == 0
            if m <= p
                if prm.use_cholesky
                    GGc = chol(speye(m) + S' * S);
                    G = GGc \ (GGc' \ speye(m));
                else
                    G = inv(speye(m) + S' * S);
                end
            else
                if prm.use_cholesky
                    GGc = chol(speye(p) + S * S');
                    G = eye(m) - S' * (GGc \ (GGc' \ speye(p))) * S;
                else
                    G = eye(m) - S' * inv(speye(p) + S * S') * S;
                end
            end
        else
            [U, L] = svd(speye(m) + S' * S);
            l = diag(L);
            % note: U and l will be reused for calculating T and inv(T)
            G = U * diag(1 ./ l) * U';
        end
        if ~isreal(G) | ~isfinite(G)
            status = 2; % fatal error (instability)
            return
        end

        %Not used
        if 1 == 0 and strcmp(method, 'EnRML') % an alternative for the EnRML 
            
            dx = Af * G * S' * (dy + HM * (x1 - xf)) / (r * sqrt(m - 1));
            dx = xf + dx - x1; % calculate dx relative to the prevous
                               % iteration - to sync with the other methods
        else
            b = G * S' * s;
        
            if ~strcmp(method, 'EnKF') | calc_stats
                if iter == 1
                    dx = Af * b;
                else
                    dx = Af * b + Af * (G * (AA * (Af' * (xf - x1))));
                end
            end
        end

        if calc_stats
            stats.inc(iter) = rmse(dx);
            stats.dy(iter) = rmse(dy);
            x2 = mean(E2')'; % Ensemble mean
            if iter == 1 % calculate forecast statics at t2
                A2 = E2 - repmat(x2, 1, m);
                stats.rmse_f = rmse(x2 - x_true);
                stats.kurt_f = mean(kurtosis(A2'));
                if strcmp(method, 'IEKF')
                    A2 = A2 * C;
                end
                stats.spread_f = mean(rmse(A2'));
                if stats.spread_f / stats.rmse_f < prm.collapseratio
                    status = 3; % fatal error (ensemble collapse)
                    return
                end
            end
            stats.rmse(iter) = rmse(x2 - x_true);
        end
        
        % exit
        %
        if rmse(dx) < r / 1000 | iter == nitermax
            x2 = mean(E2')'; % Ensemble mean
            A2 = E2 - repmat(x2, 1, m);% Deviation from Ensemble Mean

            if strcmp(method, 'IEnKF') | strcmp(method, 'EnRML')
                A2 = A2 * inv(T);
            elseif strcmp(method, 'IEKF')
                A2 = A2 * C;
            end
            
            T = sqrtm(G);
            if prm.calc_norms
                stats.normA1_f = norm(Af);
                stats.normA2_f = norm(A2);
                stats.normT = norm(inv(T));
                stats.normHA_f = norm(HA);
                stats.normHA_a = norm(HA * T);
            end
            
            if strcmp(method, 'EnKF')
                x2 = x2 + A2 * b;
            end
            
            A2 = A2 * (T * infl);
            if prm.rotate
                A2 = A2 * genU(m);
            end
            %% update model States
            E2 = repmat(x2, 1, m) + A2;
            %%

            if calc_stats % calculate analysis statics at t2
                stats.rmse_a = rmse(x2 - x_true);
                stats.spread_a = mean(rmse(A2'));
                stats.kurt_a = mean(kurtosis(A2'));
                stats.niter = iter;
            end
            
            if iter <= nitermax
                status = 0; % success
            else
                status = 1; % non-fatal error (no convergence)
            end
            
            return
        end
        
        x1 = x1 + dx;
        if strcmp(method, 'IEnKF') | strcmp(method, 'EnRML')
            if prm.epsilon == 0
                T = sqrtm(G);
            else
                lsqrt = sqrt(sqrt(1 ./ l.^ 2 + prm.epsilon));
                T = U * diag(lsqrt) * U';
            end
        end
        iter = iter + 1;
    end
    
    return
    