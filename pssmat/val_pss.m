function results = val_pss( ...
                    A_poli, Bu_poli, Bw_poli, ...
                    C_poli, Du_poli, Dw_poli, ...
                    K, w_grid, Nstart)

% function results = val_pss(...)
%
% Verifica passividade de um sistema politópico em malha fechada via
% detecção de contraexemplo no simplex dos parâmetros convexos.
%
% Estratégia eficiente:
%   1) Teste dos vértices
%   2) Amostragem grosseira no simplex
%   3) Refinamento local (fmincon) apenas quando necessário
%
% Critério:
%   Existe alpha tal que:
%
%       lambda_min((G(jw,alpha)+G(jw,alpha)*)/2) < 0
%
% => sistema NÃO passivo
%
% Inputs:
%   A_poli, Bu_poli, Bw_poli, C_poli, Du_poli, Dw_poli -> vértices
%   K        -> ganho
%   w_grid   -> grade de frequência (default: leve)
%   Nstart   -> número de refinamentos locais
%
% Outputs:
%   results.all_passive
%   results.counterexample
%   results.w
%   results.minval

    v = length(A_poli);

    % ---------------- DEFAULTS ----------------
    if nargin < 8 || isempty(w_grid)
        w_grid = logspace(-2,4,40); % mais leve
    end

    if nargin < 9 || isempty(Nstart)
        Nstart = 2*v; % bem menor que antes
    end

    Nsample = 20*v;           % amostragem leve
    threshold_refine = 1e-3;  % ativa refinamento

    results.all_passive = true;
    results.counterexample = [];
    results.w = [];
    results.minval = [];

    % ============================================================
    % 1) TESTE DOS VÉRTICES
    % ============================================================
    for i = 1:v

        A  = A_poli{i};
        Bu = Bu_poli{i};
        Bw = Bw_poli{i};
        C  = C_poli{i};
        Du = Du_poli{i};
        Dw = Dw_poli{i};

        Acl = A + Bu*K;
        Ccl = C + Du*K;

        sys = ss(Acl,Bw,Ccl,Dw);

        if ~isPassive(sys)
            results.all_passive = false;
            results.counterexample = sprintf('vertex %d',i);
            return
        end
    end

    % ============================================================
    % LOOP EM FREQUÊNCIA
    % ============================================================
    for w = w_grid

        % ========================================================
        % 2) AMOSTRAGEM NO SIMPLEX
        % ========================================================
        Alpha = rand(v, Nsample);
        Alpha = Alpha ./ sum(Alpha,1);

        vals = zeros(1,Nsample);

        for k = 1:Nsample
            vals(k) = cost_alpha(Alpha(:,k), ...
                A_poli, Bu_poli, Bw_poli, ...
                C_poli, Du_poli, Dw_poli, ...
                K, w);
        end

        [minval, idx] = min(vals);
        alpha_best = Alpha(:,idx);

        % --------------------------------------------------------
        % EARLY STOP
        % --------------------------------------------------------
        if minval < 0
            results.all_passive = false;
            results.counterexample = alpha_best;
            results.w = w;
            results.minval = minval;
            return
        end

        % ========================================================
        % 3) REFINAMENTO LOCAL
        % ========================================================
        if minval < threshold_refine

            fun = @(alpha) cost_alpha(alpha, ...
                A_poli, Bu_poli, Bw_poli, ...
                C_poli, Du_poli, Dw_poli, ...
                K, w);

            opts = optimoptions('fmincon', ...
                'Display','off', ...
                'Algorithm','sqp', ...
                'MaxIterations',50, ...
                'OptimalityTolerance',1e-5);

            Aeq = ones(1,v);
            beq = 1;
            lb = zeros(v,1);
            ub = ones(v,1);

            % multi-start leve
            for s = 1:Nstart

                if s == 1
                    alpha0 = alpha_best;
                else
                    alpha0 = rand(v,1);
                    alpha0 = alpha0 / sum(alpha0);
                end

                try
                    alpha_opt = fmincon(fun, alpha0, [], [], ...
                        Aeq, beq, lb, ub, [], opts);

                    val = fun(alpha_opt);

                    if val < 0
                        results.all_passive = false;
                        results.counterexample = alpha_opt;
                        results.w = w;
                        results.minval = val;
                        return
                    end

                catch
                    continue
                end
            end
        end

    end

end

% ================================================================
% FUNÇÃO DE CUSTO
% ================================================================
function val = cost_alpha(alpha, ...
    A_poli, Bu_poli, Bw_poli, ...
    C_poli, Du_poli, Dw_poli, ...
    K, w)

    alpha = alpha / sum(alpha);

    v = length(alpha);

    A=0; Bu=0; Bw=0; C=0; Du=0; Dw=0;

    for i=1:v
        A  = A  + alpha(i)*A_poli{i};
        Bu = Bu + alpha(i)*Bu_poli{i};
        Bw = Bw + alpha(i)*Bw_poli{i};
        C  = C  + alpha(i)*C_poli{i};
        Du = Du + alpha(i)*Du_poli{i};
        Dw = Dw + alpha(i)*Dw_poli{i};
    end

    Acl = A + Bu*K;
    Ccl = C + Du*K;

    sys = ss(Acl,Bw,Ccl,Dw);

    Gjw = squeeze(freqresp(sys,w));

    M = (Gjw + Gjw')/2;

    val = min(eig(M));

end
