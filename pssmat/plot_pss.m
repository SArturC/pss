function output = plot_pss(A,B,C,D,varargin)

    options = [];

    if nargin > 4
        if nargin == 5
            options = varargin{1};
        else
            options = struct(varargin{:});
        end
    end

    % parâmetros padrão
    if ~isfield(options,'grid_step')
        options.grid_step = 0.1;
    end
    if ~isfield(options,'frequency')
        options.frequency = logspace(-1,6,1000);
    end

    h = options.grid_step;
    w = options.frequency;

    Nv = length(A);

    % garantir formato cell para B,C,D
    if ~iscell(B), B = repmat({B},1,Nv); end
    if ~iscell(C), C = repmat({C},1,Nv); end
    if ~iscell(D), D = repmat({D},1,Nv); end

    %% ======== AMOSTRAGEM DO SIMPLEX =========
    levels = 0:h:1;
    M = length(levels);

    % gera todas as combinações possíveis
    grids = cell(1,Nv);
    [grids{:}] = ndgrid(levels);

    Alpha_full = zeros(numel(grids{1}), Nv);

    for i = 1:Nv
        Alpha_full(:,i) = grids{i}(:);
    end

    % filtrar apenas os pontos que satisfazem soma = 1
    tol = 1e-6;
    idx = abs(sum(Alpha_full,2) - 1) < tol;

    Alpha = Alpha_full(idx,:);

    poles_real = [];
    poles_imag = [];

    passive_flag = true;

    %% ======== LOOP =========
    n_samples = size(Alpha,1);
    passivity  = zeros(n_samples,1);

    for k = 1:n_samples

        alpha = Alpha(k,:);

        % combinação convexa
        Aalpha = 0;
        Balpha = 0;
        Calpha = 0;
        Dalpha = 0;

        for i = 1:Nv
            Aalpha = Aalpha + alpha(i)*A{i};
            Balpha = Balpha + alpha(i)*B{i};
            Calpha = Calpha + alpha(i)*C{i};
            Dalpha = Dalpha + alpha(i)*D{i};
        end

        H = ss(Aalpha,Balpha,Calpha,Dalpha);

        %% ======== PASSIVIDADE MIMO =========
        Gjw = freqresp(H,w);

        min_lambda = +inf;

        for iw = 1:length(w)

            G = Gjw(:,:,iw);

            M = G + G'; % G(jw) + G(jw)^*

            lambda_min = min(eig(M));

            if lambda_min < min_lambda
                min_lambda = lambda_min;
            end
        end

        passivity(k) = min_lambda;

        if min_lambda < 0
            passive_flag = false;
        end
    end

    output.pss = passive_flag;

    %% ===================== PLOT 2: PASSIVIDADE =====================
    figure
    plot(passivity,'LineWidth',2)
    grid on

    xlabel('amostra no politopo')
    ylabel('min \lambda_{min}(G(j\omega)+G(j\omega)^*)')

    yline(0,'k--')

    title('Passividade (critério geral MIMO)')

end