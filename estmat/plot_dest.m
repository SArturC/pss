function output = plot_dest(A,B,C,D,varargin)

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
    if ~isfield(options,'d')
        options.d = 0.5;
    end
    if ~isfield(options,'r')
        options.r = 5;
    end
    if ~isfield(options,'theta')
        options.theta = pi/4;
    end

    h     = options.grid_step;
    d     = options.d;
    r     = options.r;
    theta = options.theta;

    Nv = length(A); % número de vértices

    poles_real = [];
    poles_imag = [];

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

    %% ======== LOOP =========
    n_samples = size(Alpha,1);

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
        p = pole(H);

        poles_real = [poles_real; real(p)'];
        poles_imag = [poles_imag; imag(p)'];
    end

    %% ===================== PLOT =====================
    figure
    hold on
    grid on

    cmap1 = jet(n_samples);

    for i = 1:size(poles_real,2)
        for k = 1:n_samples
            scatter(poles_real(k,i),poles_imag(k,i),30,cmap1(k,:),'filled')
        end
    end

    %% ===================== D-STABILITY =====================

    % linha vertical
    xline(-d,'k--','LineWidth',2)

    % círculo (semiplano esquerdo)
    t = linspace(-pi/2,pi/2,500);
    plot(-r*cos(t), r*sin(t),'k--','LineWidth',2)

    % cone
    x = linspace(-r,0,500);
    y1 = tan(theta)*(-x);
    y2 = -tan(theta)*(-x);

    plot(x,y1,'k--','LineWidth',2)
    plot(x,y2,'k--','LineWidth',2)

    %% estética
    xlabel('Re\{p\}')
    ylabel('Im\{p\}')
    title('Polos + Região de D-estabilidade')

    legend('Polos','Re(s)=-d','|s|=r','Cone \theta')

    axis equal
%     xlim([-r 1])
%     ylim([-r r])

    output = [];

end