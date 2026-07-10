function output = nuvem(A,type,d,r,q,theta)
    if nargin < 2
        type = 'continuous';
    end
    plotD = (nargin >= 6);
    
    Nv = length(A);
    Alpha = [];
    Nrandom = 1000;
    
    for s = 1:Nrandom
        alpha = simplex_random(Nv);
        Alpha = [Alpha ; alpha];
    end
    
    Nedge = 100;
    
    for i = 1:Nv-1
        for j = i+1:Nv
            t = linspace(0,1,Nedge);
            for k = 1:length(t)
                alpha = zeros(1,Nv);
    
                alpha(i) = t(k);
                alpha(j) = 1-t(k);
    
                Alpha = [Alpha ; alpha];
            end
        end
    end
    
    if Nv >= 4
        Nface = 150;
        triples = nchoosek(1:Nv,3);
        for p = 1:size(triples,1)
            i = triples(p,1);
            j = triples(p,2);
            k = triples(p,3);
            for s = 1:Nface
                beta = simplex_random(3);
    
                alpha = zeros(1,Nv);
    
                alpha(i) = beta(1);
                alpha(j) = beta(2);
                alpha(k) = beta(3);
    
                Alpha = [Alpha ; alpha];
            end
        end
    end
    
    all_eigs    = []  ;
    worst_value = -Inf;
    worst_alpha = []  ;
    
    for p = 1:size(Alpha,1)
        alpha = Alpha(p,:);
        Aalpha = 0;
    
        for i = 1:Nv
            Aalpha = Aalpha + alpha(i)*A{i};
        end
    
        lambda = eig(Aalpha);
        all_eigs = [all_eigs ; lambda];
    
        switch lower(type)
            case 'continuous'
                metric = max(real(lambda));
            case 'discrete'
                metric = max(abs(lambda));
            otherwise
                error('Tipo deve ser continuous ou discrete')
        end
    
        if metric > worst_value
            worst_value = metric;
            worst_alpha = alpha;
        end
    end
    
    figure
    hold on
    grid on
    box on

    %----------------------------------------------------------
    % Região D
    %----------------------------------------------------------
    if plotD && strcmpi(type,'continuous')

        % reta vertical
        xline(-d,'k--','LineWidth',1);

        % círculo
        t = linspace(0,2*pi,500);

        plot(q + r*cos(t), ...
            r*sin(t), ...
            'k--', ...
            'LineWidth',1);

        % cone
        xmax = max(real(all_eigs)) + 0.5;
        xmin = min(real(all_eigs)) - 0.5;

        L = max(abs([xmin xmax])) + r + abs(q);
        x = linspace(-L,0,300);
        m = tan(theta);

        plot(x,  m*x,'k--','LineWidth',1)
        plot(x, -m*x,'k--','LineWidth',1)
    end
    
    plot(real(all_eigs),imag(all_eigs),'.','MarkerSize',10)
    
    vertex_eigs = [];
    for i = 1:Nv
        vertex_eigs = [vertex_eigs ; eig(A{i})];
    end
    plot(real(vertex_eigs), ...
        imag(vertex_eigs), ...
        'ro', ...
        'MarkerSize',12, ...
        'LineWidth',2.5, ...
        'MarkerFaceColor','none')
    
    xline(0,'k'); yline(0,'k')
    if strcmpi(type,'discrete')
        th = linspace(0,2*pi,1000);
        plot(cos(th),sin(th),'r','LineWidth',2)
    end
    xlabel('Real'); ylabel('Imag')
    if strcmpi(type,'continuous')
        title('Autovalores do Politopo')
    else
        title('Autovalores do Politopo e Círculo Unitário')
    end
    axis equal

    output.worst_value = worst_value;
    output.alpha       = worst_alpha;
    if strcmpi(type,'continuous')
        output.stable = (worst_value < 0);
    else
        output.stable = (worst_value < 1);
    end
end

function x = simplex_random(N)
    x = zeros(1,N);
    x(1) = 1-rand^(1/(N-1));
    for k = 2:N-1
        x(k) = (1-sum(x(1:k-1)))*(1-rand^(1/(N-k)));
    end
    x(N) = 1-sum(x(1:N-1));
end