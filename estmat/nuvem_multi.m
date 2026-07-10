function output = nuvem_multi(A,varargin)
    options = [];
    output = [];
    if nargin > 1
        if nargin == 2
            options = varargin{1};
        else
            options = struct(varargin{:});
        end
    end

    if ~isfield(options,'type')
        options.type = 'continuous';
    end
    if ~isfield(options,'d')
        options.d = 0;
    end
    if ~isfield(options,'r')
        options.r = 0;
    end
    if ~isfield(options,'q')
        options.q = 0;
    end
    if ~isfield(options,'theta')
        options.theta = 0;
    end
    if ~isfield(options,'vertices')
        options.vertices = length(A);
    end
    if ~isfield(options,'degrees')
        options.degrees = 1;
    end
    if ~isfield(options,'save_csv')
        options.save_csv = false;
    end
    
    if ~isfield(options,'csv_file')
        options.csv_file = "";
    end

    plotD = isfield(options,'d') && ...
            isfield(options,'r') && ...
            isfield(options,'q') && ...
            isfield(options,'theta');

    VertexPoints = generateVertices(options);
    CloudPoints  = mergeCloudPoints(generateRandomPoints(options,1000)            , ...
                                    generateEdges       (options,100,VertexPoints), ...
                                    generateFaces       (options,150,VertexPoints)     );

    CloudData  = evaluateEigenCloud(A,CloudPoints ,options);
    VertexData = evaluateEigenCloud(A,VertexPoints,options);
    
    all_eigs    = CloudData.eigs;
    vertex_eigs = VertexData.eigs;

    if options.save_csv
        csvDataCloud = table(real(CloudData.eigs), ...
                             imag(CloudData.eigs), ...
                             'VariableNames',{'Real','Imag'});
        writetable(csvDataCloud,options.csv_file);

        csvVertex = strrep(options.csv_file,'.csv','_vertices.csv');
        csvDataVertex = table(real(VertexData.eigs), ...
                              imag(VertexData.eigs), ...
                              'VariableNames',{'Real','Imag'});
        writetable(csvDataVertex,csvVertex);
    end
    
    worst_value = CloudData.worst_value;
    worst_alpha = CloudData.worst_alpha;
    
    figure
    hold on
    grid on
    box on

    plot(real(all_eigs),imag(all_eigs),'.','MarkerSize',10)

    plot(real(vertex_eigs), ...
        imag(vertex_eigs), ...
        'ro', ...
        'MarkerSize',12, ...
        'LineWidth',2.5, ...
        'MarkerFaceColor','none')
    
    xline(0,'k'); yline(0,'k')
    if strcmpi(options.type,'discrete')
        th = linspace(0,2*pi,1000);
        plot(cos(th),sin(th),'r','LineWidth',2)
    end
    xlabel('Real'); ylabel('Imag')
    if strcmpi(options.type,'continuous')
        title('Autovalores do Politopo')
    else
        title('Autovalores do Politopo e Círculo Unitário')
    end
    axis equal

    %----------------------------------------------------------
    % Região D
    %----------------------------------------------------------
    if plotD && strcmpi(options.type,'continuous')

        % reta vertical
        xline(-options.d,'k--','LineWidth',1);

        % círculo
        t = linspace(0,2*pi,500);

        plot(options.q + options.r*cos(t), ...
            options.r*sin(t), ...
            'k--', ...
            'LineWidth',1);

        % cone
        xmax = max(real(all_eigs)) + 0.5;
        xmin = min(real(all_eigs)) - 0.5;

        L = max(abs([xmin xmax])) + options.r + abs(options.q);
        x = linspace(-L,0,300);
        m = tan(options.theta);

        plot(x,  m*x,'k--','LineWidth',1)
        plot(x, -m*x,'k--','LineWidth',1)
    end

    output.worst_value = worst_value;
    output.alpha       = worst_alpha;
    if strcmpi(options.type,'continuous')
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

function mono = homogeneousBasis(alpha, degree)
    %HOMOGENEOUSBASIS Avalia todos os monômios homogêneos de um dado grau.
    %
    % mono = homogeneousBasis(alpha, degree)
    %
    % alpha  : coordenadas baricêntricas (linha ou coluna)
    % degree : grau do polinômio
    %
    % Exemplos:
    %
    % alpha = [a1 a2];
    %
    % degree = 1
    % mono =
    %   [a1
    %    a2]
    %
    % degree = 2
    % mono =
    %   [a1^2
    %    a1*a2
    %    a2^2]
    %
    % degree = 3
    % mono =
    %   [a1^3
    %    a1^2*a2
    %    a1*a2^2
    %    a2^3]

    alpha = alpha(:).';
    n = length(alpha);
    % Todos os expoentes possíveis
    expo = compositions(degree,n);
    mono = zeros(size(expo,1),1);
    for k = 1:size(expo,1)
        mono(k) = prod(alpha.^expo(k,:));
    end
end


function E = compositions(total,n)
% Todas as composições de 'total' em 'n' partes.
    if n == 1
        E = total;
        return
    end
    E = [];
    for k = total:-1:0
        T = compositions(total-k,n-1);
        E = [E;
             [k*ones(size(T,1),1) T]];
    end
end

function Aalpha = evaluateMultiSimplex(A,alpha,options)
    mono = cell(length(alpha),1);
    for s = 1:length(alpha)
        mono{s} = homogeneousBasis(alpha{s},options.degrees(s));
    end
    coef = mono{1};
    for s = 2:numel(mono)
        coef = kron(mono{s},coef);
    end
    Aalpha = zeros(size(A{1}));
    for k = 1:numel(coef)
        Aalpha = Aalpha + coef(k)*A{k};
    end
end

function CloudPoints = generateRandomPoints(options,Nrandom)
    m = length(options.vertices);
    CloudPoints = cell(m,1);
    for s = 1:m
        CloudPoints{s} = zeros(Nrandom,options.vertices(s));
        for k = 1:Nrandom
            CloudPoints{s}(k,:) = simplex_random(options.vertices(s));
        end
    end
end

function CloudPoints = generateVertices(options)
    m = length(options.vertices);
    sizes = options.vertices;
    N = prod(sizes);
    CloudPoints = cell(m,1);
    for s = 1:m
        CloudPoints{s} = zeros(N,sizes(s));
    end
    for k = 1:N
        idx = cell(1,m);
        [idx{:}] = ind2sub(sizes,k);
        for s = 1:m
            CloudPoints{s}(k,idx{s}) = 1;
        end
    end
end

function CloudPoints = generateEdges(options,Nedge,VertexPoints)
    m = numel(options.vertices);
    CloudPoints = cell(m,1);
    for s = 1:m
        CloudPoints{s} = zeros(0,options.vertices(s));
    end
    
    for free = 1:m
        Nv = options.vertices(free);
        pairs = nchoosek(1:Nv,2);
    
        for e = 1:size(pairs,1)
            i = pairs(e,1);
            j = pairs(e,2);
            t = linspace(0,1,Nedge);
    
            for p = 1:size(VertexPoints{1},1)
                alpha = getPoint(VertexPoints,p);
    
                for k = 1:length(t)
                    alpha_edge = alpha;
                    alpha_edge{free} = zeros(1,Nv);
                    alpha_edge{free}(i) = t(k);
                    alpha_edge{free}(j) = 1-t(k);
                    CloudPoints = appendPoint(CloudPoints,alpha_edge);
                end
            end
        end
    end
end

function CloudPoints = generateFaces(options,Nface,VertexPoints)
    m = numel(options.vertices);
    CloudPoints = cell(m,1);
    for s = 1:m
        CloudPoints{s} = zeros(0,options.vertices(s));
    end
    
    freePairs = nchoosek(1:m,2);

    for f = 1:size(freePairs,1)
        free1 = freePairs(f,1);
        free2 = freePairs(f,2);
    
        for p = 1:size(VertexPoints{1},1)
            alpha = getPoint(VertexPoints,p);
    
            for k = 1:Nface
                alpha_face = alpha;
                alpha_face{free1} = simplex_random(options.vertices(free1));
                alpha_face{free2} = simplex_random(options.vertices(free2));
                CloudPoints = appendPoint(CloudPoints,alpha_face);
            end
        end
    end
end

function CloudPoints = mergeCloudPoints(varargin)
    if nargin==0
        CloudPoints = {};
        return
    end
    m = length(varargin{1});
    CloudPoints = cell(m,1);
    for s = 1:m
        CloudPoints{s} = [];
    end
    for k = 1:nargin
        A = varargin{k};
        for s = 1:m
            CloudPoints{s} = [CloudPoints{s}; A{s}];
        end
    end
end

function out = evaluateEigenCloud(A,CloudPoints,options)
    all_eigs    = [];
    worst_value = -Inf;
    worst_alpha = [];
    
    for p = 1:size(CloudPoints{1},1)
        alpha = getPoint(CloudPoints,p);
    
        for s = 1:length(CloudPoints)
            alpha{s} = CloudPoints{s}(p,:);
        end
    
        Aalpha = evaluateMultiSimplex(A,alpha,options);
        lambda = eig(Aalpha);
        all_eigs = [all_eigs ; lambda];
    
        switch lower(options.type)
            case 'continuous'
                metric = max(real(lambda));
            case 'discrete'
                metric = max(abs(lambda));
            otherwise
                error('Tipo deve ser continuous ou discrete')
        end
    
        if metric>worst_value
            worst_value = metric;
            worst_alpha = alpha;
        end
    end
    
    out.eigs        = all_eigs;
    out.worst_value = worst_value;
    out.worst_alpha = worst_alpha;
end

function alpha = getPoint(CloudPoints,p)
    m = numel(CloudPoints);
    alpha = cell(m,1);
    
    for s = 1:m
        alpha{s} = CloudPoints{s}(p,:);
    end
end

function CloudPoints = appendPoint(CloudPoints,alpha)
    for s = 1:numel(alpha)
        CloudPoints{s}(end+1,:) = alpha{s};
    end
end