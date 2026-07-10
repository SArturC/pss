function output = gen_poly(order, vertices, max_iter, type)
    if nargin < 4
        type = 'continuous';
    end

    found = false;
    output = [];

    for k = 1:max_iter
        A = cell(1,vertices);
        lambda = zeros(vertices,1);
        lambda_vertex = zeros(vertices,1);
    
        for i = 1:vertices
            A{i} = randn(order);
            eigs_i = eig(A{i});
        
            switch lower(type)
                case 'continuous'
                    [~,idx] = max(real(eigs_i));
        
                    lambda(i) = real(eigs_i(idx));
                    lambda_vertex(i) = eigs_i(idx);
                case 'discrete'
                    [~,idx] = max(abs(eigs_i));
        
                    lambda(i) = abs(eigs_i(idx));
                    lambda_vertex(i) = eigs_i(idx);
                otherwise
                    error('Tipo deve ser continuous ou discrete')
            end
        end

        [lambda_max, id_vertice] = max(lambda);

        switch lower(type)
            case 'continuous'
                l_max = real(lambda_vertex(id_vertice));
                if (l_max > 0)
                    continue;
                end
            case 'discrete'
                l_max = abs(lambda_vertex(id_vertice));
                if (l_max > 1)
                    continue;
                end
        end

        for j = 1:max_iter
            alpha = simplex_random(vertices);
            Aalpha = 0;
            
            for i = 1:vertices
                Aalpha = Aalpha + alpha(i)*A{i};
            end

            switch lower(type)
                case 'continuous'
                    if ~(max(real(eig(Aalpha))) > lambda_max)
                        continue;
                    end
                case 'discrete'
                    if ~(max(abs(eig(Aalpha))) > lambda_max)
                        continue;
                    end
            end
        end

        found = true;
    
        output.A    = A;
        output.iter = k;
        output.lambda_max_vertece = lambda_vertex(id_vertice);
        break;
    end
    
    if ~found
        error('Exemplo nao encontrado.');
    end
end

% FUNÇÃO AUXILIAR

function x = simplex_random(N)
    x = zeros(1,N);
    x(1) = 1-rand^(1/(N-1));
    
    for k = 2:N-1
        x(k) = (1-sum(x(1:k-1)))*(1-rand^(1/(N-k)));
    end
    
    x(N) = 1-sum(x(1:N-1));
end
