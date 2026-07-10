function [output, history] = anl_est_c_setor(A, it_max)

V0 = eye(length(A));
fila = {struct('V', V0, 'Abar', {A}, 'it', 0)};

history.V = {};
history.it = [];
history.unstable_point = [];

found_stable = false;
hit_limit = false;
best_it = inf;

while ~isempty(fila)

    node = fila{1};
    fila(1) = [];

    % 🔹 salvar simplex visitado
    history.V{end+1} = node.V;
    history.it(end+1) = node.it;

    % 🔴 teste de vértices
    for i = 1:length(node.Abar)
%         out_v = anl_est_c({node.Abar{i}});
        lambda = eig(node.Abar{i});
        if max(real(lambda)) >= 0

            history.unstable_point = node.V(:,i);

            output.feas = -1;
            output.it   = node.it;
            return;
        end
    end

    % 🟢 teste politópico
    out = anl_est_c(node.Abar,'solver','mosek');

    if out.feas == 1
        found_stable = true;
        best_it = min(best_it, node.it);
        continue;
    end

    if node.it >= it_max
        hit_limit = true;
        continue;
    end

    % 🔁 subdivisão
    subV = subdivide_triangle(node.V);

    for k = 1:4
        Abar_new = update_vertices(subV{k}, node.Abar);

        fila{end+1} = struct( ...
            'V', subV{k}, ...
            'Abar', {Abar_new}, ...
            'it', node.it+1);
    end
end

if hit_limit
    output.feas = 0; % indeterminado
    output.it   = it_max;
else
    output.feas = 1; % realmente estável
    output.it   = best_it;
end

end

function Abar_new = update_vertices(V_new, Abar_old)
    N = length(Abar_old);
    Abar_new = cell(1,N);

    for i = 1:N
        Abar_new{i} = zeros(size(Abar_old{1}));
        for j = 1:N
            Abar_new{i} = Abar_new{i} + V_new(j,i) * Abar_old{j};
        end
    end
end

function subV = subdivide_triangle(V)
    V1 = V(:,1);
    V2 = V(:,2);
    V3 = V(:,3);

    M12 = (V1 + V2)/2;
    M23 = (V2 + V3)/2;
    M31 = (V3 + V1)/2;

    S1 = [V1,  M12, M31];
    S2 = [M12, V2,  M23];
    S3 = [M31, M23, V3];
    S4 = [M12, M23, M31];

    subV = {S1, S2, S3, S4};
end

