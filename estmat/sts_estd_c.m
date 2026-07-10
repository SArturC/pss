function output = sts_estd_c(A,B,varargin)
    options = [];

    if nargin > 2
        if nargin == 3
            options = varargin{1};
        else
            options = struct(varargin{:});
        end
    end

    if ~isfield(options,'solver')
        options.solver = 'sedumi';
    end
    if ~isfield(options,'tol')
        options.tol = 1e-7;
    end
    if ~isfield(options,'mask')
        options.mask = [];
    end
    if ~isfield(options,'lemma')
        options.lemma = 1;
    end
    if ~isfield(options,'xi')
        options.xi = 1;
    end

    if ~iscell(A)
        A = {A};
    end
    if ~iscell(B)
        B = {B};
    end

    vertices = length(A);
    order    = size(A{1},1);
    inputs   = size(B{1},2);

    output.cpusec_m = clock;

    A = rolmipvar(A,'A',vertices,1);
    B = rolmipvar(B,'B',vertices,1);
    LMIs = [];
    obj  = [];

    if options.lemma == 1
        if ~isempty(options.mask)
            maskz = options.mask;
            maskw = comp_mask(maskz);
            if isequal(maskw, maskw')
                W = masked_sdpvar_sym(order,maskw);
            else
                W = sdpvar(order,order,'diagonal');
            end
            Z = masked_sdpvar_full(inputs,order,maskz);
            W = rolmipvar(W,'W',vertices,0);
            Z = rolmipvar(Z,'Z',vertices,0);
        else
            W = rolmipvar(order,order ,'W',vertices,0);
            Z = rolmipvar(inputs,order,'Z','full',vertices,0);
        end

        LMIs = [LMIs, W >= 0];
    
        Es = A*W + W*A' + B*Z + Z'*B';
        LMIs = [LMIs, Es <= 0];

    elseif options.lemma == 5
        xi = options.xi;
        W = rolmipvar(order,order,'W','symmetric',vertices,0);
        if ~isempty(options.mask)
            maskz = options.mask;
            maskx = comp_mask(maskz);
            if isequal(maskx, maskx')
                X = masked_sdpvar_full(order,order,maskx);
            else
                X = sdpvar(order,order,'diagonal');
            end
            Z = masked_sdpvar_full(inputs,order,maskz);
            X = rolmipvar(X,'X',vertices,0);
            Z = rolmipvar(Z,'Z',vertices,0);
        else
            X = rolmipvar(order,order ,'X','symmetric',vertices,0);
            Z = rolmipvar(inputs,order,'Z','full',vertices,0);
        end

        LMIs = [LMIs, W >= 0];
    
        e11 = A*X + X'*A' + B*Z + Z'*B';
        e12 = W - X' + xi*A*X + xi*B*Z ;
        e22 = -xi*(X + X')             ;

        Es = [e11  e12 ;
              e12' e22];

        LMIs = [LMIs, Es <= 0];
    end

    output.L = 0;
    for i=1:length(LMIs,1)
        output.L = output.L + size(sdpvar(LMIs(i)),1);
    end
    output.V = size(getvariables(LMIs),2);
    
    output.cpusec_m = etime(clock,output.cpusec_m);
    sol = solvesdp(LMIs,obj,sdpsettings('verbose',0,'solver',options.solver));
    output.cpusec_s = sol.solvertime;
    
    output.feas = 0;
    output.res = min(checkset(LMIs));
    
    if sol.problem == 1
        return;
    end
    
    if isempty(obj)
        options.tol = 0;
    end

    if output.res > -options.tol
        if options.lemma == 1
            output.P    = inv(double(W));
            output.feas = 1;
            output.K    = double(Z) / double(W);
        elseif options.lemma == 5
            output.P    = inv(double(W));
            output.feas = 1;
            output.K    = double(Z) / double(X);
        end
    end
end

function M = masked_sdpvar_sym(n, mask)
    M = sdpvar(1,1);
    M = repmat(M, n, n);
    for i = 1:n
        for j = i:n
            if mask(i,j)
                v = sdpvar(1, 1, 'full');
                M(i,j) = v;
                if i ~= j
                    M(j,i) = v;
                end
            else
                M(i,j) = 0;
                M(j,i) = 0;
            end
        end
    end
end

function M = masked_sdpvar_full(n,m, mask)
    M = sdpvar(1,1);
    M = repmat(M, n, m);
    for i = 1:n
        for j = 1:m
            if mask(i,j)
                v = sdpvar(1, 1, 'full');
                M(i,j) = v;
            else
                M(i,j) = 0;
            end
        end
    end
end

function mask2 = comp_mask(mask1)
    n = size(mask1,2);
    mask2 = zeros(n);
    for j = 1:n
        for k = 1:n
            if all(mask1(:,k) <= mask1(:,j))
                mask2(k,j) = 1;
            end
        end
    end
end