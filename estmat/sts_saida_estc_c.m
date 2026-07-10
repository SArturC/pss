function output = sts_saida_estc_c(A,B,C,varargin)
    options = [];
    output  = [];
    if nargin > 1
        if nargin == 5
            options = varargin{1};
        else
            options = struct(varargin{:});
        end
    end

    if ~isfield(options,'solver')
        options.solver = 'sedumi';
    end
    if ~isfield(options,'tol')
        options.tol    = 1e-7;
    end
    if ~isfield(options,'lemma')
        options.lemma = 2;
    end
    if options.lemma == 2
        if iscell(C)
            error('This method does not support an uncertain matrix C.')
        end
    end
    if options.lemma == 3
        if ~isfield(options,'K')
            error('This method requires a state-feedback gain K.')
        end
        K = options.K;
    end
    if ~isfield(options,'mask')
        options.mask = [];
    end
    if ~iscell(A)
        A = {A};
    end
    if ~iscell(B)
        B = {B};
    end
    
    vertices = length(A)   ;
    order    = size(A{1},1);
    inputs   = size(B{1},2);
    outputs  = size(C,1);
	
    output.cpusec_m = clock;

    LMIs = [];
    obj  = [];

    if options.lemma == 2
        R1 = zeros(order-outputs,outputs);
        R2 = eye(order-outputs);
        R  = [R1 R2];
        T  = [C; R];
    
        Ab = cell(vertices,1);
        Bb = cell(vertices,1);
        for i = 1:vertices
            Ab{i} = T*A{i}/T;
            Bb{i} = T*B{i};
        end
        Ab = rolmipvar(Ab,'Ab',vertices,1);
        Bb = rolmipvar(Bb,'Bb',vertices,1);
    
        W1 = rolmipvar(outputs      ,outputs      ,'W1' ,'symmetric',vertices,0);
        W2 = rolmipvar(order-outputs,order-outputs,'W2' ,'symmetric',vertices,0);
        Z1 = rolmipvar(inputs       ,outputs      ,'Z1' ,'full'     ,vertices,0);
    
        W  = [W1 zeros(outputs,order-outputs); zeros(outputs,order-outputs)' W2];
        Z  = [Z1 zeros(inputs ,order-outputs)];
    
        LMIs = [LMIs, W >= 0];
    
        E = Ab*W + (Ab*W)' + Bb*Z + (Bb*Z)';
    elseif options.lemma == 3
        A = rolmipvar(A,'A' ,vertices,1);
        B = rolmipvar(B,'B' ,vertices,1);
        K = rolmipvar(K,'K' ,vertices,1);
        P = rolmipvar(order ,order  ,'P','symmetric',vertices,0);
        F = rolmipvar(order ,order  ,'F','full'     ,vertices,0);
        G = rolmipvar(order ,order  ,'G','full'     ,vertices,0);
        if ~isempty(options.mask)
            maskj = options.mask;
            maskh = comp_mask(maskj);
            if isequal(maskh, maskh')
                H = masked_sdpvar_full(inputs,inputs,maskh);
            else
                H = sdpvar(inputs,inputs,'diagonal');
            end
            J = masked_sdpvar_full(inputs,outputs,maskj);
            H = rolmipvar(H,'H',vertices,0);
            J = rolmipvar(J,'J',vertices,0);
        else
            H = rolmipvar(inputs,inputs ,'H',vertices,0);
            J = rolmipvar(inputs,outputs,'J','full',vertices,0);
        end

        LMIs = [LMIs, P >= 0];

        e11 = A'*F' + F*A + K'*B'*F' + F*B*K;
        e12 = P - F + A'*G' + K'*B'*G';
        e13 = F*B + C'*J' - K'*H';

        e22 = -G -G';
        e23 = G*B;

        e33 = -H - H';

        E = [e11  e12  e13 ;
             e12' e22  e23 ;
             e13' e23' e33];
    end
	    
    LMIs = [LMIs, E <= 0];

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
        output.feas = 1;
        if options.lemma == 2
            output.Lc = double(Z1)/double(W1);
        elseif options.lemma == 3
            output.Lc = double(H)\double(J);
        end
    end
end

function maskH = comp_mask(maskJ)
    m = size(maskJ,1);
    maskH = zeros(m);
    for i = 1:m
        for j = 1:m
            if all(maskJ(j,:) <= maskJ(i,:))
                maskH(i,j) = 1;
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