function output = anl_hinf_c(A,B,C,D,varargin)
    options = [];
    output = [];
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
    if ~isfield(options,'mu')
        options.mu     = 0;
    end

    if ~iscell(A)
        A = {A};
    end
    if ~iscell(B)
        B = {B};
    end
    if ~iscell(C)
        C = {C};
    end
    if ~iscell(D)
        D = {D};
    end
    
    vertices = length(A)   ;
    order    = size(A{1},1);
    inputs   = size(B{1},2);
    outputs  = size(C{1},1);
	
    output.cpusec_m = clock;

    A = rolmipvar(A,'A',vertices,1);
    B = rolmipvar(B,'B',vertices,1);
    C = rolmipvar(C,'C',vertices,1);
    D = rolmipvar(D,'D',vertices,1);
    LMIs = [];

    P = rolmipvar(order,order,'P','symmetric',vertices,0);
    if ~options.mu
        mu   = sdpvar(1);
        mur  = rolmipvar(mu,'mu',vertices,0);
    else
        mu   = options.mu;
        mur  = options.mu;
    end

    LMIs = [LMIs, P >= 0];

    h11 =  A'*P + P*A + C'*C;
    h12 =  P*B  + C'*D      ;
    h22 =  D'*D - mur*eye(inputs);
	    
    H = [h11  h12 ;
         h12' h22];
	    
    LMIs = [LMIs, H <= 0];

    obj = mu;

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
        output.P = double(P);
        output.normHinf = sqrt(double(mu));
    end
end
