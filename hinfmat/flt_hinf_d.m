function output = flt_hinf_d(A,B,C1,D11,C2,D21,varargin)
    options = [];
    output  = [];
    if nargin > 1
        if nargin == 7
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
    if ~iscell(C1)
        C1 = {C1};
    end
    if ~iscell(D11)
        D11 = {D11};
    end
    if ~iscell(C2)
        C2 = {C2};
    end
    if ~iscell(D21)
        D21 = {D21};
    end
    
    vertices   = length(A)   ;
    order      = size(A{1},1);
    inputs     = size(B{1},2);
    outputs_z  = size(C1{1},1);
    outputs_y  = size(C2{1},1);
	
    output.cpusec_m = clock;

    A    = rolmipvar(A  ,'A'  ,vertices,1);
    B    = rolmipvar(B  ,'B'  ,vertices,1);
    C1   = rolmipvar(C1 ,'C1' ,vertices,1);
    D11  = rolmipvar(D11,'D11',vertices,1);
    C2   = rolmipvar(C2 ,'C2' ,vertices,1);
    D21  = rolmipvar(D21,'D21',vertices,1);
    LMIs = [];
    obj  = [];

    Z  = rolmipvar(order    ,order    ,'Z' ,'symmetric',vertices,0);
    X  = rolmipvar(order    ,order    ,'X' ,'symmetric',vertices,0);
    F  = rolmipvar(outputs_z,order    ,'F' ,'full'     ,vertices,0);
    L  = rolmipvar(order    ,outputs_y,'L' ,'full'     ,vertices,0);
    G  = rolmipvar(order    ,order    ,'G' ,'full'     ,vertices,0);
    Df = rolmipvar(outputs_z,outputs_y,'Df','full'     ,vertices,0);
    
    if ~options.mu
        mu   = sdpvar(1);
        mur  = rolmipvar(mu,'mu',vertices,0);
        obj  = mu;
    else
        mu   = options.mu;
        mur  = options.mu;
    end

    LMIs = [LMIs, Z >= 0, X >= 0];

    h11 =  Z                 ;
    h12 =  Z                 ;
    h13 =  Z*A               ;
    h14 =  Z*A               ;
    h15 =  Z*B               ;
    h16 =  zeros(order,inputs);
    
    h22 =  X                 ;
    h23 =  X*A + L*C2 + G    ;
    h24 =  X*A + L*C2        ;
    h25 =  X*B + L*D21       ;
    h26 =  zeros(order,inputs);
    
    h33 =  Z                 ;
    h34 =  Z                 ;
    h35 =  zeros(order,inputs);
    h36 =  C1' - C2'*Df' - F';
    
    h44 =  X                 ;
    h45 =  zeros(order,inputs);
    h46 =  C1' - C2'*Df'     ;
    
    h55 =  eye(inputs)       ;
    h56 =  D11' - D21'*Df'   ;
    
    h66 =  mur*eye(outputs_z);
	    
    H = [h11  h12  h13  h14  h15  h16 ;
         h12' h22  h23  h24  h25  h26 ;
         h13' h23' h33  h34  h35  h36 ;
         h14' h24' h34' h44  h45  h46 ;
         h15' h25' h35' h45' h55  h56 ;
         h16' h26' h36' h46' h56' h66 ;];
	    
    LMIs = [LMIs, H >= 0];

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
        output.feas = 1                          ;
        U = eye(order)                           ;
        V = U'\(eye(order) - double(X)/double(Z));
        output.Af = U'\(double(G)/(V*double(Z))) ;
        output.Bf = U'\double(L)                 ;
        output.Cf = double(F)/(V*double(Z))      ;
        output.Df = double(Df)                   ;
        output.normHinf = sqrt(double(mu))       ;
    end
end
