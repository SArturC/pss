function output = sts_saida_dina_Hinf_d(A,B1,B2,C1,D11,D12,C2,D21,varargin)
    options = [];
    output  = [];
    if nargin > 1
        if nargin == 9
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
    if ~iscell(B1)
        B1 = {B1};
    end
    if ~iscell(B2)
        B2 = {B2};
    end
    if ~iscell(C1)
        C1 = {C1};
    end
    if ~iscell(D11)
        D11 = {D11};
    end
    if ~iscell(D12)
        D12 = {D12};
    end
    if ~iscell(C2)
        C2 = {C2};
    end
    if ~iscell(D21)
        D21 = {D21};
    end
    
    vertices   = length(A)   ;
    order      = size(A{1},1);
    inputs_w   = size(B1{1},2);
    inputs_u   = size(B2{1},2);
    outputs_z  = size(C1{1},1);
    outputs_y  = size(C2{1},1);
    
    output.cpusec_m = clock;
    
    A    = rolmipvar(A  ,'A'  ,vertices,1);
    B1   = rolmipvar(B1 ,'B1' ,vertices,1);
    B2   = rolmipvar(B2 ,'B2' ,vertices,1);
    C1   = rolmipvar(C1 ,'C1' ,vertices,1);
    D11  = rolmipvar(D11,'D11',vertices,1);
    D12  = rolmipvar(D12,'D12',vertices,1);
    C2   = rolmipvar(C2 ,'C2' ,vertices,1);
    D21  = rolmipvar(D21,'D21',vertices,1);
    LMIs = [];
    obj  = [];
    
    P = rolmipvar(order    ,order    ,'P','symmetric',vertices,0);
    H = rolmipvar(order    ,order    ,'H','symmetric',vertices,0);
    X = rolmipvar(order    ,order    ,'X','full'     ,vertices,0);
    Y = rolmipvar(order    ,order    ,'Y','full'     ,vertices,0);
    L = rolmipvar(inputs_u ,order    ,'L','full'     ,vertices,0);
    F = rolmipvar(order    ,outputs_y,'F','full'     ,vertices,0);
    Q = rolmipvar(order    ,order    ,'Q','full'     ,vertices,0);
    R = rolmipvar(inputs_u ,outputs_y,'R','full'     ,vertices,0);
    S = rolmipvar(order    ,order    ,'S','full'     ,vertices,0);
    J = rolmipvar(order    ,order    ,'J','full'     ,vertices,0);
    
    if ~options.mu
        mu   = sdpvar(1);
        mur  = rolmipvar(mu,'mu',vertices,0);
        obj  = mu;
    else
        mu   = options.mu;
        mur  = options.mu;
    end

    m11 =  P                     ;
    m12 =  J                     ;
    m13 =  A*X + B2*L            ;
    m14 =  A + B2*R*C2           ;
    m15 =  B1 + B2*R*D21         ;
    m16 =  zeros(order,outputs_z);
    
    m22 =  H                     ;
    m23 =  Q                     ;
    m24 =  Y*A + F*C2            ;
    m25 =  Y*B1 + F*D21          ;
    m26 =  zeros(order,outputs_z);
    
    m33 =  X + X' - P            ;
    m34 =  eye(order) + S' - J   ;
    m35 =  zeros(order,inputs_w) ;
    m36 =  X'*C1' + L'*D12'      ;
    
    m44 =  Y + Y' - H            ;
    m45 =  zeros(order,inputs_w) ;
    m46 =  C1' + C2'*R'*D12'     ;
    
    m55 =  eye(inputs_w)         ;
    m56 =  D11' + D21'*R'*D12'   ;

    m66 =  mur*eye(outputs_z)    ;
        
    M = [m11  m12  m13  m14  m15  m16 ;
         m12' m22  m23  m24  m25  m26 ;
         m13' m23' m33  m34  m35  m36 ;
         m14' m24' m34' m44  m45  m46 ;
         m15' m25' m35' m45' m55  m56 ;
         m16' m26' m36' m46' m56' m66];
        
    LMIs = [LMIs, M >= 0];
    
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
        output.feas = 1        ;
        Y = double(Y)          ;
        X = double(X)          ;
        Q = double(Q)          ;
        F = double(F)          ;
        R = double(R)          ;
        L = double(L)          ;
        U = eye(order)         ;
        V = (double(S) - Y*X)/U;

        A  = double(A) ;
        B2 = double(B2);
        C2 = double(C2);

        Ctr1 = [inv(V) -V\Y*B2; zeros(inputs_u,order) eye(inputs_u)]        ;
        Ctr2 = [Q-Y*A*X F; L R]                                             ;
        Ctr3 = [U\eye(order) zeros(order,outputs_y); -C2*X/U eye(outputs_y)];
        Ctr  = Ctr1*Ctr2*Ctr3                                               ;

        output.Ac = Ctr(1:order,1:order)           ;
        output.Bc = Ctr(1:order,...
                        order+1:order+outputs_y)   ;
        output.Cc = Ctr(order+1:order+inputs_u,...
                        1:order)                   ;
        output.Dc = Ctr(order+1:order+inputs_u,...
                        order+1:order+outputs_y)   ;
        output.normHinf = sqrt(double(mu))         ;
    end
end
