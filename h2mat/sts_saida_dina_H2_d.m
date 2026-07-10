function output = sts_saida_dina_H2_d(A,B1,B2,C1,D11,D12,C2,D21,varargin)
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
    Z = rolmipvar(outputs_z,outputs_z,'Z','symmetric',vertices,0);
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

    t11 =  Z                  ;
    t12 =  C1*X + D12*L       ;
    t13 =  C1 + D12*R*C2      ;

    t22 =  X + X' - P         ;
    t23 =  eye(order) + S' - J;

    t33 =  Y + Y' - H         ;

    Tr = [t11  t12  t13 ;
          t12' t22  t23 ;
          t13' t23' t33];

    g11 =  P                    ;
    g12 =  J                    ;
    g13 =  A*X + B2*L           ;
    g14 =  A + B2*R*C2          ;
    g15 =  B1 + B2*R*D21        ;
    
    g22 =  H                    ;
    g23 =  Q                    ;
    g24 =  Y*A + F*C2           ;
    g25 =  Y*B1 + F*D21         ;
    
    g33 =  X + X' - P           ;
    g34 =  eye(order) + S' - J  ;
    g35 =  zeros(order,inputs_w);
    
    g44 =  Y + Y' - H           ;
    g45 =  zeros(order,inputs_w);
    
    g55 =  eye(inputs_w)        ;
        
    Gr = [g11  g12  g13  g14  g15  ;
          g12' g22  g23  g24  g25  ;
          g13' g23' g33  g34  g35  ;
          g14' g24' g34' g44  g45  ;
          g15' g25' g35' g45' g55 ];

    Dr = D11 + D12*R*D21;
        
    LMIs = [LMIs, Tr >= 0, Gr >= 0, Dr == 0, trace(Z) <= mur];
    
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
        output.normH2 = sqrt(double(mu))           ;
    end
end

