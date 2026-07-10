function output = flt_h2_markov(A,E,Cz,Ez,Cy,Ey,P,mu,varargin)
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
    if ~isfield(options,'ga')
        options.ga     = 0;
    end
    if ~isfield(options,'modeDependent')
        options.modeDependent = true;
    end
    if ~isfield(options,'strictlyProper')
        options.strictlyProper = false;
    end
    
    if ~iscell(A)
        A = {A};
    end
    if ~iscell(E)
        E = {E};
    end
    if ~iscell(Cz)
        Cz = {Cz};
    end
    if ~iscell(Ez)
        Ez = {Ez};
    end
    if ~iscell(Cy)
        Cy = {Cy};
    end
    if ~iscell(Ey)
        Ey = {Ey};
    end
    
    modos     = length(A)   ;
    order     = size(A{1},1);
    inputs    = size(E{1},2);
    outputs_z = size(Cz{1},1);
    outputs_y = size(Cy{1},1);
    
    output.cpusec_m = clock;

    LMIs = [];
    obj  = [];

    if ~options.ga
        ga  = sdpvar(1);
        obj = ga;
    else
        ga  = options.ga;
        gar = options.ga;
    end

    for i = 1:modos
        Z{i} = sdpvar(order ,order ,'symmetric');
        X{i} = sdpvar(order ,order ,'symmetric');
        W{i} = sdpvar(inputs,inputs,'symmetric');
    end
    
    if options.modeDependent
        for i = 1:modos
            if options.strictlyProper
                K{i} = zeros(outputs_z,outputs_y);
            else
                K{i} = sdpvar(outputs_z,outputs_y,'full');
            end
            H{i} = sdpvar(order    ,order    ,'symmetric');
            L{i} = sdpvar(outputs_z,order    ,'full');
            F{i} = sdpvar(order    ,outputs_y,'full');
            M{i} = sdpvar(order    ,order    ,'full');
        end
    else
        if options.strictlyProper
            K0 = zeros(outputs_z,outputs_y);
        else
            K0 = sdpvar(outputs_z,outputs_y,'full');
        end
        H0 = sdpvar(order    ,order    ,'symmetric');
        L0 = sdpvar(outputs_z,order    ,'full');
        F0 = sdpvar(order    ,outputs_y,'full');
        M0 = sdpvar(order    ,order    ,'full');
    
        for i = 1:modos
            K{i} = K0;
            H{i} = H0;
            L{i} = L0;
            F{i} = F0;
            M{i} = M0;
        end
    end

    Zp = cell(modos,1);
    Xp = cell(modos,1);
    Mu = 0;

    for i = 1:modos
        first = true;
        for j = 1:modos
            if first
                Zp{i} = P(i,j)*Z{j};
                Xp{i} = P(i,j)*X{j};
                first = false;
            else
                Zp{i} = Zp{i} + P(i,j)*Z{j};
                Xp{i} = Xp{i} + P(i,j)*X{j};
            end
        end

        t11 =  W{i}                        ;
        t21 =  Zp{i}*E{i}                  ;
        t31 =  H{i}*E{i} + F{i}*Ey{i}      ;
        t41 =  Ez{i} - K{i}*Ey{i}          ;
    
        t22 =  Zp{i}                       ;
        t32 =  zeros(order    ,order)      ;
        t42 =  zeros(outputs_z,order)      ;
    
        t33 =  H{i} + H{i}' + Zp{i} - Xp{i};
        t43 =  zeros(outputs_z,order)      ;

        t44 =  eye(outputs_z)              ;
    
        T = [t11  t21' t31' t41' ;
             t21  t22  t32' t42' ;
             t31  t32  t33  t43' ;
             t41  t42  t43  t44 ];
    
        g11 =  Z{i}                         ;
        g21 =  Z{i}                         ;
        g31 =  Zp{i}*A{i}                   ;
        g41 =  H{i}*A{i} + F{i}*Cy{i} + M{i};
        g51 =  Cz{i} - K{i}*Cy{i} + L{i}    ;
        
        g22 =  X{i}                         ;
        g32 =  Zp{i}*A{i}                   ;
        g42 =  H{i}*A{i} + F{i}*Cy{i}       ;
        g52 =  Cz{i} - K{i}*Cy{i}           ;
        
        g33 =  Zp{i}                        ;
        g43 =  zeros(order,order)           ;
        g53 =  zeros(outputs_z,order)       ;
        
        g44 =  H{i} + H{i}' + Zp{i} - Xp{i} ;
        g54 =  zeros(outputs_z,order)       ;
        
        g55 =  eye(outputs_z)               ;
            
        G = [g11  g21' g31' g41' g51'  ;
             g21  g22  g32' g42' g52'  ;
             g31  g32  g33  g43' g53'  ;
             g41  g42  g43  g44  g54'  ;
             g51  g52  g53  g54  g55  ];

        Mu = Mu + mu(i)*trace(W{i});
            
        LMIs = [LMIs, T >= 0, G >= 0];
    end

    LMIs = [LMIs, ga >= Mu];
    
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
        output.feas = 1   ;
        Af = cell(modos,1);
        Bf = cell(modos,1);
        Cf = cell(modos,1);
        Df = cell(modos,1);

        for i = 1:modos
            Af{i} = -double(H{i})\double(M{i});
            Bf{i} = -double(H{i})\double(F{i});
            Cf{i} = -double(L{i});
            Df{i} =  double(K{i});
        end

        output.Af = Af;
        output.Bf = Bf;
        output.Cf = Cf;
        output.Df = Df;
        output.normH2 = sqrt(double(ga));
    end
end
