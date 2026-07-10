function output = anl_h2_markov(A,E,Cz,Ez,Pr,mu,varargin)
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
    if ~isfield(options,'ga')
        options.ga     = 0;
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
    
    modos     = length(A)   ;
    order     = size(A{1},1);
    inputs    = size(E{1},2);
    outputs_z = size(Cz{1},1);
    
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
        W{i} = sdpvar(inputs,inputs,'symmetric');
        P{i} = sdpvar(order ,order ,'symmetric');
    end

    Pp = cell(modos,1);
    Mu   =  0;

    for i = 1:modos
        first = true;
        for j = 1:modos
            if first
                Pp{i} = Pr(i,j)*P{j};
                first = false;
            else
                Pp{i} = Pp{i} + Pr(i,j)*P{j};
            end
        end

        t11 =  W{i}                  ;
        t21 =  Pp{i}*E{i}            ;
        t31 =  Ez{i}                 ;
        t22 =  Pp{i}                 ;
        t32 =  zeros(outputs_z,order);
        t33 =  eye(outputs_z)        ;
    
        T = [t11  t21' t31' ;
             t21  t22  t32' ;
             t31  t32  t33 ];
    
        g11 =  P{i}                  ;
        g21 =  Pp{i}*A{i}            ;
        g31 =  Cz{i}                 ;
        g22 =  Pp{i}                 ;
        g32 =  zeros(outputs_z,order);
        g33 =  eye(outputs_z)        ;
            
        G = [g11  g21' g31' ;
             g21  g22  g32' ;
             g31  g32  g33 ];

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
        output.feas   = 1               ;
        output.normH2 = sqrt(double(ga));
    end
end
