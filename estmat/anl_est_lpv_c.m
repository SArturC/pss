function output = anl_est_lpv_c(A,varargin)
    options = [];
    output = [];
    if nargin > 1
        if nargin == 2
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
    if ~isfield(options,'bounds')
        options.bounds     = 0;
    end

    if ~iscell(A)
        A = {A};
    end
    
    vertices = length(A)   ;
    order    = size(A{1},1);
	
    output.cpusec_m = clock;

    A = rolmipvar(A,'A',vertices,1);
    obj  = [];
    LMIs = [];

    F = rolmipvar(order,order,'F','full',vertices,0);
    G = rolmipvar(order,order,'G','full',vertices,1);
    P = rolmipvar(order,order,'P','symmetric',vertices,1);
    Pdiff = diff(P,'Pdiff',options.bounds);

    LMIs = [LMIs, P >= 0];

    e11 =  Pdiff + F*A + A'*F';
    e12 =  P - F + A'*G';
    e22 = -G - G';
	    
    E = [e11  e12 ;
         e12' e22];
	    
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
        output.P = double(P);
    end
end

