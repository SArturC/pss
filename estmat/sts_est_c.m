function output = sts_est_c(A,B,varargin)
% function output = sts_est_c(A,Bu,varargin)
%
% Synthesizes a state-feedback controller that renders a continuous-time
% linear system stable in the presence of polytopic uncertainty
% using linear matrix inequalities (LMIs). The LMIs are formulated using
% YALMIP and can be solved with any SDP solver supported by YALMIP
% (SeDuMi/MOSEK can be used with default options).
%
% The function computes a state-feedback gain K such that
%
%       u = Kx
%
% ensures the closed-loop system is Lyapunov stable.
%
% Additionally, an optional bound on the norm of the gain K can be imposed
% through a convex relaxation. This introduces auxiliary variables and
% bilinear terms, which are handled by fixing a scalar parameter (gnorm)
% externally.
%
% Inputs:
%   A,Bu -> Vertices of the uncertain polytope (cell arrays)
%           describing the system
%
%           x_dot = A x + Bu u
%
%   (Note: disturbance/output channels are not considered in this function.)
%
% Optional parameters:
%
%   tol         -> Feasibility tolerance for the LMIs.
%
%   solver      -> SDP solver used to solve the LMIs.
%
%   technique   -> LMI formulation used for controller synthesis:
%                  'qdr'  : quadratic Lyapunov relaxation
%                  'fins' : full-information null-space formulation
%
%   xi          -> Set of scalar parameters used in the 'fins' technique
%                  (default: logspace(-3,3,7)).
%
%   knorm       -> Upper bound on the squared norm of the controller gain K.
%                  If knorm > 0, additional LMIs are introduced to enforce
%                  a constraint of the type:
%
%                      ||K||^2 <= knorm
%
%                  via a convex relaxation.
%
%   gnorm       -> Fixed scalar parameter used to linearize the norm constraint.
%                  This parameter removes bilinearity in the LMIs by fixing
%                  part of the constraint.
%
%                  IMPORTANT:
%                  The pair (knorm, gnorm) introduces a non-convex coupling.
%                  Therefore, gnorm must be selected externally.
%
%                  The recommended approach is to define a grid of values
%                  for gnorm and call this function for each value until
%                  a feasible solution is found.
%
% Outputs:
%
%   output.cpusec_s -> CPU time required to solve the LMIs (seconds).
%   output.cpusec_m -> CPU time required to assemble the LMIs (seconds).
%   output.Pv       -> Lyapunov matrix obtained from the synthesis problem.
%   output.K        -> State-feedback gain matrix.
%   output.Knorm    -> Norm of gain K.
%   output.V        -> Number of scalar decision variables in the problem.
%   output.L        -> Number of LMI rows used in the optimization problem.
%   output.res      -> Minimum primal residual returned by the LMI solver.
%   output.feas     -> Feasibility flag (1 if feasible, 0 otherwise).
%
% Example: system with 2 states and 2 vertices
%
%   A{1}  = [-0.9 0.2; -0.5 -1.9];
%   A{2}  = [-1.1 0.1; -0.4 -2.0];
%
%   Bu{1} = [1;0];
%   Bu{2} = [1;0];
%
%   output = sts_est_c(A,Bu)
%
% Example with norm constraint:
%
%   knorm = 10;
%   Ggrid = logspace(-2,2,20);
%
%   for g = Ggrid
%       out = sts_est_c(A,Bu,...
%           'knorm',knorm,'gnorm',g);
%       if out.feas
%           break
%       end
%   end
%
% Date: 27/03/2026
% Author: a298980@dac.unicamp.br

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
    if ~isfield(options,'technique')
        options.technique = 'qdr';
    end
    if ~isfield(options,'xi')
        options.xi = 1;
    end
    if ~isfield(options,'knorm')
        options.knorm = 0;
    end
    if ~isfield(options,'gnorm')
        options.gnorm = 0;
    end
    if ~isfield(options,'degP')
        options.degP = 0;
    end
    
    if ~iscell(A)
        A = {A};
    end
    if ~iscell(B)
        B = {B};
    end

    vertices = length(A);
    order    = size(A{1},1);
    inputs_u = size(B{1},2);

    output.cpusec_m = clock;

    A = rolmipvar(A,'A',vertices,1);
    B = rolmipvar(B,'B',vertices,1);
    LMIs = [];
    obj  = [];

    switch options.technique

        case 'qdr'
            W = rolmipvar(order   ,order,'W','symmetric',vertices,0);
            Z = rolmipvar(inputs_u,order,'Z','full'     ,vertices,0);
            LMIs = [LMIs, W >= 0];

            Es = A*W + W*A' + B*Z + Z'*B';
            LMIs = [LMIs, Es <= 0];
            
            if options.knorm > 0
                kz   = sdpvar(1);
                kzr  = rolmipvar(kz,'kz',vertices,0);
                kw   = options.gnorm;
            	
                Zlim = [-kzr*eye(order) Z'; Z -eye(inputs_u)];
                Wlim = [kw*eye(order) eye(order); eye(order) W];
                
                LMIs = [LMIs, Zlim <= 0, Wlim >= 0, kz*kw*kw <= options.knorm];
                obj  = obj - kz;
            end

        case 'fins'
            W = rolmipvar(order   ,order,'W','symmetric',vertices,options.degP);
            X = rolmipvar(order   ,order,'X','full'     ,vertices,0           );
            Z = rolmipvar(inputs_u,order,'Z','full'     ,vertices,0           );
            xi = rolmipvar(options.xi,'xi',vertices,0);
            LMIs = [LMIs, W >= 0];

            Lam = A*X + B*Z  ;
                    
            es11 = Lam + Lam'    ;
            es12 = W -X' + xi*Lam;
            es22 = -xi*(X + X')  ;

            Es = [es11  es12 ;
                  es12' es22];

            LMIs = [LMIs, Es <= 0];

            if options.knorm > 0
	            beta = rolmipvar(1,1,'beta','full',vertices,0);
	            mu   = options.gnorm;
	        
                Klim = [X+X'-mu*eye(order) Z'; Z beta*eye(inputs_u)];
                
                LMIs = [LMIs, Klim >= 0, beta*(1/mu) <= options.knorm];
                obj = obj - double(beta);
            end
    end
    
    output.L = 0;
    for i=1:size(LMIs,1)
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
        W = double(W);
        Z = double(Z);
        
        output.W = W; % pesquisar como extrair os vertices para poder inverter
        output.feas = 1;
        
        switch(options.technique)
            case 'qdr'
                output.K = Z / W;
            case 'fins'
                X = double(X);
                output.K = Z / X;
                output.xi = double(xi);
        end
        
        output.Knorm = norm(output.K, 2);
    end
end
