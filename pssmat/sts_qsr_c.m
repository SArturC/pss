function output = sts_qsr_c(varargin)
% function output = sts_qsr_c(varargin)
%
% Synthesizes a state-feedback controller that renders a continuous-time
% linear system QSR-dissipative in the presence of polytopic uncertainty
% using linear matrix inequalities (LMIs).
%
% The LMIs are constructed using ROLMIP (which internally relies on YALMIP)
% and can be solved with any SDP solver supported by YALMIP (SeDuMi is used
% by default).
%
% The function computes a state-feedback gain K such that:
%
%       u = Kx
%
% ensures the closed-loop system satisfies a specified QSR-dissipativity
% condition with respect to disturbance input w and output y.
%
% -------------------------------------------------------------------------
% REQUIREMENTS
% -------------------------------------------------------------------------
% - ROLMIP
% - YALMIP
% - SDP solver
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   output = sts_qsr_c(A,Bu,Bw,C,Du,Dw)
%   output = sts_qsr_c(A,Bu,Bw,C,Du,Dw,'param',value,...)
%
%   output = sts_qsr_c()
%   output = sts_qsr_c('param',value,...)
%
%   options.param = value;
%   output = sts_qsr_c(options);
%
% If system matrices are not provided, a built-in example is used.
% Optional parameters can still be specified in this case.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   A,Bu,Bw,C,Du,Dw -> System matrices at the vertices of the polytope:
%
%       x_dot = A x + Bu u + Bw w
%           y = C x + Du u + Dw w
%
%   They can be provided as:
%     - Cell arrays (multi-vertex case)
%     - Numeric matrices (single-vertex case)
%
% Optional parameters (name-value pairs or struct):
%
%   'mode'        -> Dissipativity condition:
%                    'passive' : standard passivity test (default)
%                    'input'   : maximize input passivity index (delta)
%                    'output'  : maximize output passivity index (epsilon)
%                    'strict'  : maximize both indices
%
%   'tol'         -> Feasibility tolerance (default: 1e-7). Used as lower
%                    bound on scalar decision variables in optimization.
%
%   'solver'      -> SDP solver (default: 'sedumi')
%
%   'technique'   -> Synthesis method:
%                    'qdr'  : quadratic relaxation (default)
%                    'fins' : full-information null-space
%
%   'xi'          -> Scalar parameter used in the 'fins' formulation
%                    (default: 1)
%
%                    IMPORTANT:
%                    This parameter introduces non-convexity when treated
%                    as a decision variable. Therefore, xi must be fixed
%                    externally.
%
%                    A practical approach is to define a grid of values for
%                    xi and call the function repeatedly until a feasible
%                    solution is found.
%
%   'knorm'       -> Upper bound on ||K||^2 (default: 0 = inactive)
%
%   'gnorm'       -> Fixed scalar used to convexify the norm constraint
%                    (default: 0)
%
%                    This parameter removes bilinearity in the LMIs by
%                    fixing part of the constraint.
%
%   'degP'        -> Degree of Lyapunov matrix (default: 0)
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%   output.feas     -> Feasibility flag
%   output.W        -> Inverse of the Lyapunov matrix
%   output.K        -> State-feedback gain
%   output.Knorm    -> Norm of K
%   output.del      -> Input passivity index (if applicable)
%   output.eps      -> Output passivity index (if applicable)
%   output.cpusec_s -> Solver time
%   output.cpusec_m -> Assembly time
%   output.V        -> Number of decision variables
%   output.L        -> Number of LMI rows
%   output.res      -> Minimum residual
%
% -------------------------------------------------------------------------
% EXAMPLES
% -------------------------------------------------------------------------
% Built-in example is used when system matrices are omitted.
%
% 1) Default synthesis
%
%   output = sts_qsr_c();
%
% 2) Input passivity synthesis
%
%   output = sts_qsr_c('mode','input');
%
% 3) Using options structure
%
%   options.mode = 'strict';
%   options.solver = 'sdpt3';
%   output = sts_qsr_c(options);
%
% 4) Grid search over the xi parameter
%
%   Xigrid = logspace(-3,3,7);
%
%   for xi = Xigrid
%       out = sts_qsr_c('technique','fins','xi',xi);
%       if out.feas
%           break
%       end
%   end
%
% 5) Norm-constrained synthesis (grid search in gnorm)
%
%   knorm = 10;
%   Ggrid = logspace(-2,2,20);
%
%   for g = Ggrid
%       out = sts_qsr_c('knorm',knorm,'gnorm',g);
%       if out.feas
%           break
%       end
%   end
%
% 6) Custom system
%
%   A{1}  = [-1.95  0.71 ;
%             0.02 -4.18];
%   A{2}  = [-0.97 -2.27 ;
%            -0.39 -0.74];
% 
%   Bu{1} = [-0.23 ;
%            -2.02];
%   Bu{2} = [-2.29 ;
%            -0.69];
% 
%   Bw{1} = [-1.26 ;
%            -1.83];
%   Bw{2} = [-1.45 ;
%             1.88];
% 
%   C{1}  = [-0.21 -0.74];
%   C{2}  = [-0.09 -0.50];
% 
%   Du{1} = -0.75;
%   Du{2} = -0.50;
% 
%   Dw{1} = 4.54;
%   Dw{2} = 0.77;
%
%   output = sts_qsr_c(A,Bu,Bw,C,Du,Dw,'mode','input');
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% - Numeric matrices are automatically converted to single-vertex polytopes.
%
% - The parameters xi (in 'fins') and gnorm (in norm-constrained synthesis)
%   must be fixed externally, as they introduce non-convexity in the LMIs.
%   A grid search over these scalar parameters is recommended to recover
%   feasible solutions.
%
% - Passivity indices are obtained by maximizing their corresponding
%   scalar decision variables.
%
% -------------------------------------------------------------------------
% Date: 27/03/2026
% Author: a298980@dac.unicamp.br

    if nargin == 0
        A = []; Bu = []; Bw = []; C = []; Du = []; Dw = [];
        args = {};

    elseif nargin >= 1 && (~iscell(varargin{1}) && ~isnumeric(varargin{1}))
        A = []; Bu = []; Bw = []; C = []; Du = []; Dw = [];
        args = varargin;

    else
        if nargin < 6
            error('Provide A,Bu,Bw,C,Du,Dw or call sts_qsr_c() for demo.');
        end

        A  = varargin{1};
        Bu = varargin{2};
        Bw = varargin{3};
        C  = varargin{4};
        Du = varargin{5};
        Dw = varargin{6};

        args = varargin(7:end);
    end

    if isempty(A)
        A{1}  = [-1.95  0.71 ;
                  0.02 -4.18];
        A{2}  = [-0.97 -2.27 ;
                 -0.39 -0.74];

        Bu{1} = [-0.23 ;
                 -2.02];
        Bu{2} = [-2.29 ;
                 -0.69];

        Bw{1} = [-1.26 ;
                 -1.83];
        Bw{2} = [-1.45 ;
                  1.88];

        C{1}  = [-0.21 -0.74];
        C{2}  = [-0.09 -0.50];

        Du{1} = -0.75;
        Du{2} = -0.50;

        Dw{1} = 4.54;
        Dw{2} = 0.77;
    end

    options = [];

    if ~isempty(args)
        if length(args) == 1 && isstruct(args{1})
            options = args{1};
        else
            options = struct(args{:});
        end
    end

    if ~isfield(options,'solver')
        options.solver = 'sedumi';
    end
    if ~isfield(options,'tol')
        options.tol = 1e-7;
    end
    if ~isfield(options,'mode')
        options.mode = 'passive';
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
    if ~iscell(Bu)
        Bu = {Bu};
    end
    if ~iscell(Bw)
        Bw = {Bw};
    end
    if ~iscell(C)
        C = {C};
    end
    if ~iscell(Du)
        Du = {Du};
    end
    if ~iscell(Dw)
        Dw = {Dw};
    end

    vertices = length(A);
    order    = size(A{1},1);
    inputs_u = size(Bu{1},2);
    inputs_w = size(Bw{1},2);
    outputs  = size(C{1},1);

    output.cpusec_m = clock;
    
    A  = rolmipvar(A ,'A' ,vertices,1);
    Bu = rolmipvar(Bu,'Bu',vertices,1);
    Bw = rolmipvar(Bw,'Bw',vertices,1);
    C  = rolmipvar(C ,'C' ,vertices,1);
    Du = rolmipvar(Du,'Du',vertices,1);
    Dw = rolmipvar(Dw,'Dw',vertices,1);
    LMIs = [];

    switch(options.mode)
        case 'passive'
            obj         = [ ]                ;
            options.tol =  0                 ;
            Q    = zeros  (outputs,outputs)  ;
            S    = 0.5*eye(outputs,inputs_w) ;
            R    = zeros  (inputs_w,inputs_w);
        case 'input'
            del  = sdpvar(1)                  ;
            LMIs = [LMIs, del >= options.tol] ;
            obj  = -del                       ;
            Q    =  zeros  (outputs,outputs)  ;
            S    =  0.5*eye(outputs,inputs_w) ;
            R    = -del*eye(inputs_w,inputs_w);
        case 'output'
            eps  = sdpvar(1)                  ;
            LMIs = [LMIs, eps >= options.tol] ;
            obj  = -eps                       ;
            Q    = -eps*eye(outputs,outputs)  ;
            S    =  0.5*eye(outputs,inputs_w) ;
            R    =  zeros  (inputs_w,inputs_w);
        case 'strict'
            eps  = sdpvar(1)                  ;
            del  = sdpvar(1)                  ;
            LMIs = [LMIs, eps >= options.tol] ;
            LMIs = [LMIs, del >= options.tol] ;
            obj  = -(eps + del)               ;
            Q    = -eps*eye(outputs,outputs)  ;
            S    =  0.5*eye(outputs,inputs_w) ;
            R    = -del*eye(inputs_w,inputs_w);
    end

    switch options.technique
        case 'qdr'
            W = rolmipvar(order   ,order,'W','symmetric',vertices,0);
            Z = rolmipvar(inputs_u,order,'Z','full'     ,vertices,0);
            LMIs = [LMIs, W >= 0];

            Lam = A*W + Bu*Z;
            Sig = C*W + Du*Z;

            ps11 =  Lam' + Lam         ;
            ps12 =  Bw - Sig'*S        ;
            ps13 =  Sig'*Q             ;
            ps22 = -(Dw'*S + S'*Dw + R);
            ps23 =  Dw'*Q              ;
            ps33 =  Q                  ;
                
            if strcmp(options.mode,'passive') || strcmp(options.mode,'input')
                Ps = [ps11  ps12 ;
                      ps12' ps22];
            else
                Ps = [ps11  ps12  ps13  ;
                      ps12' ps22  ps23  ;
                      ps13' ps23' ps33 ];
            end
                
            LMIs = [LMIs, Ps <= 0];
            
            if options.knorm > 0
                kz   = sdpvar(1)                    ;
                kzr  = rolmipvar(kz,'kz',vertices,0);
                kw   = options.gnorm                ;
            	
                Zlim = [-kzr*eye(order) Z'; Z -eye(inputs_u)]  ;
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

            Lam = A*X + Bu*Z;
            Sig = C*X + Du*Z;

            ps11 = xi*(Lam' + Lam);
            ps12 = W + Lam - xi*X';
            ps13 = Bw - xi*Sig'*S ;
            ps14 = xi*Sig'*Q      ;

            ps22 = -(X + X');
            ps23 = -Sig'*S  ;
            ps24 =  Sig'*Q  ;

            ps33 = -(Dw'*S + S'*Dw + R);
            ps34 =   Dw'*Q             ;

            ps44 = Q;

            if strcmp(options.mode,'passive') || strcmp(options.mode,'input')
                Ps = [ps11  ps12  ps13  ;
                      ps12' ps22  ps23  ;
                      ps13' ps23' ps33 ];
            else
                Ps = [ps11  ps12  ps13  ps14  ;
                      ps12' ps22  ps23  ps24  ;
                      ps13' ps23' ps33  ps34  ;
                      ps14' ps24' ps34' ps44 ];
            end

            LMIs = [LMIs, Ps <= 0];
            
            if options.knorm > 0
                beta  = sdpvar(1)                        ;
	            betar = rolmipvar(beta,'beta',vertices,0);
	            mu    = options.gnorm                    ;
	        
                Klim = [X+X'-mu*eye(order) Z'; Z betar*eye(inputs_u)];
                
                LMIs = [LMIs, Klim >= 0, beta*(1/mu) <= options.knorm];
                obj = obj - beta;
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
        
        output.W = W;
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
        
        switch(options.mode)
            case 'passive'
            case 'input'
                output.del = double(del);
            case 'output'
                output.eps = double(eps);
            case 'strict'
                output.del = double(del);
                output.eps = double(eps);
        end
    end
end
