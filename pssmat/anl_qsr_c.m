function output = anl_qsr_c(varargin)
% function output = anl_qsr_c(varargin)
%
% Analyzes the QSR-dissipativity of a continuous-time linear system affected
% by polytopic uncertainty using linear matrix inequalities (LMIs).
%
% The LMIs are constructed using ROLMIP (which internally relies on YALMIP)
% and can be solved with any SDP solver supported by YALMIP (SeDuMi is used
% by default).
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
%   output = anl_qsr_c(A,B,C,D)
%   output = anl_qsr_c(A,B,C,D,'param',value,...)
%
%   output = anl_qsr_c()
%   output = anl_qsr_c('param',value,...)
%
%   options.param = value;
%   output = anl_qsr_c(options);
%
% If A,B,C,D are not provided, the function automatically runs a built-in
% example (2-state system with 2 vertices). Optional parameters can still
% be specified in this case.
%
% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
%   A,B,C,D   -> System matrices at the vertices of the polytope.
%                They can be provided as:
%                  - Cell arrays (multi-vertex case), or
%                  - Numeric matrices (single-vertex case)
%
% Optional parameters (passed as name-value pairs):
%
%   'mode'    -> Type of dissipativity test:
%                'passive' : standard passivity test
%                'input'   : computes input passivity index (delta)
%                'output'  : computes output passivity index (epsilon)
%                'strict'  : computes both indices
%
%   'tol'     -> Feasibility tolerance for the LMIs (default: 1e-7). This
%                tolerance is enforced on scalar decision variables (e.g.,
%                passivity indices) when the problem involves optimization.
%
%   'solver'  -> SDP solver used by YALMIP (default: 'sedumi')
%
%   'degP'    -> Degree of the parameter-dependent Lyapunov matrix
%                (default: 0)
%
% -------------------------------------------------------------------------
% OUTPUTS
% -------------------------------------------------------------------------
%   output.feas     -> Feasibility flag (1 if feasible, 0 otherwise)
%   output.P        -> Lyapunov matrix solution
%   output.del      -> Input passivity index (if applicable)
%   output.eps      -> Output passivity index (if applicable)
%   output.cpusec_s -> Solver CPU time (seconds)
%   output.cpusec_m -> LMI assembly time (seconds)
%   output.V        -> Number of scalar decision variables
%   output.L        -> Number of LMI rows
%   output.res      -> Minimum primal residual
%
% -------------------------------------------------------------------------
% EXAMPLES
% -------------------------------------------------------------------------
% The function uses a built-in example whenever the system matrices A,B,C,D
% are not provided. In this case, any additional parameters are still applied
% normally.
%
% 1) Run built-in example (default settings: passive test)
%
%   output = anl_qsr_c();
%
% 2) Run built-in example computing the input passivity index
%
%   output = anl_qsr_c('mode','input');
%
% 3) Run built-in example with custom solver and tolerance
%
%   output = anl_qsr_c('mode','strict','solver','sdpt3','tol',1e-6);
%
% 4) Custom system with 2 vertices
%
%   A = cell(2,1);
%   B = cell(2,1);
%   C = cell(2,1);
%   D = cell(2,1);
%
%   A{1} = [-1.99  0.98 ;
%           -0.28 -1.77];
%   A{2} = [-1.30  0.46 ;
%           -0.49  0.09];
%
%   B{1} = [-1.26 ;
%           -1.83];
%   B{2} = [-1.45 ;
%            1.88];
%
%   C{1} = [-0.32 0.15];
%   C{2} = [-0.16 0.01];
%
%   D{1} = 4.54;
%   D{2} = 1.77;
%
%   output = anl_qsr_c(A,B,C,D,'mode','input');
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
% - If numeric matrices are provided instead of cell arrays, they are
%   automatically converted to single-vertex polytopes.
%
% - The input and output passivity indices are obtained by maximizing their
%   respective decision variables.
%
% -------------------------------------------------------------------------
% Date: 08/04/2026
% Author: silva.artur@protonmail.com

    if nargin == 0
        A = []; B = []; C = []; D = [];
        args = {};

    elseif nargin >= 1 && (~iscell(varargin{1}) && ~isnumeric(varargin{1}))
        A = []; B = []; C = []; D = [];
        args = varargin;

    else
        if nargin < 4
            error('You must provide A,B,C,D or call anl_qsr_c() for demo.');
        end

        A = varargin{1};
        B = varargin{2};
        C = varargin{3};
        D = varargin{4};
        args = varargin(5:end);
    end

    if isempty(A)
        A = cell(2,1);
        B = cell(2,1);
        C = cell(2,1);
        D = cell(2,1);

        A{1} = [-1.99  0.98 ;
                -0.28 -1.77];
        A{2} = [-1.30  0.46 ;
                -0.49  0.09];

        B{1} = [-1.26 ;
                -1.83];
        B{2} = [-1.45 ;
                 1.88];

        C{1} = [-0.32 0.15];
        C{2} = [-0.16 0.01];

        D{1} = 4.54;
        D{2} = 1.77;
    end
    
    options = struct();

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
        options.tol    = 1e-7;
    end
    if ~isfield(options,'mode')
        options.mode   = 'passive';
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

    switch(options.mode)
        case 'passive'
            obj  = [ ]                       ;
            Q    = zeros  (outputs,outputs)  ;
            S    = 0.5*eye(outputs,inputs)   ;
            R    = zeros  (inputs,inputs)    ;
        case 'input'
            del  = sdpvar(1)                 ;
            LMIs = [LMIs, del >= options.tol];
            obj  = -del                      ;
            Q    =  zeros  (outputs,outputs) ;
            S    =  0.5*eye(outputs,inputs)  ;
            R    = -del*eye(inputs,inputs)   ;
        case 'output'
            eps  = sdpvar(1)                 ;
            LMIs = [LMIs, eps >= options.tol];
            obj  = -eps                      ;
            Q    = -eps*eye(outputs,outputs) ;
            S    =  0.5*eye(outputs,inputs)  ;
            R    =  zeros  (inputs,inputs)   ;
        case 'strict'
            eps  = sdpvar(1)                 ;
            del  = sdpvar(1)                 ;
            LMIs = [LMIs, eps >= options.tol];
            LMIs = [LMIs, del >= options.tol];
            obj  = -(eps + del)              ;
            Q    = -eps*eye(outputs,outputs) ;
            S    =  0.5*eye(outputs,inputs)  ;
            R    = -del*eye(inputs,inputs)   ;
    end

    P = rolmipvar(order,order,'P','symmetric',vertices,options.degP);
    LMIs = [LMIs, P >= 0]        ;
            
    Qc = C'*Q*C                  ;
    Sc = C'*Q*D + C'*S           ;
    Rc = D'*Q*D + D'*S + S'*D + R;

    ps11 =  A'*P + P*A - Qc;
    ps12 =  P*B  - Sc      ;
    ps22 = -Rc             ;
	    
    Ps   = [ps11  ps12 ;
            ps12' ps22];
	    
    LMIs = [LMIs, Ps <= 0];

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
        output.P = double(P);
        output.feas = 1;
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
