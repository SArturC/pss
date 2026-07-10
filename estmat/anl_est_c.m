function output = anl_est_c(A,varargin)
% function output = anl_est_c(A,varargin)
%
% Esta função verifica a estabilidade robusta de sistemas lineares
% contínuos sujeitos a incerteza politópica utilizando desigualdades
% matriciais lineares (LMIs) formuladas com YALMIP.
%
% O sistema é descrito por um conjunto de vértices da matriz de estado:
%
%        x_dot = A(α)x,   A(α) ∈ conv{A1, A2, ..., AN}
%
% onde conv{·} denota o fecho convexo.
%
% Diferentes condições (Teorema 1 e Lemas 1–4) podem ser selecionadas
% através do parâmetro 'technique'.
%
% -------------------------------------------------------------------------
% ENTRADA:
% -------------------------------------------------------------------------
%   A : célula contendo os vértices {A1, A2, ..., AN}, onde cada Ai ∈ R^{n×n}
%
% PARÂMETROS OPCIONAIS (passados como pares 'chave',valor):
%
%   'solver'    : solver SDP utilizado pelo YALMIP
%                 (default: 'sedumi')
%
%   'technique' : condição de estabilidade a ser testada:
%                 't1'  - Teorema 1 (Lyapunov quadrático)
%                 'l1'  - Lema 1 (relaxação com variáveis auxiliares)
%                 'l2'  - Lema 2 (relaxação com variável V)
%                 'l3'  - Lema 3 (condição estrita)
%                 'l4'  - Lema 4 (relaxação com parâmetros politópicos)
%
%   'degP'      : grau do polinômio de Lyapunov dependente de parâmetros
%                 (usado em rolmipvar)
%                 (default: 0)
%
% -------------------------------------------------------------------------
% SAÍDA:
% -------------------------------------------------------------------------
%   output.feas      : 1 se o sistema é estável, 0 caso contrário
%   output.V         : número de variáveis escalares do problema LMI
%   output.L         : número total de linhas das LMIs
%   output.res       : menor resíduo primal das restrições
%   output.cpusec_s  : tempo de solução (segundos)
%   output.cpusec_m  : tempo de montagem das LMIs (segundos)
%   output.P         : matriz de Lyapunov (se factível)
%
% -------------------------------------------------------------------------
% DESCRIÇÃO DAS CONDIÇÕES:
% -------------------------------------------------------------------------
%   T1:  Busca uma função de Lyapunov quadrática comum:
%        A'P + PA < 0
%
%   L1–L4: Relaxações LMI que reduzem o conservadorismo introduzindo
%          variáveis auxiliares e/ou dependência paramétrica.
%
% -------------------------------------------------------------------------
% OBSERVAÇÕES:
% -------------------------------------------------------------------------
% - Requer YALMIP e um solver SDP instalado (SeDuMi, MOSEK, etc.).
% - A função utiliza rolmipvar para modelagem de incertezas politópicas.
% - O critério de factibilidade é baseado no resíduo das LMIs.
%
% -------------------------------------------------------------------------
% Autor: Artur Cesar da Silva - 298980
% Data: 07/04/2026
% -------------------------------------------------------------------------

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
    if ~isfield(options,'technique')
        options.technique   = 't1';
    end
    if ~isfield(options,'degP')
        options.degP = 0;
    end

    if ~iscell(A)
        A = {A};
    end
    
    vertices = length(A)   ;
    order    = size(A{1},1);
	
    output.cpusec_m = clock;

    A = rolmipvar(A,'A',vertices,1);
    LMIs = [];

    switch options.technique
        case 't1'
            P = rolmipvar(order,order,'P','symmetric',vertices,0);
            
            LMIs = [LMIs, P >= 0];

            E = A'*P + P*A;

            LMIs = [LMIs, E <= 0];

        case 'l1'
            P  = rolmipvar(order,order,'P','symmetric',vertices,options.degP);
            X1 = rolmipvar(order,order,'X1','full',vertices,0);
            X2 = rolmipvar(order,order,'X2','full',vertices,0);
            
            LMIs = [LMIs, P >= 0];

            e11 =  X1*A + A'*X1' ;
            e12 =  P - X1 + A'*X2;
            e22 = -X2 - X2'      ;
            
            E = [e11  e12 ;
                 e12' e22];

            LMIs = [LMIs, E <= 0];
            
        case 'l2'
            P = rolmipvar(order,order,'P','symmetric',vertices,options.degP);
            V = rolmipvar(order,order,'V','full',vertices,0);

            LMIs = [LMIs, P >= 0];

            e11 = -V - V'            ;
            e12 =  V'*A + P          ;
            e13 =  V'                ;
            e22 = -P                 ;
            e23 =  zeros(order,order);
            e33 = -P                 ;
            
            E = [e11  e12  e13 ;
                 e12' e22  e23 ;
                 e13' e23' e33];

            LMIs = [LMIs, E <= 0];
            
        case 'l3'
            P = rolmipvar(order,order,'P','symmetric',vertices,options.degP);
            LMIs = [LMIs, P >= 0];

            E = A'*P + P*A;

            LMIs = [LMIs, E <= -eye(order)];
            
        case 'l4'
            P  = rolmipvar(order,order,'P','symmetric',vertices,options.degP);
            X1 = rolmipvar(order,order,'X1','full',vertices,options.degP);
            X2 = rolmipvar(order,order,'X2','full',vertices,options.degP);
            
            LMIs = [LMIs, P >= 0];

            e11 =  X1*A + A'*X1' ;
            e12 =  P - X1 + A'*X2;
            e22 = -X2 - X2'      ;
            
            E = [e11  e12 ;
                 e12' e22];

            LMIs = [LMIs, E <= 0];
    end
    
    output.L = 0;
    for i=1:length(LMIs,1)
        output.L = output.L + size(sdpvar(LMIs(i)),1);
    end
    output.V = size(getvariables(LMIs),2);

    output.cpusec_m = etime(clock,output.cpusec_m);
    sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver',options.solver));
    output.cpusec_s = sol.solvertime;
    
    output.feas = 0;
    output.res = min(checkset(LMIs));
    
    if sol.problem == 1
        return;
    end
    
    if output.res > 0
        output.P = double(P);
        output.feas = 1;
    end
end
