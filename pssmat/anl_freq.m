function output = anl_freq(A,B,C,D,varargin)
% function output = anl_freq(A,B,C,D,varargin)
%
% Analyzes frequency-domain passivity and closed-loop stability of a 
% continuous-time linear system in the presence of polytopic uncertainty.
%
% The function evaluates the system along a convex combination of the 
% polytope vertices and computes:
%
%   - The minimum real part of the frequency response G(jw)
%   - The poles of a closed-loop transformation H/(1 + H*H)
%
% This allows assessing whether the polytope is passive and whether the 
% closed-loop system remains stable for all convex combinations.
%
% Inputs:
%   A,B,C,D -> Vertices of the uncertain polytope (cell arrays)
%              describing the system
%
%          x_dot = A x + B u
%              y = C x + D u
%
%   Optional parameters:
%
%     alpha      -> Grid of convex coefficients used to sweep the polytope
%                   (default: linspace(0,1,100))
%
%     frequency  -> Frequency grid used to evaluate G(jw)
%                   (default: logspace(-2,4,2000))
%
% Outputs:
%
%   output.pss -> Passivity flag (1 if the polytope is passive, 0 otherwise)
%                 Verified numerically as:
%                 min Re{G(jw)} >= 0 for all alpha and w
%
%   output.est -> Stability flag (1 if the closed-loop system is stable,
%                 0 otherwise)
%                 Verified numerically as:
%                 Re{lambda(H/(1+H*H))} < 0 for all alpha
%
% The function also generates:
%
%   - A plot of pole trajectories as alpha varies
%   - A plot of the minimum real part of G(jw) versus alpha
%
% Example: system with 2 vertices
%
%   A{1} = [0 1; -0.64 -0.64];
%   A{2} = [0 1; -2.00 -1.00];
%
%   B = [0.5;0.5];
%   C = [1 0];
%   D =  0.1;
%
%   output = anl_freq(A,B,C,D)
%
% Date: 19/03/2026
% Author: a298980@dac.unicamp.br

    options = [];

    if nargin > 4
        if nargin == 5
            options = varargin{1};
        else
            options = struct(varargin{:});
        end
    end

    if ~isfield(options,'alpha')
        options.alpha = linspace(0,1,100);
    end
    if ~isfield(options,'frequency')
        options.frequency = logspace(-2,4,2000);
    end

    alpha = options.alpha;
    w     = options.frequency;

    A1 = A{1};
    A2 = A{2};

    poles_real = [];
    poles_imag = [];
    passivity  = zeros(length(alpha),1);

    stable_flag   = true;
    passive_flag  = true;

    for k = 1:length(alpha)

        a = alpha(k);

        Aalpha = a*A1 + (1-a)*A2;

        H = ss(Aalpha,B,C,D);

        Hcl = H/(1 + H*H);
        p   = pole(Hcl);

        poles_real = [poles_real; real(p)'];
        poles_imag = [poles_imag; imag(p)'];

        if any(real(p) >= 0)
            stable_flag = false;
        end

        Gjw = squeeze(freqresp(H,w));
        passivity(k) = min(real(Gjw));

        if passivity(k) < 0
            passive_flag = false;
        end
    end

    output.pss = passive_flag;
    output.est = stable_flag;

    figure
    hold on
    grid on

    cmap1 = jet(length(alpha));

    for i = 1:size(poles_real,2)
        for k = 1:length(alpha)
            scatter(poles_real(k,i),poles_imag(k,i),30,cmap1(k,:),'filled')
        end
    end

    xlabel('Re\{p\}')
    ylabel('Im\{p\}')
    title('Trajetória dos polos do sistema politópico')

    colorbar
    colormap(jet)

    figure
    plot(alpha,passivity,'LineWidth',2)
    grid on

    xlabel('\alpha')
    ylabel('min Re\{G(j\omega)\}')
    yline(0,'k--')

    title('Passividade ao longo do politopo')
end