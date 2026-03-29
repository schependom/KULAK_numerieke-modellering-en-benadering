function plotres(x,r,w)
%PLOTRES plots residuals in combination with weights
% Inputs
% x     Vector holding the points in which to plot residuals and weights
% r     Residuals that are to be plotted
% w     Vector with weights that are to be plotted

    % Check input
    %
    % Als lengte van residuvector of gewichtenvector
    % niet gelijk is aan (~=) de lengte van de inputvector,
    % gooi dan een error.
    %
    if any([length(r),length(w)] ~= length(x))
        error('plotres:input', ...
                'r and w should have the same number of elements as x');  
    end

    % Meerdere grafieken op dezelfde figuur
    hold on     

    % Plot residuen
    %   
    %   r   <-> rood
    %   +   <-> plusjes als punten
    %   -   <-> volle lijn
    %   --  <-> stippellijn
    plot(x,r,'r+-')
    
    % x is gewoon een linspace met daarin het discrete
    % interval waarin we willen plotten:
    %
    %   e.g. x = [-1, -0.9, -0.8, ..., 1]
    %          = [-1  -0.8  -0.8  ...  1]
    % 
    % (Vector kan zowel zonder als met komma's).
 
    % Plot de x-as:
    % Trek een lijn tussen de punten (x_1, 0) en (x_2, 0)
    %
    %   x-coordinaten: x_1, x_2
    %   y-coordinaten: 0,   0
    %
    % De k staat voor zwart ('blac-K')
    plot([x(1) x(end)], [0 0], 'k')

    % Alpha = afstand van x-as (0) tot maximale residu (max(abs(r)))
    scalefactor = max(abs(r));
    
    % Plot de herschaalde gewichten
    %
    %   g   <-> green
    %   --  <-> stippelllijn
    plot(x,w/max(w)*scalefactor,'g--')

    % Symmetrische ylim zorgt voor verticale centrering van de x-as
    ylim([-scalefactor scalefactor])
    
    hold off
end
