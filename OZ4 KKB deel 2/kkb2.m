function c = kkb2(x, f, w, M)
%KKB2 
%
% Discrete weighted least-squares approximation of degree n for the
% function that takes the values f in the points x for the given weigths w.
%
% This method uses an efficient implementation of the normal equations to 
% solve the system of equations for LARGE NUMBER OF POINTS.
%
% Inputs
%   x     Vector holding the points in which function values are given
%   f     Function values in the given points
%   w     Vector with weights in given points
%   M     Degree of the approximating polynomial
% 
% Outputs
%   c     Coefficients of the approximating degree n polynomial
   
    % Make sure that all vectors are column vectors
    x = x(:); f = f(:); w = w(:);

    % We zoeken bij een gewogen kleinste kwadraten benadering
    % de vector c zodat 
    % 
    %   Gc = F
    %
    % waarbij
    %
    %   G_ij    = <e_i, e_j>_w 
    %           = gewogen inproduct van basisVTen e_i en e_j
    %           = (M+1) x (M+1)
    %
    %   c       = vector met M+1 coefficienten van de
    %             te benaderen veelterm y_M = sum (c_m e_m)
    %
    %   F       = vector met M+1 inproducten tussen
    %             te benaderen functiewaarden en de basisVTen
    %
    %                   F_i = <f, e_i>_w    (0 <= i <= M)
    %                   
    %   f       = [f_1, ..., f_N]
    %           = vector met N functiewaarden die we willen benaderen
    %
    %
    % We benaderen de vector f met functiewaarden 
    %
    %       -> we verkrijgen y_LS
    %
    %       -> y_LS = evaluatie van onze benaderende VT y_M
    %               = [y_M(x_1), ..., y_M(x_N)]
    %           
    % 
    % De deelruimte waarin we zullen benaderen (en waarin y_LS dus
    % zal liggen, is 
    %
    %       De vectorruimte met GEEVALUEERDE (!) veeltermen 
    %       van graad<=M in de gegeven punten x_1, ... x_N.
    %       
    %
    % Als 
    %
    %   {e_0, ..., e_M} = {1, ..., x^M} 
    %
    % de basis is voor
    % de vectorruimte waarin we benaderen (zie hierboven), dan is
    %
    %   De Grammatrix G met G_ij = <e_i, e_j>_w 
    %   de HANKEL-matrix,
    %
    %   De matrix A met de evaluaties van e_0, ..., e_M in de
    %   punten x_1, ..., x_N de VANDERMONDE matrix.
    %
    %   => ZEER SLECHT GECONDITIONEERD PROBLEEM
    %   => ZEER INSTABIELE METHODE om KK probleem op te lossen
    %
    % We houden dit in het achterhoofd maar zullen toch deze methode
    % gebruiken om het probleem op te lossen.
    %
    % Een alternatieve methode is bvb. 
    %
    %   1. via de QR-factorisatie van A     A=QR
    %   2. via de SVD van A                 A=USV^T
    %
    %
    % Het oplossen van de NORMAALVGLen is EQUIVALENT aan het
    % OPLOSSEN VAN EEN OVERGEDETERMINEERD STELSEL
    %
    %   { c_0 + c_1 x + c_2 x^2 + ... + c_M x^M = f_1}
    %   ...
    %   { c_0 + c_1 x + c_2 x^2 + ... + c_M x^M = f_N}
    %
    %
    % Het is duidelijk dat dit stelsel pas overgedetermineerd is als
    % het aantal vergelijkingen (N) groter is dan het aantal onbekenden,
    % en dus het aantal coefficienten (M+1):
    %
    %       overgedetermineerd <=> N > M+1
    %
    % We gooien daarom een foutmelding <=>  N <= M+1
    
    % N is het aantal datapunten (x_1,f_1) ... (x_N, f_N)
    % dat we discreet willen benaderen
    N = length(x);

    % Aantal datapunten N moet minstens gelijk zijn aan het aantal 
    % coefficienten (M+1) om een unieke oplossing te hebben.
    %
    % We laten overgedetermineerd & interpolatie toe
    %
    % N > M+1  -> Overgedetermineerd (Kleinste Kwadraten)
    % N = M+1  -> Interpolatie (Uniek bepaald, residu = 0)
    % N < M+1  -> Ondergedetermineerd (Oneindig veel oplossingen) -> ERROR
    %
    % --------
    % OM OPGAVE 1B TE DOEN RUNNEN, VOEREN WE DEZE CHECK NIET UIT.
    % MATLAB ZAL NOG STEEDS EEN OPLOSSING TERUGGEVEN (OOK AL ZIJN ER
    % ONEINDIG VEEL OPLOSSINGEN). DAT IS DAN DE MINIMAL NORM SOLUTION
    % EN HET RESIDU ZAL DAN OOK (BIJNA) 0 ZIJN: 10e-12.
    % --------

    % Check of alle input geldig is
    if any([length(f),length(w)]~= length(x))
        error('plotres:input', ...
                'r and w should have the same number of elements as x');  
    end
    
    % =====================================================================
    % EFFICIËNTE OPBOUW VIA DE HANKEL EIGENSCHAP
    % =====================================================================

    % HOE WERKT EEN HANKEL MATRIX?
    %
    % De Gram-matrix G = A^T * W * A heeft afmeting (M+1) x (M+1)
    % %.
    % Het element op rij i en kolom j is het inproduct:
    %       G_ij = som( w_n * x_n^(i-1) * x_n^(j-1) ) 
    %            = som( w_n * x_n^(i+j-2) )
    %
    % We zien dat G_ij ENKEL afhangt van de totale macht (i+j-2).
    % Hierdoor is de matrix constant over de anti-diagonalen.
    %
    % Hoewel G in totaal (M+1)^2 elementen heeft, zijn er maar (2M + 1) 
    % UNIEKE sommen, lopend van macht 0 tot en met macht 2M.
    
    %%%%
    % STAP 1: BEREKEN ENKEL DE UNIEKE MACHTEN
    %
    % In plaats van de zware matrix A te vermenigvuldigen, berekenen we 
    % hier enkel de unieke machten x^0 t.e.m. x^2M voor alle N punten.
    %
    %   - x(:, ones(1,2*M)) kopieert de vector x gewoon 2M keer naast elkaar.
    %
    %   - cumprod(..., 2) neemt het product per rij (dimensie 2).
    %     Dit genereert in één klap: [x^1, x^2, x^3, ..., x^2M].
    %
    %   - ones(N, 1) voegt de kolom voor x^0 toe aan het begin.
    %
    % Resultaat: X is een N x (2M+1) matrix.
    %
    X = [ones(N, 1), cumprod(x(:,ones(1,2*M)), 2)];
    
    %%%%
    % STAP 2: BEREKEN DE (2M+1) UNIEKE SOMMEN
    %
    % We moeten deze machten nu vermenigvuldigen met de gewichten (w_n) en
    % optellen.
    %
    % De wiskundig trage manier is: W = diag(w) en dan A'*W*A. Dit bouwt
    % een gigantische matrix W vol met nullen in het geheugen op.
    % 
    % De slimme manier: we vermenigvuldigen de vector w' (1 x N) direct 
    % met matrix X (N x 2M+1).
    % 
    % Dit levert direct een vector wX (1 x 2M+1) met EXACT de 2M+1 unieke 
    % sommen die we nodig hebben, zonder W ooit op te bouwen!
    wX = w'*X;

    % --> (1 x N) * (N x (2M+1))
    % --> (1 x 2M+1)
    
    %%%%
    % STAP 3: VUL DE GRAM-MATRIX IN
    %
    % We hebben nu alle bouwstenen in de vector wX,
    % die de 2M+1 unieke sommen (= GEWOGEN INPRODUCTEN) bevat
    % uit de Grammiaan matrix.
    %
    % Merk op dat voor de Grammiaan matrix geldt
    %
    %   - G = A^T * W * A
    %   - Dimensie (M+1) x (M+1)
    % 
    % Waarbij
    %
    %   - A     =   [ |  ...    |   ]
    %               [a_1 ... a_{M+1}]
    %               [ |  ...    |   ]
    %
    %   - a_i   =   [ e_i(x_1) ]
    %               [    ...   ]
    %               [ e_i(x_N) ]
    %
    %
    % hankel(C, R) schikt deze elementen in een (M+1) x (M+1) vierkant.
    %   - C = de eerste kolom van de matrix (machten x^0 t.e.m. x^M)
    %   - R = de laatste rij van de matrix (machten x^M t.e.m. x^2M)
    %
    % hankel() voert GEEN BEREKENINGEN uit, maar legt de waarden uit wX
    % simpelweg op de juiste anti-diagonalen.
    %
    AtWA = hankel(wX(1:M+1), wX(M+1:end)');
    
    % Form nx1 vector b for the right side of the equation
    % Only n weighted sums are needed here
    %
    %   [w_1 f_1, ..., w_N f_N] * [e_0(x_1), ..., e_{M}(x_1)]
    %                             [          ...            ]
    %                             [e_0(x_N), ..., e_{M}(x_N)]
    %
    % Voor het rechterlid b = A^T * W * f hebben we geen machten tot 2M 
    % nodig, enkel de machten t.e.m. M (de basisveeltermen).
    %
    % We nemen dus de eerste M+1 kolommen van onze eerder gemaakte matrix X.
    % (w.*f)' is de efficiënte vector-vermenigvuldiging i.p.v. (W*f).
    %
    b = (w.*f)'*X(:,1:M+1);
    
    % Los het stelsel (A^T*W*A)c = (A^T*W*f) op via de backslash operator.
    % Omdat we AtWA expliciet hebben opgebouwd, lost dit exact het
    % normaalstelsel op
    c = AtWA \ b';
end
