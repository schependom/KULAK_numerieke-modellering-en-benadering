%%% Test difference between kkb2 and kkb2simple

%% Generate test data
% De uit te testen waarden voor N
%
% Ns = [10, 20, ..., 990, 1000]
%    = [10  20  ...  990  1000]
%
Ns = 10:10:1000;

% Hou de timings bij voor KKB2 en KKB2Simple
t1=[];
t2=[];

% Loop over N datapunten
for N = Ns

    % Maak N datapunten (x_n, f_n) \in R^N
    x = linspace(-1,1,N); 
    f = randn(1,N); 

    % Maak N gewichten \in R^N
    w = rand(1,N);
    
    %
    % kkb2 (De efficiënte Hankel-methode)
    %

    % Start timing
    tic;

    % Average over 100 runs
    for k = 1:100
        c1 = kkb2(x,f,w,6);
        %
        % c1 \in R^{ (M+1) x 1 }
        %
        %    --> coefficienten van y_M(x)
        %    --> y_M(x) = benaderende veelterm van graad M
        %    -->        = \sum_m^M c_m x^m
    end

    % Log time
    t1 = [t1, toc];
    
    %
    % kkb2simple (De brute-force methode)
    %

    % Start timing
    tic;

    % Average over 100 runs
    for k = 1:100
        c2 = kkb2simple(x,f,w,6);
    end

    % Log time
    t2 = [t2, toc];

end

%% Plot comparison figure

close all;
figure;
hold on;

% Gebruik semilogy voor een logaritmische y-as (runtijd groeit exponentieel)
semilogy(Ns, t1, 'b-', 'LineWidth', 1.5);
semilogy(Ns, t2, 'r-', 'LineWidth', 1.5);

legend('kkb2 (Hankel)', 'kkb2simple (Brute-force)', 'Location', 'northwest');
xlabel('Aantal datapunten (N)');
ylabel('Tijd in seconden (Log-schaal)');
title('Performantievergelijking: KKB2 vs KKB2Simple');
grid on;

% =========================================================================
% CONCLUSIE VAN DEZE TEST:
% =========================================================================
%
% Waarom snijden deze twee curves elkaar op de plot?
% 
% 1. Lage N (KKB2Simple is sneller): 
%
%    Voor kleine datasets is de brute-force wiskunde (diag(w) en A'*W*A) 
%    zodanig klein dat de processor dit in hardware vrijwel instantaan 
%    oplost. De "slimme" Hankel-methode (KKB2) is hier trager door 
%    'overhead': het aanroepen van functies als cumprod() en hankel() 
%    kost relatief meer administratietijd dan het eigenlijke rekenwerk.
%
% 2. Hoge N (KKB2 is exponentieel sneller):
%
%    Naarmate N groeit, ontploft het geheugengebruik van KKB2Simple 
%    (asymptotische complexiteit O(N^2)). De diag(w) matrix vult zich 
%    met miljoenen nullen waar de computer zinloos mee gaat rekenen. 
%    KKB2 schaalt veel vlakker (O(N)) omdat het met w'*X enkel de strikt 
%    noodzakelijke 2M+1 unieke elementen berekent en vector-vermenigvuldiging
%    gebruikt. Geen overbodige nullen, geen verspild geheugen.
% =========================================================================