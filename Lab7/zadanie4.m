function [integration_error, Nt, ft_5, xr, yr, yrmax] = zadanie4()
    % Numeryczne całkowanie metodą Monte Carlo.
    %
    %   integration_error - wektor wierszowy. Każdy element integration_error(1,i)
    %       zawiera błąd całkowania obliczony dla liczby losowań równej Nt(1,i).
    %       Zakładając, że obliczona wartość całki dla Nt(1,i) próbek wynosi
    %       integration_result, błąd jest definiowany jako:
    %       integration_error(1,i) = abs(integration_result - reference_value),
    %       gdzie reference_value to wartość referencyjna całki.
    %
    %   Nt - wektor wierszowy zawierający liczby losowań, dla których obliczano
    %       wektor błędów całkowania integration_error.
    %
    %   ft_5 - gęstość funkcji prawdopodobieństwa dla n=5
    %
    %   [xr, yr] - tablice komórkowe zawierające informacje o wylosowanych punktach.
    %       Tablice te mają rozmiar [1, length(Nt)]. W komórkach xr{1,i} oraz yr{1,i}
    %       zawarte są współrzędne x oraz y wszystkich punktów zastosowanych
    %       do obliczenia całki przy losowaniu Nt(1,i) punktów.
    %
    %   yrmax - maksymalna dopuszczalna wartość współrzędnej y losowanych punktów

    Nt = 5:50:10^4;
    reference_value = 0.0473612919396179; % wartość referencyjna całki
    integration_error = zeros(1, length(Nt));
    xr = cell(1, length(Nt));
    yr = cell(1, length(Nt));

    ft_5 = calculate_ft(5);
    yrmax = ft_5 + 0.01;

    for i = 1:length(Nt)
        [integration_result, xr{1,i}, yr{1,i}] = monte_carlo_integration(@calculate_ft, Nt(1,i), yrmax);
        integration_error(1,i) = abs(integration_result - reference_value);
    end

    loglog(Nt, integration_error);
    xlabel('Number of random points (N)');
    ylabel('Integration error');
    title('Montle Carlo');
    saveas(gcf, 'zadanie4.png');
end

function ft = calculate_ft(t)
    mu = 10; 
    sigma = 3; 
    ft = (1 / (sigma * sqrt(2 * pi))) * exp(-((t - mu)^2) / (2 * sigma^2));
end

function [integral, xr, yr] = monte_carlo_integration(f, N, yrmax)
    a = 0; 
    b = 5; 
    S = (b - a) * yrmax;
    
    xr = a + (b - a) * rand(1, N);
    yr = yrmax * rand(1, N);
    fx = zeros(1, N);
    N1 = 0;
    for i = 1:N
        fx(i) = f(xr(i));
        if yr(i) < fx(i)
            N1 = N1 + 1;
        end
    end

    integral = N1/N * S;
end