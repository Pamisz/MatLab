function [integration_error, Nt, ft_5, integral_1000] = zadanie2()
    % Numeryczne całkowanie metodą trapezów.
    % Nt - wektor zawierający liczby podprzedziałów całkowania
    % integration_error - integration_error(1,i) zawiera błąd całkowania wyznaczony
    %   dla liczby podprzedziałów równej Nt(i). Zakładając, że obliczona wartość całki
    %   dla Nt(i) liczby podprzedziałów całkowania wyniosła integration_result,
    %   to integration_error(1,i) = abs(integration_result - reference_value),
    %   gdzie reference_value jest wartością referencyjną całki.
    % ft_5 - gęstość funkcji prawdopodobieństwa dla n=5
    % integral_1000 - całka od 0 do 5 funkcji gęstości prawdopodobieństwa
    %   dla 1000 podprzedziałów całkowania
    
    Nt = 5:50:10^4;
    reference_value = 0.0473612919396179; % wartość referencyjna całki
    integration_error = zeros(size(Nt));

    ft_5 = calculate_ft(5);

    integral_1000 = trapezoidal_integration(@calculate_ft, 1000);

    for i = 1:length(Nt)
        integration_result = trapezoidal_integration(@calculate_ft, Nt(i));
       
        integration_error(i) = abs(integration_result - reference_value);
    end

    loglog(Nt, integration_error);
    xlabel('Number of subintervals (N)');
    ylabel('Integration error');
    title('Trapezoidal method');
    saveas(gcf, 'zadanie2.png');

end

function ft = calculate_ft(t)
    mu = 10; 
    sigma = 3; 
    ft = (1 / (sigma * sqrt(2 * pi))) * exp(-((t - mu)^2) / (2 * sigma^2));
end

function integral = trapezoidal_integration(f, N)
    a = 0; 
    b = 5;
    h = (b - a) / N; 
    integral = 0;
    for i = 1:N
        x = a + (i-1)*h;
        x_next = a + i * h;
        integral = integral + ((f(x)+f(x_next))/2)*h;
    end
end