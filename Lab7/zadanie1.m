function [integration_error, Nt, ft_5, integral_1000] = zadanie1()
    % Numeryczne całkowanie metodą prostokątów.
    % Nt - wektor zawierający liczby podprzedziałów całkowania
    % integration_error - integration_error(1,i) zawiera błąd całkowania wyznaczony
    %   dla liczby podprzedziałów równej Nt(i). Zakładając, że obliczona wartość całki
    %   dla Nt(i) liczby podprzedziałów całkowania wyniosła integration_result,
    %   to integration_error(1,i) = abs(integration_result - reference_value),
    %   gdzie reference_value jest wartością referencyjną całki.
    % ft_5 - gęstość funkcji prawdopodobieństwa dla n=5
    
    reference_value = 0.0473612919396179; % Reference value of the integral

    Nt = 5:50:10^4; 
    integration_error = zeros(1, length(Nt)); 

    ft_5 = calculate_ft(5);

    integral_1000 = rectangle_integration(@calculate_ft, 1000);

    for i = 1:length(Nt)
        integration_result = rectangle_integration(@calculate_ft, Nt(i));
        integration_error(i) = abs(integration_result - reference_value);
    end

    % Plotting
    loglog(Nt, integration_error);
    xlabel('Number of subintervals (N)');
    ylabel('Integration error');
    title('Rectangle method');
    saveas(gcf, 'zadanie1.png');
end

function ft = calculate_ft(t)
    mu = 10; 
    sigma = 3; 
    ft = (1 / (sigma * sqrt(2 * pi))) * exp(-((t - mu)^2) / (2 * sigma^2));
end

function integral = rectangle_integration(f, N)
    a = 0; 
    b = 5; 
    h = (b - a) / N; 
    integral = 0;
    for i = 1:N
        x = a + (i-1)*h;
        x_next = a + i * h;
        integral = integral + f((x+x_next)/2)*h;
    end
end