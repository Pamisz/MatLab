function [matrix_condition_numbers, max_coefficients_difference_1, max_coefficients_difference_2] = zadanie3()
% Zwracane są trzy wektory wierszowe:
% matrix_condition_numbers - współczynniki uwarunkowania badanych macierzy Vandermonde
% max_coefficients_difference_1 - maksymalna różnica między referencyjnymi a obliczonymi współczynnikami wielomianu,
%       gdy b zawiera wartości funkcji liniowej
% max_coefficients_difference_2 - maksymalna różnica między referencyjnymi a obliczonymi współczynnikami wielomianu,
%       gdy b zawiera zaburzone wartości funkcji liniowej
    N = 5:40;
    matrix_condition_numbers = zeros(1, length(N));
    for i = 1:length(N)
        V = vandermonde_matrix(N(i));
        matrix_condition_numbers(i) = cond(V);
    end

    a1 = randi([20,30]);
    max_coefficients_difference_1 = zeros(1, length(N));
    for i = 1:length(N)
        ni = N(i);
        V = vandermonde_matrix(ni);
        b = linspace(0,a1,ni)';
        reference_coefficients = [ 0; a1; zeros(ni-2,1) ]; 
        calculated_coefficients = V \ b;
        max_coefficients_difference_1(i) = max(abs(calculated_coefficients-reference_coefficients));
    end

    max_coefficients_difference_2 = zeros(1, length(N));
    for i = 1:length(N)
        ni = N(i);
        V = vandermonde_matrix(ni);
        b = linspace(0,a1,ni)' + rand(ni,1)*1e-10; 
        reference_coefficients = [ 0; a1; zeros(ni-2,1) ]; 
        calculated_coefficients = V \ b;
        max_coefficients_difference_2(i) = max(abs(calculated_coefficients-reference_coefficients));
    end


    figure;

    subplot(3,1,1);
    semilogy(N, matrix_condition_numbers, '-o');
    title('Vandermonde');
    xlabel('Rozmiar macierzy Vandermonde (N)');
    ylabel('Współczynnik uwarunkowania');
    grid on;

    subplot(3,1,2);
    plot(N, max_coefficients_difference_1, '-o');
    title('Wartości funkcji liniowej');
    xlabel('Rozmiar macierzy Vandermonde (N)');
    ylabel('Maksymalna różnica');
    grid on;

    subplot(3,1,3);
    plot(N, max_coefficients_difference_2, '-o');
    title('Zaburzone wartości funkcji liniowej');
    xlabel('Rozmiar macierzy Vandermonde (N)');
    ylabel('Maksymalna różnica');
    grid on;

    saveas(gcf, 'zadanie3.png');

for i = 1:length(N)
    ni = N(i);
    V = vandermonde_matrix(ni);
    
    b = linspace(0,a1,ni)' + rand(ni,1)*1e-10;
    reference_coefficients = [ 0; a1; zeros(ni-2,1) ]; 
    
    calculated_coefficients = V \ b;
    
    max_coefficients_difference_2(i) = max(abs(calculated_coefficients-reference_coefficients));
end

end

function V = vandermonde_matrix(N)
x = linspace(0, 1, N);
    V = zeros(N, N);
    for i = 1:N
        V(:, i) = x.^(i-1);
    end
end