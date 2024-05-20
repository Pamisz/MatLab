function [country, source, degrees, x_coarse, x_fine, y_original, y_yearly, y_approximation, mse, msek] = zadanie4(energy)
% Głównym celem tej funkcji jest wyznaczenie danych na potrzeby analizy dokładności aproksymacji wielomianowej.
% 
% energy - struktura danych wczytana z pliku energy.mat
% country - [String] nazwa kraju
% source  - [String] źródło energii
% x_coarse - wartości x danych aproksymowanych
% x_fine - wartości, w których wyznaczone zostaną wartości funkcji aproksymującej
% y_original - dane wejściowe, czyli pomiary produkcji energii zawarte w wektorze energy.(country).(source).EnergyProduction
% y_yearly - wektor danych rocznych
% y_approximation - tablica komórkowa przechowująca wartości nmax funkcji aproksymujących dane roczne.
%   - nmax = length(y_yearly)-1
%   - y_approximation{i} stanowi aproksymację stopnia i
%   - y_approximation{i} stanowi wartości funkcji aproksymującej w punktach x_fine
% mse - wektor mający nmax wierszy: mse(i) zawiera wartość błędu średniokwadratowego obliczonego dla aproksymacji stopnia i.
%   - mse liczony jest dla aproksymacji wyznaczonej dla wektora x_coarse
% msek - wektor mający (nmax-1) wierszy: msek zawiera wartości błędów różnicowych zdefiniowanych w treści zadania 4
%   - msek(i) porównuj aproksymacje wyznaczone dla i-tego oraz (i+1) stopnia wielomianu
%   - msek liczony jest dla aproksymacji wyznaczonych dla wektora x_fine

country = 'France';
source = 'Nuclear';
degrees = 1:4;
x_coarse = [];
x_fine = [];
y_original = [];
y_yearly = [];
y_approximation = [];
mse = [];
msek = [];

% Sprawdzenie dostępności danych
if isfield(energy, country) && isfield(energy.(country), source)
    % Przygotowanie danych do aproksymacji
    dates = energy.(country).(source).Dates;
    y_original = energy.(country).(source).EnergyProduction;

    % Obliczenie danych rocznych
    n_years = floor(length(y_original) / 12);
    y_cut = y_original(end-12*n_years+1:end);
    y4sum = reshape(y_cut, [12 n_years]);
    y_yearly = sum(y4sum,1)';


    N = length(y_yearly);
    P = (N-1)*10+1;
    x_coarse = linspace(-1, 1, N)';
    x_fine = linspace(-1, 1, P)';
    
    mse = zeros(N-1, 1);
    msek = zeros(N-2, 1);
    % Pętla po wielomianach różnych stopni
      for i = 1:N-1
            % Wykonanie aproksymacji wielomianowej
            p = my_polyfit(x_coarse, y_yearly, i);
            % Obliczenie wartości aproksymacji dla x_coarse
            y_approximation{i} = polyval(p, x_fine);
            % Obliczenie błędu średniokwadratowego
            mse(i) = mean((y_yearly - polyval(p, x_coarse)).^2);
      end
      
      for i = 1:N-2
        % Obliczenie błędu różnicowego
        msek(i) = mean((y_approximation{i} - y_approximation{i+1}).^2);
      end
      % Tworzenie trzech wykresów w jednym oknie
    figure;
    
    % Pierwszy wykres
    subplot(3, 1, 1);
    plot(x_coarse, y_yearly, 'b-', 'LineWidth', 1.5); % Dane oryginalne
    hold on;
    for i = 1:4
        plot(x_fine, y_approximation{degrees(i)}, '--', 'LineWidth', 1.2); % Aproksymacje dla wybranych stopni
    end
    hold off;
    xlabel('x');
    ylabel('y');
    title('Dane oryginalne oraz aproksymacje wielomianowe');
    legend('Dane oryginalne', 'Stopień 1', 'Stopień 2', 'Stopień 3', 'Stopień 4', 'Location', 'best');
    
    % Drugi wykres
    subplot(3, 1, 2);
    semilogy(1:length(mse), mse, 'r-', 'LineWidth', 1.5); % Błędy średniokwadratowe
    xlabel('Stopień wielomianu');
    ylabel('Błąd średniokwadratowy');
    title('Błąd średniokwadratowy w funkcji stopnia wielomianu');
    
    % Trzeci wykres
    subplot(3, 1, 3);
    semilogy(1:length(msek), msek, 'g-', 'LineWidth', 1.5); % Błędy różnicowe
    xlabel('Indeks');
    ylabel('Błąd różnicowy');
    title('Zbieżność funkcji aproksymujących');
    
    % Zapisanie wykresów do pliku
    saveas(gcf, 'zadanie4.png');

else
    disp(['Dane dla (country=', country, ') oraz (source=', source, ') nie są dostępne.']);
end

end

function p = my_polyfit(x, y, deg)
    % Obliczanie macierzy Vandermonde'a
    A = zeros(length(x), deg + 1);
    for i = 1:length(x)
        for j = 1:(deg + 1)
            A(i, j) = x(i)^(deg + 1 - j);
        end
    end
    
    % Obliczanie współczynników wielomianu
    p = (A' * A) \ (A' * y);
end