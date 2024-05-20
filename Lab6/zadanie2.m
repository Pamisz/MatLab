function [country, source, degrees, y_original, y_movmean, y_approximation, mse] = zadanie2(energy)
% Głównym celem tej funkcji jest wyznaczenie aproksymacji wygładzonych danych o produkcji energii elektrycznej w wybranym kraju i z wybranego źródła energii.
% Wygładzenie danych wykonane jest poprzez wywołanie funkcji movmean.
% Wybór kraju i źródła energii należy określić poprzez nadanie w tej funkcji wartości zmiennym typu string: country, source.
% Dopuszczalne wartości tych zmiennych można sprawdzić poprzez sprawdzenie zawartości struktury energy zapisanej w pliku energy.mat.
% 
% energy - struktura danych wczytana z pliku energy.mat
% country - [String] nazwa kraju
% source  - [String] źródło energii
% degrees - wektor zawierający cztery stopnie wielomianu dla których wyznaczono aproksymację
% y_original - dane wejściowe, czyli pomiary produkcji energii zawarte w wektorze energy.(country).(source).EnergyProduction
% y_movmean - średnia 12-miesięczna danych wejściowych, y_movmean = movmean(y_original,[11,0]);
% y_approximation - tablica komórkowa przechowująca cztery wartości funkcji aproksymującej wygładzone dane wejściowe. y_approximation stanowi aproksymację stopnia degrees(i).
% mse - wektor o rozmiarze 4x1: mse(i) zawiera wartość błędu średniokwadratowego obliczonego dla aproksymacji stopnia degrees(i).
country = 'France';
source =  'Nuclear';
degrees = [1, 5, 10, 15];
y_original = [];
y_movmean = [];
y_approximation = {};
mse = [];

% Sprawdzenie dostępności danych
if isfield(energy, country) && isfield(energy.(country), source)
    % Przygotowanie danych do aproksymacji
    dates = energy.(country).(source).Dates;
    y_original = energy.(country).(source).EnergyProduction;
    y_movmean = movmean(y_original,[11,0]);

    x = linspace(-1,1,length(y_original))';

    % Pętla po wielomianach różnych stopni
    for i = 1:length(degrees)
        % Wyznaczanie współczynników wielomianu aproksymującego
        p = polyfit(x, y_movmean, degrees(i));
    
        % Wyznaczanie wartości aproksymowanej funkcji na podstawie wyznaczonego wielomianu
        y_approximation{i} = polyval(p, x);
    
        % Wyznaczanie błędu średniokwadratowego
        mse(i) = mean((y_movmean - y_approximation{i}).^2);
    end

    figure;
    subplot(2,1,1);
    hold on;
    plot(x, y_original, 'k-', 'DisplayName', 'Original Data', 'Color', 'R');
    plot(x, y_movmean, 'k-', 'DisplayName', 'Moving Average Data');
    for i = 1:length(degrees)
        plot(x, y_approximation{i}, 'DisplayName', ['Degree ', num2str(degrees(i))]);
    end
    hold off;
    xlabel('x');
    ylabel('Energy Production');
    title('Moving Average Energy Production Approximations');
    legend('Location', 'best');
    
    subplot(2,1,2);
    bar(mse);
    xlabel('Degree of Polynomial');
    ylabel('Mean Squared Error');
    title('Mean Squared Error for Different Degrees');
    set(gca, 'XTick', 1:length(degrees));
    set(gca, 'XTickLabel', degrees);
    saveas(gcf, 'zadanie2.png');

else
    disp(['Dane dla (country=', country, ') oraz (source=', source, ') nie są dostępne.']);
end

end