a = 1;
b = 50;
ytolerance = 1e-12;
max_iterations = 100;

[time_bisection, ~, iterations_bisection, xtab_bisection, xdif_bisection] = bisection_method(a, b, max_iterations, ytolerance, @rocket_velocity);
[time_secant, ~, iterations_secant, xtab_secant, xdif_secant] = secant_method(a,b,max_iterations,ytolerance, @rocket_velocity);
figure

subplot(2,1,1)
plot(1:length(xtab_bisection), xtab_bisection, '-o', 'DisplayName', 'Bisekcja');
hold on;
plot(1:length(xtab_secant), xtab_secant, '-o', 'DisplayName', 'Sieczne');
hold off;
xlabel('Iteracja');
ylabel('Przybliżenie');
title('Zmiany w kolejnych przybliżeniach');
legend;
grid on;

subplot(2, 1, 2); 
x_bisection_diff = 1:length(xdif_bisection);
x_secant_diff = 1:length(xdif_secant);
semilogy(x_bisection_diff, xdif_bisection, '-o', 'DisplayName', 'Bisekcja');
hold on;
semilogy(x_secant_diff, xdif_secant, '-o', 'DisplayName', 'Sieczne');
hold off;
xlabel('Iteracja');
ylabel('Różnica (skala logarytmiczna)');
title('Różnice w kolejnych przybliżeniach');
legend;
grid on;

saveas(gcf, '4.6.png');

function velocity_delta = rocket_velocity(t)
M = 750; % [m/s]
g = 1.622; %[m/s^2]
m0 = 150000; %[kg]
u = 2000; %[m/s]
q = 2700; %[kg/s]

if t <= 0
    error("Nieprawidłowa wartość czasu.")
end

velocity_delta = u * log(m0/(m0-q*t))-g*t;
velocity_delta = velocity_delta - M;

end

function [xsolution,ysolution,iterations,xtab,xdif] = secant_method(a,b,max_iterations,ytolerance,fun)
xsolution = NaN;
ysolution = NaN;
iterations = 0;
xtab = [];
xdif = [];
    while iterations <= max_iterations
        iterations = iterations + 1;
        c = b - (fun(b)*(b-a)/(fun(b)-fun(a)));

        fc = abs(fun(c));

        xtab = [xtab;c];


        if length(xtab) > 1
            xdif = [xdif; abs(xtab(end) - xtab(end - 1))];
        end

        if fc < ytolerance || abs(b - a) < ytolerance
            xsolution = c;
            ysolution = fc;
            return;
        end
        a = b;
        b = c;
    end

    xsolution = c;
    ysolution = fc;
    warning('Nie osiągnięto wymaganego kryterium tolerancji przed zakończeniem iteracji.');
end

function [xsolution,ysolution,iterations,xtab,xdif] = bisection_method(a,b,max_iterations,ytolerance,fun)
    xtab = [];
    xdif = [];
    iterations = 0;
    xsolution = NaN;
    ysolution = NaN;

    if fun(a) * fun(b) > 0
        error('Funkcja musi mieć przeciwny znak na końcach przedziału.');
    end

    while iterations <= max_iterations
        iterations = iterations + 1;
        c = (a + b) / 2;

        fc = abs(fun(c));

        xtab = [xtab;c]


        if length(xtab) > 1
            xdif = [xdif; abs(xtab(end) - xtab(end - 1))];
        end

        if fc < ytolerance || abs(b - a) < ytolerance
            xsolution = c;
            ysolution = fc;
            return;
        else if fun(a) * fun(c) < 0
            b = c;
        else 
            a = c;
        end
   
    end

    xsolution = c;
    ysolution = fc;
    warning('Nie osiągnięto wymaganego kryterium tolerancji przed zakończeniem iteracji.');
    end
end