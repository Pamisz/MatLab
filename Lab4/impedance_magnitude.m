a = 1;
b = 50;
ytolerance = 1e-12;
max_iterations = 100;
[omega_bisection, ~, iterations_bisection, xtab_bisection, xdif_bisection] = bisection_method(a, b, max_iterations, ytolerance, @impedance_difference);
[omega_secant, ~, iterations_secant, xtab_secant, xdif_secant] = secant_method(a,b,max_iterations,ytolerance, @impedance_difference);
figure

subplot(2, 1, 1); 

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
saveas(gcf, '4.4.png');


function impedance_delta = impedance_difference(omega)
    R = 525;  
    C = 7e-5; 
    L = 3;    
    M = 75;

    if omega <= 0
        error('Omega musi być większa od zera.');
    end

    impedance = 1 / sqrt(1/(R^2) + (omega*C - 1/(omega*L))^2);

    impedance_delta = impedance - M;
end

function [xsolution,ysolution,iterations,xtab,xdif] = bisection_method(a,b,max_iterations,ytolerance,fun)

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
