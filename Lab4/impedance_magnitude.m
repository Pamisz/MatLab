a = 1;
b = 50;
ytolerance = 1e-12;
max_iterations = 100;
[omega_bisection, ~, iterations_bisection, xtab_bisection, xdif_bisection] = bisection_method(a, b, max_iterations, ytolerance, @impedance_difference);
[omega_secant, ~, iterations_secant, xtab_secant, xdif_secant] = secant_method(a,b,max_iterations,ytolerance, @impedance_difference);
figure

subplot(2, 1, 1); 
plot(xtab_bisection, 'DisplayName', 'Bisekcja');
hold on;
plot(xtab_secant, 'DisplayName', 'Sieczne');
hold off;
xlabel('Iteracja');
ylabel('Przybliżenie');
title('Zmiany w kolejnych przybliżeniach');
legend;
grid on;

subplot(2, 1, 2); 
semilogy(xdif_bisection, 'DisplayName', 'Bisekcja');
hold on;
semilogy(xdif_secant, 'DisplayName', 'Sieczne');
hold off;
xlabel('Iteracja');
ylabel('Różnica (skala logarytmiczna)');
title('Różnice w kolejnych przybliżeniach');
legend;
grid on;
saveas(gcf, '4.4.png');




