a = 1;
b = 60000;
ytolerance = 1e-12;
max_iterations = 100;

[n_bisection, ~, iterations_bisection, xtab_bisection, xdif_bisection] = bisection_method(a, b, max_iterations, ytolerance, @estimate_execution_time);
[n_secant, ~, iterations_secant, xtab_secant, xdif_secant] = secant_method(a,b,max_iterations,ytolerance, @estimate_execution_time);

figure

subplot(2,1,1)
plot(1:length(xtab_bisection), xtab_bisection, 'DisplayName', 'Bisekcja');
hold on;
plot(1:length(xtab_secant), xtab_secant, 'DisplayName', 'Sieczne');
hold off;
xlabel('Iteracja');
ylabel('Przybliżenie');
title('Zmiany w kolejnych przybliżeniach');
legend;
grid on;

subplot(2, 1, 2); 
x_bisection_diff = 1:length(xdif_bisection);
x_secant_diff = 1:length(xdif_secant);
semilogy(x_bisection_diff, xdif_bisection, 'DisplayName', 'Bisekcja');
hold on;
semilogy(x_secant_diff, xdif_secant, 'DisplayName', 'Sieczne');
hold off;
xlabel('Iteracja');
ylabel('Różnica (skala logarytmiczna)');
title('Różnice w kolejnych przybliżeniach');
legend;
grid on;
saveas(gcf, '4.8.png');