function plot_problem_5(N, time_Jacobi, time_Gauss_Seidel, iterations_Jacobi, iterations_Gauss_Seidel)
    subplot(2, 1, 1);
    plot(N, time_Jacobi, '-o', 'DisplayName', 'Jacobi');
    hold on;
    plot(N, time_Gauss_Seidel, '-o', 'DisplayName', 'Gauss-Seidel');
    hold off;
    title('Time complexity');
    xlabel('size');
    ylabel('time [s]');
    legend('Location', 'eastoutside');

    subplot(2, 1, 2);
    bar(N, [iterations_Jacobi', iterations_Gauss_Seidel']);
    title('Iterations');
    xlabel('size');
    ylabel('iterations');
    legend('Jacobi', 'Gauss-Seidel', 'Location', 'eastoutside');

    saveas(gcf, 'zadanie5.png');
end
