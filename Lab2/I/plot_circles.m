function plot_circles(a, circles, index_number)
    rectangle('Position', [0, 0, a, a], 'EdgeColor', 'b', 'LineWidth', 2); 
    hold on;
    for i  = 1:length(circles)
        R = circles(1,i);
        X = circles(2,i);
        Y = circles(3,i);
        plot_circle(R, X, Y);
    end
    hold off;
    axis equal;
    axis([0 a 0 a])
    grid on;
    print -dpng zadanie1.png
end