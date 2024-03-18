function [circles, index_number, circle_areas, rand_counts, counts_mean] = generate_circles(a, r_max, n_max)
index_number = 193249; % numer Twojego indeksu

    circles = zeros(3, n_max); 
    circle_areas = zeros(1,n_max);

    rand_counts = zeros(n_max, 1);
    counts_mean = zeros(n_max, 1);
    
    i = 1;
    k = n_max+1;
    while i ~= k
        R = rand * r_max + eps;
        X = rand * (a-2*R) + R + eps;
        Y = rand * (a-2*R) + R + eps;
        
        rand_counts(i) = rand_counts(i)+1; 

        intersects = false;
        for j = 1:i-1
            if norm([X; Y] - circles(2:3, j)) < R + circles(1, j)
                intersects = true;
                break;
            end
        end
        
        if ~intersects
            circles(:, i) = [R; X; Y];
            circle_area = pi * R^2;
            if i == 1
                circle_areas(i) = circle_area;
            else
                circle_areas(i) = circle_areas(i-1) + circle_area;
            end

            for j = 1:i
                counts_mean(i) = counts_mean(i) + rand_counts(j);
            end
            counts_mean(i) = counts_mean(i) ./ i;
            i = i+1;
        end
            
    end
end



