function [lake_volume, x, y, z, zmin] = zadanie5()
    % Funkcja zadanie5 wyznacza objętość jeziora metodą Monte Carlo.
    %
    %   lake_volume - objętość jeziora wyznaczona metodą Monte Carlo
    %
    %   x - wektor wierszowy, który zawiera współrzędne x wszystkich punktów
    %       wylosowanych w tej funkcji w celu wyznaczenia obliczanej całki.
    %
    %   y - wektor wierszowy, który zawiera współrzędne y wszystkich punktów
    %       wylosowanych w tej funkcji w celu wyznaczenia obliczanej całki.
    %
    %   z - wektor wierszowy, który zawiera współrzędne z wszystkich punktów
    %       wylosowanych w tej funkcji w celu wyznaczenia obliczanej całki.
    %
    %   zmin - minimalna dopuszczalna wartość współrzędnej z losowanych punktów

    N = 1e6; 
    x = 100 * rand(1, N); 
    y = 100 * rand(1, N); 
    zmin = -50; 
    z = zmin + (0 - zmin) * rand(1, N); 

    points_inside = 0; 
    for i = 1:N
        if z(i) > get_lake_depth(x(i), y(i))
            points_inside = points_inside + 1; 
        end
    end

    volume_bound = 100 * 100 * 50;
    lake_volume = points_inside / N * volume_bound
end