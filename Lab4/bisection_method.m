function [xsolution,ysolution,iterations,xtab,xdif] = bisection_method(a,b,max_iterations,ytolerance,fun)
% a - lewa granica przedziału poszukiwań miejsca zerowego
% b - prawa granica przedziału poszukiwań miejsca zerowego
% max_iterations - maksymalna liczba iteracji działania metody bisekcji
% ytolerance - wartość abs(fun(xsolution)) powinna być mniejsza niż ytolerance
% fun - nazwa funkcji, której miejsce zerowe będzie wyznaczane
%
% xsolution - obliczone miejsce zerowe
% ysolution - wartość fun(xsolution)
% iterations - liczba iteracji wykonana w celu wyznaczenia xsolution
% xtab - wektor z kolejnymi kandydatami na miejsce zerowe, począwszy od xtab(1)= (a+b)/2
% xdiff - wektor wartości bezwzględnych z różnic pomiędzy i-tym oraz (i+1)-ym elementem wektora xtab; xdiff(1) = abs(xtab(2)-xtab(1));


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