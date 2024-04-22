function [xsolution,ysolution,iterations,xtab,xdif] = secant_method(a,b,max_iterations,ytolerance,fun)
% a - lewa granica przedziału poszukiwań miejsca zerowego (x0=a)
% b - prawa granica przedziału poszukiwań miejsca zerowego (x1=b)
% max_iterations - maksymalna liczba iteracji działania metody siecznych
% ytolerance - wartość abs(fun(xsolution)) powinna być mniejsza niż ytolerance
% fun - nazwa funkcji, której miejsce zerowe będzie wyznaczane
%
% xsolution - obliczone miejsce zerowe
% ysolution - wartość fun(xsolution)
% iterations - liczba iteracji wykonana w celu wyznaczenia xsolution
% xtab - wektor z kolejnymi kandydatami na miejsce zerowe, począwszy od x2
% xdiff - wektor wartości bezwzględnych z różnic pomiędzy i-tym oraz (i+1)-ym elementem wektora xtab; xdiff(1) = abs(xtab(2)-xtab(1));

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