function [A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Jacobi(N)
% A - macierz z równania macierzowego A * x = b
% b - wektor prawej strony równania macierzowego A * x = b
% M - macierz pomocnicza opisana w instrukcji do Laboratorium 3 – sprawdź wzór (5) w instrukcji, który definiuje M jako M_J.
% bm - wektor pomocniczy opisany w instrukcji do Laboratorium 3 – sprawdź wzór (5) w instrukcji, który definiuje bm jako b_{mJ}.
% x - rozwiązanie równania macierzowego
% err_norm - norma błędu residualnego wyznaczona dla rozwiązania x; err_norm = norm(A*x-b)
% time - czas wyznaczenia rozwiązania x
% iterations - liczba iteracji wykonana w procesie iteracyjnym metody Jacobiego
% index_number - Twój numer indeksu
index_number = 193249;
L1 = 9;

[A,b] = generate_matrix(N, L1);


L = tril(A, -1);
U = triu(A, 1);
D = diag(diag(A));

M = -(L + U) / D;
bm = D \ b;

x = ones(N, 1);
max_iterations = 1000; 
tolerance = 1e-12; 
tic;
for iterations = 1:max_iterations
    x_new = M * x + bm;
    residual = A * x_new - b;
    err_norm = norm(residual);
    if err_norm < tolerance
        break;
    end
    x = x_new;
end
time = toc;
end
