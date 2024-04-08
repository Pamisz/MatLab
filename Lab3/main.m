%%               3.1
N = 100;
[A,b,x,time_direct,err_norm,index_number] = solve_direct(N);
%%               3.2
N = 1000:1000:8000;
n = length(N);
vtime_direct = ones(1,n);

for i = 1:n
    tic
    solve_direct(N(i));
    vtime_direct(i) = toc;
end
plot_direct(N,vtime_direct);
%%                3.3
N = 100;
[A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Jacobi(N);
%%                3.4
N = 100;
[A,b,M,bm,x,err_norm,time,iterations,index_number] = solve_Gauss_Seidel(N);
%%                3.5
N = 1000:1000:8000;
n = length(N);
time_Jacobi = ones(1,n);
time_Gauss_Seidel = 2*ones(1,n);
iterations_Jacobi = 40*ones(1,n);
iterations_Gauss_Seidel = 40*ones(1,n);

for i = 1:n
    [~,~,~,~,~,~,time_Jacobi(i),iterations_Jacobi(i),~] = solve_Jacobi(N(i));
    [~,~,~,~,~,~,time_Gauss_Seidel(i),iterations_Gauss_Seidel(i),~] = solve_Gauss_Seidel(N(i));

end
plot_problem_5(N,time_Jacobi,time_Gauss_Seidel,iterations_Jacobi,iterations_Gauss_Seidel);
%% DIRECT

data = load('filtr_dielektryczny.mat');
A = data.A;
b = data.b;

x = A\b;
err_norm_direct = norm(A* x - b)
%% JACOBI
data = load('filtr_dielektryczny.mat');
A = data.A;
b = data.b;

L = tril(A, -1);
U = triu(A, 1);
D = diag(diag(A));

M = -inv(D) * (L + U);
bm = D \ b;

x = ones(size(b));
max_iterations = 1000; 
tolerance = 1e-12; 
errors = zeros(1, max_iterations);
for iterations = 1:max_iterations
    x_new = M * x + bm;
    residual = A * x_new - b;
    err_norm_jacobi = norm(residual);
     errors(iterations) = err_norm_jacobi;
    if err_norm_jacobi < tolerance
        break;
    end
    x = x_new;
end
err_norm_jacobi

plot(1:iterations, errors(1:iterations));
xlabel('Iteracje');
ylabel('Błąd');
title('Zbieżność metody Jacobiego');
saveas(gcf, 'Jacobi.png');
%% GAUSS
data = load('filtr_dielektryczny.mat');
A = data.A;
b = data.b;

M = -(D + L) \ U;
bm = (D + L) \ b;

x = ones(size(b));
tolerance = 1e-12; 
errors = zeros(1, max_iterations);
for iterations = 1:max_iterations
    x_new = M * x + bm;
    residual = A * x_new - b;
    err_norm_gauss_seidel = norm(residual);
    errors(iterations) = err_norm_gauss_seidel;
    if err_norm_gauss_seidel < tolerance
        break;
    end
    x = x_new;
end
err_norm_gauss_seidel
plot(1:iterations, errors(1:iterations));
xlabel('Iteracje');
ylabel('Błąd');
title('Zbieżność metody Gaussa');
saveas(gcf, 'Gauss.png');

