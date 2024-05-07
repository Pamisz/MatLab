import numpy as np
import matplotlib.pyplot as plt
import math
import time


def doting(vec1, vec2):
    if (len(vec1) != len(vec2)):
        raise ValueError("Wektory mają różną długość!")

    result = 0.0
    for i in range(len(vec1)):
        result += vec1[i] * vec2[i]
    return result


def jacobi(A, b, tolerance=1e-9, max_iterations=1000):
    N = len(b)
    x_prev = np.zeros_like(b)
    iter_count = 0
    res_norms = []

    start = time.time()
    for k in range(max_iterations):
        x = np.zeros_like(b)
        for i in range(N):
            s1 = doting(A[i, :i], x_prev[:i])
            s2 = doting(A[i, i + 1:], x_prev[i + 1:])
            x[i] = (b[i] - s1 - s2) / A[i, i]

        res = np.array([doting(A[j], x) for j in range(N)]) - b
        norm_res = np.linalg.norm(res, 2)
        res_norms.append(norm_res)

        if norm_res < tolerance:
            break

        x_prev = x.copy()
        iter_count += 1

    end = time.time()
    return x, iter_count, res_norms, (end - start)


def gauss_seidel(A, b, tolerance=1e-9, max_iterations=1000):
    N = len(b)
    x = np.zeros_like(b)
    iter_count = 0
    res_norms = []

    start = time.time()
    for k in range(max_iterations):
        for i in range(N):
            s1 = doting(A[i, :i], x[:i])
            s2 = doting(A[i, i + 1:], x[i + 1:])
            x[i] = (b[i] - s1 - s2) / A[i, i]

        res = np.array([doting(A[j], x) for j in range(N)]) - b
        norm_res = np.linalg.norm(res, 2)
        res_norms.append(norm_res)

        if norm_res < tolerance:
            break

        iter_count += 1

    end = time.time()
    return x, iter_count, res_norms, (end - start)


def lu_decomposition(A):
    n = len(A)
    L = [[0] * n for _ in range(n)]
    U = [[0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i, n):
            s = sum(L[i][k] * U[k][i] for k in range(i))
            U[i][j] = A[i][j] - s

        L[i][j] = 1
        for j in range(i + 1, n):
            s = sum(L[j][k] * U[k][i] for k in range(i))
            if U[i][j] == 0:
                continue
            L[j][i] = (A[j][i] - s) / U[i][j]

    return L, U


def forward_substitution(L, b):
    n = len(L)
    y = [0] * n

    for i in range(n):
        sum_ = sum(L[i][j] * y[j] for j in range(i))
        y[i] = (b[i] - sum_)

    return y


def backward_substitution(U, y):
    n = len(U)
    x = [0] * n

    for i in range(n - 1, -1, -1):
        sum_ = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - sum_) / U[i][i]

    return x


def solve_lu(A, b):
    L, U = lu_decomposition(A)
    y = forward_substitution(L, b)
    x = backward_substitution(U, y)
    return x


def matrix_vector_product(A, x):
    n = len(A)
    result = [0] * n
    for i in range(n):
        result[i] = sum(A[i][j] * x[j] for j in range(len(x)))
    return result


def norm_euclidean(vec):
    return math.sqrt(sum(x * x for x in vec))


#ZADANIE A
print("ZADANIE A")
N = 949
a1 = 7
a2 = -1
a3 = -1

A = np.zeros((N, N))
for i in range(N):
    A[i, i] = a1
    if i > 0:
        A[i, i - 1] = a2
        A[i - 1, i] = a2
    if i > 1:
        A[i, i - 2] = a3
        A[i - 2, i] = a3

b = np.array([np.sin(n * 4) for n in range(N)])
print("Stworzono układ równań.")

#ZADANIE B
print("\nZADANIE B")
x_jacobi, jacobi_iters, jacobi_res_norms, jacobi_time = jacobi(A, b)
x_gauss, gauss_iters, gauss_res_norms, gauss_time = gauss_seidel(A, b)

print("Metoda Jacobiego:\nLiczba iteracji: ", jacobi_iters, "\nCzas operacji: ", jacobi_time)
print("\nMetoda Gaussa-Seidla:\nLiczba iteracji: ", gauss_iters, "\nCzas operacji: ", gauss_time)

plt.figure()
plt.plot(jacobi_res_norms, label="Jacobi")
plt.plot(gauss_res_norms, label="Gauss-Seidel")
plt.yscale("log")
plt.xlabel("Iteracja")
plt.ylabel("Norma residuum")
plt.title("Zmiana normy residum w kolejnych iteracjach")
plt.legend()
plt.savefig("WYKRES B")
plt.show()

#ZADANIE C
print("\nZADANIE C")
N = 949
a1 = 3
a2 = -1
a3 = -1

A = np.zeros((N, N))
for i in range(N):
    A[i, i] = a1
    if i > 0:
        A[i, i - 1] = a2
        A[i - 1, i] = a2
    if i > 1:
        A[i, i - 2] = a3
        A[i - 2, i] = a3
        
b = np.array([np.sin(n * 4) for n in range(N)])

x_jacobi, jacobi_iters, jacobi_res_norms, jacobi_time = jacobi(A, b)
x_gauss, gauss_iters, gauss_res_norms, gauss_time = gauss_seidel(A, b)

print("\nMetoda Jacobiego:\nLiczba iteracji: ", jacobi_iters, "\nCzas operacji: ", jacobi_time)
print("\nMetoda Gaussa-Seidla:\nLiczba iteracji: ", gauss_iters, "\nCzas operacji: ", gauss_time)

plt.figure()
plt.plot(jacobi_res_norms, label="Jacobi")
plt.plot(gauss_res_norms, label="Gauss-Seidel")
plt.yscale("log")
plt.xlabel("Iteracja")
plt.ylabel("Norma residuum")
plt.title("Zmiana normy residum w kolejnych iteracjach")
plt.legend()
plt.savefig("WYKRES C")
plt.show()

#ZADANIE D
print("\nZADANIE D")
N = 949
a1 = 3
a2 = -1
a3 = -1

A = np.zeros((N, N))
for i in range(N):
    A[i, i] = a1
    if i > 0:
        A[i, i - 1] = a2
        A[i - 1, i] = a2
    if i > 1:
        A[i, i - 2] = a3
        A[i - 2, i] = a3

b = np.array([np.sin(n * 4) for n in range(N)])

solv = solve_lu(A, b)
Ax = matrix_vector_product(A, solv)
residuum = [Ax[i] - b[i] for i in range(len(b))]
norm_residuum = norm_euclidean(residuum)
print("Norma residuum metody bezpośredniej LU: ", norm_residuum)

#ZADANIE E
print("\nZADANIE E")
a1 = 7
Ns = [100, 500, 1000, 2000, 3000]
jacobi_times = []
gauss_times = []
lu_times = []

for n in Ns:
    A_temp = np.zeros((n, n))
    for i in range(n):
        A_temp[i, i] = a1
        if i > 0:
            A_temp[i, i - 1] = a2
            A_temp[i - 1, i] = a2
        if i > 1:
            A_temp[i, i - 2] = a3
            A_temp[i - 2, i] = a3

    b_temp = np.array([np.sin(k * 4) for k in range(n)])

    start_time = time.time()
    jacobi(A_temp, b_temp)
    jacobi_times.append(time.time() - start_time)

    start_time = time.time()
    gauss_seidel(A_temp, b_temp)
    gauss_times.append(time.time() - start_time)

    start_time = time.time()
    solve_lu(A_temp, b_temp)
    lu_times.append(time.time() - start_time)

plt.figure()
plt.plot(Ns, jacobi_times, label="Jacobi")
plt.plot(Ns, gauss_times, label="Gauss-Seidel")
plt.plot(Ns, lu_times, label="LU")
plt.xlabel("Rozmiar macierzy (N)")
plt.ylabel("Czas rozwiązania (s)")
plt.title("Porównanie czasu rozwiązania dla różnych metod")
plt.legend()
plt.savefig("WYKRES E")
plt.show()
