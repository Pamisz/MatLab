import math

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def lagrange_basis(x, i, x_points):
    result = 1
    for j in range(len(x_points)):
        if j != i:
            result *= (x - x_points[j]) / (x_points[i] - x_points[j])
    return result


def lagrange_interpolation(x_points, y_points, x_values):
    interpolated_values = []
    for x in x_values:
        interpolated_value = 0
        for i in range(len(x_points)):
            interpolated_value += y_points[i] * lagrange_basis(x, i, x_points)
        interpolated_values.append(interpolated_value)
    return interpolated_values


def get_Chebyshev_nodes(N, points, heights):
    k = np.arange(N)
    chebyshev_nodes = np.cos((np.pi * k) / (N - 1))
    scaled_nodes = 0.5 * (chebyshev_nodes + 1) * (len(points) - 1)
    indices = np.round(scaled_nodes).astype(int)
    nodes = [points[i] for i in indices]
    values = [heights[i] for i in indices]

    return nodes, values


def get_nodes(N, points, heights):
    if N > len(points):
        raise ValueError("N cannot be greater than the length of the points list")

    if N == 1:
        return [points[0]], [heights[0]]

    indices = np.linspace(0, len(points) - 1, N, dtype=int)
    nodes = [points[i] for i in indices]
    selected_heights = [heights[i] for i in indices]

    return nodes, selected_heights


def tridiagonal_solver(a, b, c, d):
    n = len(d)
    w = np.zeros(n - 1, float)
    g = np.zeros(n, float)
    p = np.zeros(n, float)

    w[0] = c[0] / b[0]
    g[0] = d[0] / b[0]

    for i in range(1, n - 1):
        w[i] = c[i] / (b[i] - a[i - 1] * w[i - 1])

    for i in range(1, n):
        g[i] = (d[i] - a[i - 1] * g[i - 1]) / (b[i] - a[i - 1] * w[i - 1])

    p[-1] = g[-1]

    for i in range(n - 2, -1, -1):
        p[i] = g[i] - w[i] * p[i + 1]

    return p


def cubic_spline_coefficients(x, y):
    n = len(x) - 1
    h = [x[i + 1] - x[i] for i in range(n)]

    alpha = [0] * n
    for i in range(1, n):
        alpha[i] = (3 / h[i] * (y[i + 1] - y[i]) - 3 / h[i - 1] * (y[i] - y[i - 1]))

    l = [1] + [0] * n
    mu = [0] * n
    z = [0] * (n + 1)

    for i in range(1, n):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1]
        mu[i] = h[i] / l[i]
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i]

    l[-1] = 1
    z[-1] = 0

    b = [0] * n
    c = [0] * (n + 1)
    d = [0] * n

    for j in range(n - 1, -1, -1):
        c[j] = z[j] - mu[j] * c[j + 1]
        b[j] = (y[j + 1] - y[j]) / h[j] - h[j] * (c[j + 1] + 2 * c[j]) / 3
        d[j] = (c[j + 1] - c[j]) / (3 * h[j])

    return y[:-1], b, c[:-1], d


def cubic_spline_interpolation(x_points, y_points, x_values, linspace):
    a, b, c, d = cubic_spline_coefficients(x_points, y_points)
    n = len(x_points) - 1
    interpolated_values = []

    for x in x_values:
        for i in range(n):
            if linspace:
                if x_points[i] <= x <= x_points[i + 1]:
                    dx = x - x_points[i]
                    interpolated_values.append(a[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3)
                    break
            else:
                if x_points[i] >= x >= x_points[i + 1]:
                    dx = x - x_points[i]
                    interpolated_values.append(a[i] + b[i] * dx + c[i] * dx ** 2 + d[i] * dx ** 3)
                    break

    return interpolated_values

def draw_plot(title, linspace, points, heights, larange, plot_name):
    plt.figure(figsize=(10, 8))
    plt.ylabel("Height (m)")
    plt.xlabel("Distance (m)")
    plt.title(title)
    plt.plot(points, heights, label="Original data ", color="blue")
    x = 0
    y = 0
    interpolated_values = 0
    size = [5, 10, 15, 20]
    for s in size:
        if linspace:
            x, y = get_nodes(s, points, heights)
            plt.scatter(x, y, marker='o', label=f"Linspace nodes - {s}")
        else:
            x, y = get_Chebyshev_nodes(s, points, heights)
            plt.scatter(x, y, marker='o', label=f"Chebyshev nodes - {s}")
        if larange:
            interpolated_values = lagrange_interpolation(x, y, points)
        else:
            interpolated_values = cubic_spline_interpolation(x, y, points, linspace)

        plt.plot(points, interpolated_values, label=f"Interpolation for {s} nodes")

    plt.legend()
    plt.savefig(plot_name, format='png')
    #plt.show()


def draw_plots_for_tour(path, tour_nr):
    file = pd.read_csv(path)
    points = file["Dystans (m)"]
    heights = file["Wysokość (m)"]
    draw_plot(f"Lagrange method, linspace nodes, tour {tour_nr}", 1, points, heights, 1, f"Plot{tour_nr}a")
    draw_plot(f"Lagrange method, Chebyshev nodes,tour {tour_nr}", 0, points, heights, 1, f"Plot{tour_nr}b")
    draw_plot(f"Spline method, linspace nodes, tour {tour_nr}", 1, points, heights, 0, f"Plot{tour_nr}c")
    draw_plot(f"Spline method, Chebyshev nodes, tour {tour_nr}", 0, points, heights, 0, f"Plot{tour_nr}d")


#main
draw_plots_for_tour('2018_paths/SpacerniakGdansk.csv', 1)
draw_plots_for_tour('2018_paths/WielkiKanionKolorado.csv', 2)
draw_plots_for_tour('2018_paths/MountEverest.csv', 3)
draw_plots_for_tour('2018_paths/GlebiaChallengera.csv', 4)
draw_plots_for_tour('2018_paths/Unsyncable_ride.csv', 5)



