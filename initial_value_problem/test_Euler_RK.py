from Euler_RK import *
import matplotlib.pyplot as plt

# ======================================================
# Equations
# ======================================================

# Linear system

# y'' - 2y' + y = 0
# y(0) = 5
# y'(0) = 7
def f1(x, y):
    return np.array([y[1], 2 * y[1] - y[0]])

# Nonlinear system

# y''y'y = (3cos3x - x^3/3 + 5x + 1)(-9sin3x - x^2 + 5)(-27cos3x - 2x)
# y(0) = 4
# y'(0) = 5
def f2(x, y):
    a = 3 * np.cos(3 * x) - x**3 / 3 + 5 * x + 1
    b = -9 * np.sin(3 * x) - x**2 + 5
    c = -27 * np.cos(3 * x) - 2 * x
    return np.array([y[1], a * b * c / y[1] / y[0]])

# ======================================================
# Exact solutions
# ======================================================

def sol1(x):
    return 5 * np.exp(x) + 2 * x * np.exp(x)

def sol2(x):
    return 3 * np.cos(3 * x) - x**3 / 3 + 5 * x + 1

# ======================================================
# Test function
# ======================================================
def test(f, sol, x, y0):
    
    y_1 = euler_method(f, x, y0)
    y_2 = runge_kutta_2_method(f, x, y0)
    y_3 = runge_kutta_4_method(f, x, y0)

    plt.scatter(x, sol(x), label = "exact")
    plt.plot(x, y_1[0], label = "Eiler")
    plt.plot(x, y_2[0], label = "RK2")
    plt.plot(x, y_3[0], label = "RK4")
    
    plt.legend()
    plt.show()

    return 0

# ======================================================
# Tests
# ======================================================

def main():
    NODES = 10
    
    x1 = np.linspace(0, 1, NODES)
    test(f1, sol1, x1, (5, 7))
    
    x2 = np.linspace(0, 1, NODES)
    test(f2, sol2, x2, (4, 5))

main()
