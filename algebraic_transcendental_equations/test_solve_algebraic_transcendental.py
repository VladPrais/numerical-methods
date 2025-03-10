from solve_algebraic_transcendental import *
import matplotlib.pyplot as plt
import numpy as np

# ======================================================
# Equations
# ======================================================

# transcendental equations

# x - cos(x^2) = 0
def f1(x):
    return x - np.cos(x**2)

# for iteration method
# x = cos(x^2)
def u1(x):
    return np.cos(x**2)

# x^2 - e^(x^3) = 0
def f2(x):
    return 3 * x - np.exp(x**3)

# for iteration method
# x = (e^(x^3))^0.5
def u2(x):
    return np.exp(x**3) / 3

# algebraic equations

# x^2 - 2x + 1 = 0
def f3(x):
    return x**2 - 2 * x + 1

# for iteration method
# x = (1 + x^2) / 2
def u3(x):
    return (1 + x**2) / 2

# ======================================================
# Test function
# ======================================================
def test(f, u, x0, a, b):
    
    x = np.linspace(a, b)
    
    plt.plot(x, x)
    plt.plot(x, u(x))

    s1 = iteration_method(x0, u)
    s2 = newton_method(x0, f)
    s3 = false_position_method(x0, x0 + 0.5, f)

    plt.scatter(s1, s1, label = "iteration method")
    plt.scatter(s2, s2, label = "Newton method")
    plt.scatter(s3, s3, label = "false position method\n")
    
    plt.legend()
    plt.show()

    return 0

# ======================================================
# Tests
# ======================================================

def main():

    print("\nTest1")
    test(f1, u1, 0, 0, 2)
    print("\nTest2")
    test(f2, u2, 0, 0, 1)
    print("\nTest3")
    test(f3, u3, -1, 0, 2)

main()