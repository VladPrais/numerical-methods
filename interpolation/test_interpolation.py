import numpy as np
import matplotlib.pyplot as plt
import random
from interpolation import *

# ======================================================
# Function for interpolation
# ======================================================

# f(x) = cos(3x)
def f1(x):
    return np.cos(3 * x)

# f(x) = e^(-x)sinx
def f2(x):
    return np.exp(-x) * np.sin(x)

# f(x) = 2x^2 - x + 1
def f3(x):
    return -3 * x**3 + 2 * x**2 - x + 1


# ======================================================
# Test function
# ======================================================

def test(y, a, b, max_nodes):

    # step for split points
    h = 0.1

    for i in range(2, max_nodes+1):
        # create the random nodes
        x0 = np.array([random.uniform(a, b) for _ in range(i)])
        
        # function values in these nodes
        y0 = y(x0)
    
        x = np.arange(a, b + h, h)

        plt.scatter(x0, y0, label = "nodes")
        plt.scatter(x, lagrange_polynomial(x0, y0, x), label = "polynomial", marker = "x")
        plt.plot(x, y(x), label = "func")

        plt.legend()
        plt.title(f"Interpolation using {i} nodes")
        plt.show()

# ======================================================
# Tests
# ======================================================

def main():

    test(f1, 0, 3, 7)
    test(f2, 0, 3, 10)
    test(f3, -1, 1, 4)

main()
