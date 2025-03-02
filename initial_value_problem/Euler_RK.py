import numpy as np

def euler_method(f, x, y0):
    x = np.asarray(x, dtype=float)

    xs = x.size

    y = np.zeros((xs, len(y0)))
    y[0] = y0
    for i in range(xs - 1):
        h = x[i + 1] - x[i]
        y[i + 1] = y[i] + h * f(x[i], y[i])

    return y.T

def runge_kutta_2_method(f, x, y0):
    x = np.asarray(x, dtype=float)
    
    xs = x.size
    
    y = np.zeros((xs, len(y0)))
    y[0] = y0
    
    for i in range(xs - 1):
        h = x[i + 1] - x[i]
        
        y[i + 1] = y[i] + h * f(x[i] + h/2, y[i] + h/2 * f(x[i], y[i]))
	
    return y.T

def runge_kutta_4_method(f, x, y0):
    x = np.asarray(x, dtype=float)
    
    xs = x.size
    
    y = np.zeros((xs, len(y0)))
    y[0] = y0
    
    for i in range(xs - 1):
        h = x[i + 1] - x[i]
        k1 = f(x[i], y[i])
        k2 = f(x[i] + h/2, y[i] + h/2 * k1)
        k3 = f(x[i] + h/2, y[i] + h/2 * k2)
        k4 = f(x[i] + h, y[i] + h * k3)
        y[i + 1] = y[i] + h/6 * (k1 + 2 * k2 + 2 * k3 + k4)
        
    return y.T
