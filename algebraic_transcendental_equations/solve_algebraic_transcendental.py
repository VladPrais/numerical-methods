def iteration_method(x0, y, eps = 1e-4):
    k = 0
    
    while True:
        y0 = y(x0)
        k += 1

        if abs(y0 - x0) < eps:
            break
            
        x0 = y0

    print(f"{k} iterations using iteration method")
    return x0

def newton_method(x0, y, eps = 1e-4):
    h = 1e-3
    k = 0
    
    while True:
        dy = (y(x0 + h) - y(x0)) / h
        xi = x0 - y(x0) / dy
        k += 1
        
        if abs(xi - x0) < eps:
            break

        x0 = xi
        
    print(f"{k} iterations using Newton method")
    return xi

def false_position_method(x0, x1, y, eps = 1e-4):
    k = 0
    
    while True:
        x2 = (y(x1) * x0 - y(x0) * x1) / (y(x1) - y(x0))
        k += 1
        
        if abs(x1 - x2) < eps:
            break

        x0, x1 = x1, x2
        
    print(f"{k} iterations using false position method")
    return x2