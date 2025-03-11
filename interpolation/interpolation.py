def lagrange_polynomial(x, y, x0):
    N = len(x)
    S = 0
    for i in range(N):
        mul = 1
        for j in range(N):
            if i != j:
                mul *= (x0 - x[j]) / (x[i] - x[j])
        S += y[i] * mul
    return S