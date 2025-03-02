This repository is a little clue, that can help you to undestand how to
- solve **differential equations** using Euler and Runge-Kutta methods.
### Solving initial value problem
Let's write the formulas, assuming that:
$$h=x_{n+1}-x_n$$

##### Euler method:
$$y_{n+1} = y_n + hf(x_n, y_n)$$

##### The second order Runge-Kutta method:
$$y_{n+1} = y_n + hf(x_n + \frac{h}{2}, y_n + \frac{h}{2}f(x_n,y_n))$$

##### The fourth order Runge-Kutta method:
$$y_{n+1} = y_n + \frac{h}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$
Where
$$k_1 = f(x_n,y_n)$$
$$k_2 = f\left( x_n + \frac{h}{2},y+\frac{h}{2}k_{1} \right)$$
$$k_3 = f\left( x_n + \frac{h}{2},y+\frac{h}{2}k_{2} \right)$$
$$k_{4} = f(x_n + h, y_n + hk_{3})$$