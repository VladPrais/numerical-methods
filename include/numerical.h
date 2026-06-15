#ifndef NUMERICAL_H
#define NUMERICAL_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

namespace odeint
{

template <typename FuncType, typename InitCond, typename Coefs, size_t N>
void runge_kutta_core(FuncType f, double x, const InitCond& y, InitCond& y_next, double h, const std::array<double, N>& alpha, const std::array<Coefs, N>& beta, const std::array<double, N>& gamma, std::array<InitCond, N>& k)
{
	size_t y_size = y.size();
	InitCond y_temp;
	double temp;
	y_next = y;

	for(size_t i = 0; i < N; i++)
	{
		y_temp = y;
		for(size_t j = 0; j < i; j++)
		{
			temp = h * beta[i][j];
			for(size_t c = 0; c < y_size; c++)
			{
				y_temp[c] += temp * k[j][c];
			}
		}

		k[i] = f(x + h * alpha[i], y_temp);

		temp = h * gamma[i];
		for(size_t j = 0; j < y_size; j++)
		{
			y_next[j] += temp * k[i][j];
		}
	}
}

/* FuncType must be `Range(double, Range)` structure` for all methods. */

template <typename FuncType, typename Domain, typename InitCond>
std::vector<InitCond> euler_method(FuncType f, const Domain &x, const InitCond &y0)
{
	size_t x_size = x.size(), y_size = y0.size();
	double x0, x1, h;
	std::vector<InitCond> sol(x_size, y0);
	sol[0] = y0;
	InitCond dy;

	for(size_t i = 0; i < x_size - 1; i++)
	{
		x0 = x[i];
		x1 = x[i + 1];
		h = x1 - x0;
		dy = f(x0, sol[i]);
		
		for(size_t j = 0; j < y_size; j++)
		{
			sol[i + 1][j] = sol[i][j] + h * dy[j];
		}
	}

	return sol;
}

template <typename FuncType, typename InitCond>
std::vector<InitCond> rk1_classic(FuncType f, double a, double b, size_t n, const InitCond& y0)
{
	double x, h = (b - a) / static_cast<double>(n - 1);

	std::array<double, 1> alpha {0.0};
	std::array<std::array<double, 1>, 1> beta {0.0};
	std::array<double, 1> gamma {1.0};
	std::array<InitCond, 1> k;

	std::vector<InitCond> sol(n);
	sol[0] = y0;

	for(size_t i = 0; i < n - 1; i++)
	{
		x = a + i * h;

		runge_kutta_core(f, x, sol[i], sol[i + 1], h, alpha, beta, gamma, k);
	}
	return std::move(sol);
}

template <typename FuncType, typename Domain, typename InitCond>
std::vector<InitCond> rk2_method(FuncType f, const Domain &x, const InitCond &y0)
{
	size_t x_size = x.size(), y_size = y0.size();
	double x0, x1, h;
	std::vector<InitCond> sol(x_size, y0);
	InitCond temp, k1, k2;
	sol[0] = y0;

	for(size_t i = 0; i < x_size - 1; i++)
	{
		x0 = x[i];
		x1 = x[i + 1];
		h = x1 - x0;
		k1 = f(x0, sol[i]);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * k1[j];
		}

		k2 = f(x0 + h, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			sol[i + 1][j] = sol[i][j] + h / 2.0 * (k1[j] + k2[j]);
		}
	}

	return sol;
}

template <typename FuncType, typename InitCond>
std::vector<InitCond> rk2_classic(FuncType f, double a, double b, size_t n, const InitCond& y0)
{
	double x, h = (b - a) / static_cast<double>(n - 1);

	std::array<double, 2> alpha {0.0, 1.0};
	std::array<std::array<double, 2>, 2> beta {0.0, 0.0, 1.0, 0.0};
	std::array<double, 2> gamma {1.0 / 2.0, 1.0 / 2.0};
	std::array<InitCond, 2> k;

	std::vector<InitCond> sol(n);
	sol[0] = y0;

	for(size_t i = 0; i < n - 1; i++)
	{
		x = a + i * h;

		runge_kutta_core(f, x, sol[i], sol[i + 1], h, alpha, beta, gamma, k);
	}
	return std::move(sol);
}

template <typename FuncType, typename Domain, typename InitCond>
std::vector<InitCond> rk4_method(FuncType f, const Domain &x, const InitCond &y0)
{
	size_t x_size = x.size(), y_size = y0.size();
	double x0, x1, h;
	std::vector<InitCond> sol(x_size, y0);
	InitCond temp, k1, k2, k3, k4;
	sol[0] = y0;

	for(size_t i = 0; i < x_size - 1; i++)
	{
		x0 = x[i];
		x1 = x[i + 1];
		h = x1 - x0;
		k1 = f(x0, sol[i]);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h / 2.0 * k1[j];
		}
		k2 = f(x0 + h / 2.0, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h / 2.0 * k2[j];
		}
		k3 = f(x0 + h / 2.0, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * k3[j];
		}
		k4 = f(x0 + h, temp);

		for(int j = 0; j < y_size; j++)
		{
			sol[i + 1][j] = sol[i][j] + h / 6.0 * (k1[j] + 2.0 * k2[j] + 2.0 * k3[j] + k4[j]);
		}
	}

	return sol;
}

template <typename FuncType, typename InitCond>
std::vector<InitCond> rk4_classic(FuncType f, double a, double b, size_t n, const InitCond& y0)
{
	double x, h = (b - a) / static_cast<double>(n - 1);

	std::array<double, 4> alpha {0.0, 1.0 / 2.0, 1.0 / 2.0, 1.0};
	std::array<std::array<double, 3>, 4> beta {{
		{0.0, 0.0, 0.0},
		{1.0 / 2.0, 0.0, 0.0},
		{0.0, 1.0 / 2.0, 0.0},
		{0.0, 0.0, 1.0}
	}};
	std::array<double, 4> gamma {1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0};
	std::array<InitCond, 4> k;

	std::vector<InitCond> sol(n);
	sol[0] = y0;

	for(size_t i = 0; i < n - 1; i++)
	{
		x = a + i * h;

		runge_kutta_core(f, x, sol[i], sol[i + 1], h, alpha, beta, gamma, k);
	}
	return std::move(sol);
}

template <typename FuncType, typename InitCond>
std::pair<std::vector<double>, std::vector<InitCond>> rk4_adaptive(FuncType f, double a, double b, const InitCond& y0, double atol=1e-7, double rtol=1e-5)
{
	size_t y_size = y0.size();
	double x = a, h = 0.1, max_tol, max_diff, norm_diff, temp;

	std::array<double, 7> alpha {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0};
	const std::array<std::array<double, 6>, 7> beta {{
		{0.0, 0.0, 0.0, 0.0, 0.0},
		{1.0 / 5.0, 0.0, 0.0, 0.0, 0.0},
		{3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0},
		{44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0},
		{19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0},
		{9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0},
		{35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0}
	}};
	std::array<double, 7> gamma {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0};
	std::array<double, 7> diff {-71.0 / 57600.0 , 0.0, 71.0 / 16695.0, -71.0 / 1920.0, 17253.0 / 339200.0, -22.0 / 525.0, 1.0 / 40.0};
	std::array<InitCond, 7> k;

	std::vector<double> x_list;
	std::vector<InitCond> sol;
	InitCond y_next;

	x_list.push_back(a);
	sol.push_back(y0);

	while(std::abs(b - x) > 1e-2)
	{
		if(x + h >= b)
		{
			h = b - x;
		}

		runge_kutta_core(f, x, sol.back(), y_next, h, alpha, beta, gamma, k);

		max_diff = 0.0;
		norm_diff = 0.0;
		for(size_t c = 0; c < y_size; c++)
		{
			max_diff = 0.0;
			for(size_t j = 0; j < k.size(); j++)
			{
				max_diff += diff[j] * k[j][c];
			}
			max_diff = std::max(std::abs(max_diff * h), 1e-9);

			max_tol = atol + rtol * std::max(std::abs(sol.back()[c]), std::abs(y_next[c]));
			temp = max_diff / max_tol;
			norm_diff = norm_diff > temp ? norm_diff : temp;
		}

		if(norm_diff < 1.0)
		{
			x += h;
			x_list.push_back(x);
			sol.push_back(y_next);
		}
		h = 0.9 * h * std::pow(1.0 / norm_diff, 0.2);
	}

	auto result = std::make_pair(x_list, sol);

	return result;
}

/*
template <typename FuncType, typename Domain, typename InitCond>
std::vector<InitCond> rk45_method(FuncType f, const Domain &x, const InitCond &y0, double atol=1e-6)
{
	size_t y_size = y0.size();
	double x0 = x[0], x1 = x[1], h = 0.1;
	std::vector<InitCond> sol;
	sol.push_back(y0);
	InitCond temp, k1, k2, k3, k4, k5, k6, k7;

	const std::array<double, 6> alpha {0.0, 1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0};
	const std::array<std::array<double, 6>, 6> beta {{
		{0.0, 0.0, 0.0, 0.0, 0.0},
		{1.0 / 5.0, 0.0, 0.0, 0.0, 0.0},
		{3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0},
		{44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0},
		{19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0, 0.0},
		{9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0}
	}};
	const std::array<double, 6> y4c {35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0};
	const std::array<double, 7> y5c {5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0};
//	const std::array<double, 7> dc {-71.0 / 57600.0, 0.0, 71.0 / 16695.0, -71.0 / 1920.0, 17253.0 / 339200.0, -22.0 / 525.0, 1.0 / 40.0};

	size_t i = 0;

	while(x0 + h <= x1)
	{
		k1 = f(x0, sol[i]);
		k2 = f(x0 + alpha[0] * h, sol[i] + h * beta[1][0] * k1);
		k3 = f(x0 + alpha[1] * h, sol[i] + h * beta[2][0] * k1 + beta[2][1] * k2);
		k4 = f(x0 + alpha[2] * h, sol[i] + h * beta[3][0] * k1 + beta[3][1] * k2 + beta[3][2] * k3);
		k5 = f(x0 + alpha[3] * h, sol[i] + h * beta[4][0] * k1 + beta[4][1] * k2 + beta[4][2] * k3 + beta[4][3] * k4);
		k6 = f(x0 + alpha[4] * h, sol[i] + h * beta[5][0] * k1 + beta[5][1] * k2 + beta[5][2] * k3 + beta[5][3] * k4 + beta[5][4] * k5);

		temp = sol[i] + h * (y4c[0] * k1 + y4c[1] * k2 + y4c[2] * k3 + y4c[3] * k4 + y4c[4] * k5 + y4c[5] * k6);

		k7 = f(x0 + alpha[5] * h, temp);

		sol[i + 1]= sol[i] + h * (y5c[0] * k1 + y5c[1] * k2 + y5c[2] * k3 + y5c[3] * k4 + y5c[4] * k5 + y5c[5] * k6 + y5c[6] * k7);
	}

	while(x0 + h <= x1)
	{
		k1 = f(x0, sol[i]);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * beta[1][0] * k1[j];
		}
		k2 = f(x0 + h / 5.0, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * (beta[2][0] * k1[j] + beta[2][1] * k2[j]);
		}
		k3 = f(x0 + h * 3.0 / 10.0, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * (beta[3][0] * k1[j] + beta[3][1] * k2[j] + beta[3][2] * k3[j]);
		}
		k4 = f(x0 + h * 4.0 / 5.0, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * (beta[4][0] * k1[j] + beta[4][1] * k2[j] + beta[4][2] * k3[j] + beta[4][3] * k4[j]);
		}
		k5 = f(x0 + h * 8.0 / 9.0, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * (beta[5][0] * k1[j] + beta[5][1] * k2[j] + beta[5][2] * k3[j] + beta[5][3] * k4[j] + beta[5][4] * k5[j]);
		}
		k6 = f(x0 + h, temp);

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * (y4c[0] * k1[j] + y4c[1] * k2[j] + y4c[2] * k3[j] + y4c[3] * k4[j] + y4c[4] * k5[j] + y4c[5] * k6[j]);
		}
		k7 = f(x0 + h, temp);
		for(size_t j = 0; j < y_size; j++)
		{
			k7[j] = sol[i][j] + h * (y5c[0] * k1[j] + y5c[1] * k2[j] + y5c[2] * k3[j] + y5c[3] * k4[j] + y5c[4] * k5[j] + y5c[5] * k6[j] + y5c[6] * k7[j]);
		}
		for(size_t j = 0; j < y_size; j++)
		{

		}

		for(size_t j = 0; j < y_size; j++)
		{
			temp[j] = std::abs(dc[0] * k1[j] + dc[1] * k2[j] + dc[2] * k3[j] + dc[3] * k4[j] + dc[4] * k5[j] + dc[5] * k6[j] + dc[6] * k7[j]);
		}

		double abs_v = 0.0;
		for(auto v: temp)
		{
			abs_v += v * v;
		}
		abs_v = std::sqrt(abs_v);
		h = 0.9 * h * std::pow((atol / abs_v), 1.0 / 5.0);

		if(abs_v < atol)
		{
			x0 += h;
			i++;
		}

	}
	return sol;
}
*/

} // end of namespace

template <typename FuncType>
double simple_iteration_method(FuncType f, double x0, double eps=1e-6, int n=1e3)
{
	double x1;
	for(int i = 0; i < n; i++)
	{
		x1 = f(x0);
		if(std::abs(x1 - x0) < eps)
		{
			break;
		}
		x0 = x1;
	}
	return x1;
}

template <typename FuncType>
double newton_method(FuncType f, double x0, double eps=1e-6, int n=1e3)
{
	double x1, df, h = 1e-5;
	for(int i = 0; i < n; i++)
	{
		df = (f(x0 + h) - f(x0)) / h;
		x1 = x0 - f(x0) / df;
		if(std::abs(x1 - x0) < eps)
		{
			break;
		}
		x0 = x1;
	}
	return x1;
}

template <typename FuncType>
double secant_method(FuncType f, double x0, double x1, double eps=1e-6, int n=1e3)
{
	double x2;
	for(int i = 0; i < n; i++)
	{
		//assert(f(x1) - f(x0) > 1e-6);
		x2 = (f(x1) * x0 - f(x0) * x1) / (f(x1) - f(x0));

		if(std::abs(x2 - x1) < eps)
		{
			break;
		}
		x0 = x1;
		x1 = x2;
	}
	return x2;
}

template <typename FuncType>
double rectangle_method(FuncType& f, double a, double b, double h)
{
	double r = 0.0;

	for(double k = a; k < b; k += h)
	{
		r += h * f(k);
	}

	return r;
}

#endif
