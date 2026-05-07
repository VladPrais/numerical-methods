#ifndef DSOLVE_H
#define DSOLVE_H

#include <vector>

template <typename FuncType, typename Domain, typename Range>
std::vector<Range> euler_method(FuncType f, const Domain &x, const Range &y0)
{
	double x_size = x.size(), y_size = y0.size();
	double x0, x1, h;
	std::vector<Range> sol(x_size, y0);
	Range dy;
	sol[0] = y0;

	for(int i = 0; i < x_size - 1; i++)
	{
		x0 = x[i];
		x1 = x[i + 1];
		h = x1 - x0;
		dy = f(x0, sol[i]);
		
		for(int j = 0; j < y_size; j++)
		{
			sol[i + 1][j] = sol[i][j] + h * dy[j];
		}
	}

	return sol;
}

template <typename FuncType, typename Domain, typename Range>
std::vector<Range> rk2_method(FuncType f, const Domain &x, const Range &y0)
{
	double x_size = x.size(), y_size = y0.size();
	double x0, x1, h;
	std::vector<Range> sol(x_size, y0);
	Range temp, k1, k2;
	sol[0] = y0;

	for(int i = 0; i < x_size - 1; i++)
	{
		x0 = x[i];
		x1 = x[i + 1];
		h = x1 - x0;
		k1 = f(x0, sol[i]);

		for(int j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h * k1[j];
		}

		k2 = f(x0 + h, temp);

		for(int j = 0; j < y_size; j++)
		{
			sol[i + 1][j] = sol[i][j] + h / 2.0 * (k1[j] + k2[j]);
		}
	}

	return sol;
}

template <typename FuncType, typename Domain, typename Range>
std::vector<Range> rk4_method(FuncType f, const Domain &x, const Range &y0)
{
	double x_size = x.size(), y_size = y0.size();
	double x0, x1, h;
	std::vector<Range> sol(x_size, y0);
	Range temp, k1, k2, k3, k4;
	sol[0] = y0;

	for(int i = 0; i < x_size - 1; i++)
	{
		x0 = x[i];
		x1 = x[i + 1];
		h = x1 - x0;
		k1 = f(x0, sol[i]);

		for(int j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h / 2.0 * k1[j];
		}
		k2 = f(x0 + h / 2.0, temp);

		for(int j = 0; j < y_size; j++)
		{
			temp[j] = sol[i][j] + h / 2.0 * k2[j];
		}
		k3 = f(x0 + h / 2.0, temp);

		for(int j = 0; j < y_size; j++)
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

#endif
