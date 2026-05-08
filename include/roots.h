#ifndef ROOTS_H
#define ROOTS_H

#include <cmath>

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

#endif
