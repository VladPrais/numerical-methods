#ifndef ROOTS_H
#define ROOTS_H

#include <cmath>

template <typename FuncType>
double simple_iteration_method(FuncType f, double x0, double eps=1e-6)
{
	double x1 = f(x0);
	while(std::abs(x1 - x0) > eps)
	{
		x0 = x1;
		x1 = f(x0);
	}
	return x1;
}

template <typename FuncType>
double newton_method(FuncType f, double x0, double eps=1e-6)
{
	double h = 1e-5, df = (f(x0 + h) - f(x0)) / h;
	double x1 = x0 - f(x0) / df;
	while(std::abs(x1 - x0) > eps)
	{
		x0 = x1;
		df = (f(x0 + h) - f(x0)) / h;
		x1 = x0 - f(x0) / df;
	}
	return x1;
}

template <typename FuncType>
double secant_method(FuncType f, double x0, double x1, double eps=1e-6)
{
	double x2 = (f(x1) * x0 - f(x0) * x1) / (f(x1) - f(x0));
	while(std::abs(x2 - x1) > eps)
	{
		x0 = x1;
		x1 = x2;
		x2 = (f(x1) * x0 - f(x0) * x1) / (f(x1) - f(x0));
	}
	return x2;
}

#endif
