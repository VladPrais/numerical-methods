#include <math.h>

#define EPS		1e-5
#define MAX_STEP	300

double iteration_method(double (*f)(double), double x0)
{
	double x1;
	for(int i = 0; i < MAX_STEP; ++i){
		x1 = f(x0);

		if(fabs(x1 - x0) < EPS){
			break;
		}

		x0 = x1;
	}
	return x1;
}

double newton_method(double (*f)(double), double x0)
{
	double df_dx, h = 1e-3, x1;

	for(int i = 0; i < MAX_STEP; ++i){
		df_dx = (f(x0 + h) - f(x0)) / h;
		x1 = x0 - f(x0) / df_dx;

		if(fabs(x1 - x0) < EPS){
			break;
		}

		x0 = x1;
	}
	return x1;
}

double false_position_method(double (*f)(double), double x0, double x1)
{
	double x2;

	for(int i = 0; i < MAX_STEP; ++i){
		x2 = (f(x1) * x0 - f(x0) * x1) / (f(x1) - f(x0));

		if(fabs(x2 - x1) < EPS){
			break;
		}

		x0 = x1;
		x1 = x2;
	}
	return x1;
}
