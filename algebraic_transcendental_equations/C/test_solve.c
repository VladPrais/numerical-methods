#include <stdio.h>
#include <math.h>
#include "solve.h"

double w1(double);
double w2(double);
double w3(double);
double f1(double);
double f2(double);
double f3(double);
void pprint(double, double, double);

int main(void)
{
	double x0, x1;
	double sol1, sol2, sol3;

	x0 = 0.0;
	x1 = 0.1;
	
	sol1 = iteration_method(w1, x0);
	sol2 = newton_method(f1, x0);
	sol3 = false_position_method(f1 ,x0, x1);
	
	pprint(sol1, sol2, sol3);

	x0 = 0.0;
	x1 = 0.1;

	sol1 = iteration_method(w2, x0);
	sol2 = newton_method(f2, x0);
	sol3 = false_position_method(f2, x0, x1);
	
	pprint(sol1, sol2, sol3);

	x0 = 0.5;
	x1 = 0.6;

	sol1 = iteration_method(w3, x0);
	sol2 = newton_method(f3, x0);
	sol3 = false_position_method(f3, x0, x1);
	
	pprint(sol1, sol2, sol3);

	return 0;
}

// pretty print
void pprint(double sol1, double sol2, double sol3)
{
	printf("Iteration method\tNewton method\tFalse position method\
		\n%f\t\t%f\t%f\n\n", sol1, sol2, sol3);
}

// for iteration method
double w1(double x)
{
	return cos(x);
}
double w2(double x)
{
	return exp(-x) * cos(x);
}
double w3(double x)
{
	return sqrt(x) + 1;
}

// for Newton and false position methods
double f1(double x)
{
	return x - cos(x);
}
double f2(double x)
{
	return x - exp(-x) * cos(x);
}
double f3(double x)
{
	return x - sqrt(x) - 1;
}
