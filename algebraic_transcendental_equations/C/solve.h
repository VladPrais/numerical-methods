#ifndef __SOLVE_H__
#define __SOLVE_H__

double iteration_method(double (*)(double), double);
double newton_method(double (*)(double), double);
double false_position_method(double (*)(double), double, double);

#endif
