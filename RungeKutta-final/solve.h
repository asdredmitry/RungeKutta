#ifndef SOLVE_H
#define SOLVE_H
#include "help.h"


void fillArrays(void);

double f1(double t, double y1, double y2);
double f2(double t, double y1, double y2);


pair runge_kutta(double h, double y1, double y2, double t);

double getH(double h, double err, double tol);

pair findSolution(double t0, double y10, double y20, double t, double tol);

void solve_fixed(double n, double l, double r, double y1, double y2, vector * t, vector * s1, vector * s2);

void solve_change(double l, double r, double y1, double y2, double tol, vector * t, vector * s1, vector * s2);

pair find_val(double t0, vector * t, vector * s1, vector * s2);
#endif
