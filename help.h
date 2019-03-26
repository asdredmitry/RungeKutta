#ifndef HELP_H
#define HELP_H
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

double * c;
double * b;
double ** a;

typedef struct 
{
    double y1;
    double y2;
}pair;

double norm(pair tmp1, pair tmp2);

typedef struct 
{
    double * data;
    int i;
    int n;
}vector;

void init(vector * v, int n);
void reall(vector * v);
void push_back(vector * v, double val);
void free_vec(vector * v);


double max(double x, double y);
double min(double x, double y);

void write_data(vector * tr, vector * s1r, vector * s2r, vector * tl, vector * s1l, vector * s2l);
void write_data_1(vector * t, vector * s1, vector * s2);

double norm_full(vector * t, vector * s1);
double check_func(double t);

#endif