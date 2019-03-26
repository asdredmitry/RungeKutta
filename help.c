#include "help.h"
double norm(pair tmp1, pair tmp2)
{
    return max(fabs(tmp1.y1 - tmp2.y1), fabs(tmp1.y2 - tmp2.y2));
}
void init(vector * v, int n)
{
    v->data = (double *)malloc(n*sizeof(double));
    v->i = 0;
    v->n = n;
}
void reall(vector * v)
{
    double * tmp;
    tmp = (double *)malloc(sizeof(double)*(2*v->n));
    memcpy(tmp, v->data, (v->n)*sizeof(double));
    v->n *= 2;
    free(v->data);
    v->data = tmp;
}
void push_back(vector * v, double val)
{
    if(v->i == v->n)
        reall(v);
    v->data[v->i] = val;
    v->i++;
}
void free_vec(vector * v)
{
    if(v->data != NULL)
        free(v->data);
    v->data = NULL;
    v->i = 0;
    v->n = 0;
}

double max(double x, double y)
{
    return (x > y) ? x : y;
}
double min(double x, double y)
{
    return (x > y) ? y : x;
}
void write_data(vector * tr, vector * s1r, vector * s2r, vector * tl, vector * s1l, vector * s2l)
{
    int i;
    FILE * output = fopen("data.dat", "w");
    if(!output)
    {
        printf("Cannot open file \n");
        exit(EXIT_FAILURE);
    }
    for(i = tl->i - 1; i >= 0; i--)
    {
        fprintf(output, "%lf %lf %lf \n", tl->data[i], s1l->data[i], s2l->data[i]);
    }
    for(i = 0; i < tr->i; i++)
    {
        fprintf(output, "%lf %lf %lf \n", tr->data[i], s1r->data[i], s2r->data[i]);
    }
    fclose(output);
}
void write_data_1(vector * t, vector * s1, vector * s2)
{
    int i;
    FILE * output = fopen("data.dat", "w");
    if(!output)
    {
        printf("Cannot open file \n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < t->i; i++)
    {
        fprintf(output, "%lf %lf %lf \n", t->data[i], s1->data[i], s2->data[i]);
    }
    fclose(output);
}