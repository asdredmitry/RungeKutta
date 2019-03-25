#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
const double lol = 0;
const double omega = 2;
double alpha = 0.2;
const double EPS = 0.00001;
const double minh = 0.00001;
const double c[] = {0, 0.5, 2.0/3.0, 1.0/3.0, 5.0/6.0, 1.0/6.0, 1.0};
const double b[] = {13.0/200.0, 0.0, 11.0/40.0, 11.0/40.0, 4.0/25.0, 4.0/25.0, 13.0/200.0};
double ** a;
typedef struct
{
    double y1;
    double y2;
}pair;
typedef struct 
{
    double * data;
    int i; 
    int n;
}vector;
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
void fillArrays()
{
    a[0][0] = 1.0/2.0;
    a[1][0] = 2.0/9.0;
    a[1][1] = 4.0/9.0;
    a[2][0] = 7.0/36.0;
    a[2][1] = 2.0/9.0;
    a[2][2] = -1.0/12.0;
    a[3][0] = -35.0/144.0;
    a[3][1] = -55.0/36.0;
    a[3][2] = 35.0/48.0;
    a[3][3] = 15.0/8.0;
    a[4][0] = -1.0/360.0;
    a[4][1] = -11.0/36.0;
    a[4][2] = -1.0/8.0;
    a[4][3] = 1.0/2.0;
    a[4][4] = 1.0/10.0;
    a[5][0] = -41.0/260.0;
    a[5][1] = 22.0/13.0;
    a[5][2] = 43.0/156.0;
    a[5][3] = -118.0/39.0;
    a[5][4] = 32.0/195.0;
    a[5][5] = 80.0/39.0;
}
double f1(double t, double y1, double y2)
{
    return y2;
}
double f2(double t, double y1, double y2) 
{
    return -y1;
    //return -y1*y1*y1 + alpha*cos(t);
    //return (1. + 3.*cos(2.*t))/2. - y2,
    //return -(1 + alpha*y1*y1)*y1 + cos(t);
    //return -sin(t);
}
double checkSol(double t)
{
    return sin(t);
}
pair rungeKutta(double h, double y1, double y2, double t)
{ 
    double tmpt, tmpy1, tmpy2;
    double k1[7];
    double k2[7];
    pair output;
    int i, j;
    k1[0] = f1(t, y1, y2);
    k2[0] = f2(t, y1, y2);
    for(i = 1; i < 7; i++)
    {
        tmpy1 = y1;
        tmpy2 = y2;
        tmpt = t + c[i + 1]*h;
        for(j = 0; j < i + 1; j++)
        {
            tmpy1 += a[i - 1][j]*h*k1[i];
            tmpy2 += a[i - 1][j]*h*k2[i];
        }
        k1[i] = f1(tmpt, tmpy1, tmpy2);
        k2[i] = f2(tmpt, tmpy1, tmpy2);
    }
    output.y1 = 0.;
    output.y2 = 0.;
    for(i = 0;i < 7; i++)
    {
        output.y1 += b[i]*k1[i];
        output.y2 += b[i]*k2[i];
    }
    output.y1 *= h;
    output.y2 *= h;
    output.y1 += y1;
    output.y2 += y2;
    return output;
}   
double norm(pair tmp1, pair tmp2)
{
    return max(fabs(tmp1.y1 - tmp2.y1), fabs(tmp1.y2 - tmp2.y2));
    //return sqrt(pow(tmp1.y1 - tmp2.y1,2) + pow(tmp2.y2 - tmp1.y2, 2));
}
pair findSolution(double t0, double y10, double y20, double t, double tol)
{
    double h, err; 
    pair tmp, tmp1, tmp2, tmp3;
    double fac, facmax, facmin;
    fac = 0.8;
    facmax = 1.5;
    facmin = 0;
    h = 0.001;
    tmp.y1 = y10;
    tmp.y2 = y20;
    while(t0 < t)
    {
        tmp1 = rungeKutta(h, tmp.y1, tmp.y2, t0);
        tmp2 = rungeKutta(h, tmp1.y1, tmp1.y2, t0 + h);
        tmp3 = rungeKutta(2*h, tmp.y1, tmp.y2, t0);
        err = norm(tmp2, tmp3)/(pow(2, 7) - 1);
        if(err < tol)
        {
            if(h > (t - t0))
                return rungeKutta(h, tmp.y1, tmp.y2, t0);
            else if(2*h > (t - t0))
            {}
        }
    }
}
double getH(double h, double err, double tol)
{
    double fac, facmax, facmin;
    fac = 0.8;
    facmax = 1.5;
    facmin = 0;
    h = h*min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));    
    return h;
}
pair findVal(double t0, double y10, double y20, double t, double tol)
{
    pair tmp, tmp1, tmp2, tmp3;
    double h;
    double err, fac, facmax, facmin;
}
void solveFixed(double n, double l, double r, double y1, double y2)
{
    double h;
    pair tmp;
    FILE * out = fopen("check.dat", "w");
    h = (r - l)/n;
    tmp.y1 = y1;
    tmp.y2 = y2;
    fprintf(out, "%lf %lf %lf \n",l ,tmp.y1, tmp.y2);
    while(l <= r)
    {
        tmp = rungeKutta(h, tmp.y1, tmp.y2, l);
        l += h;
        fprintf(out, "%lf %lf %lf \n", l, tmp.y1, tmp.y2);
    }
    fclose(out);
}
void solveChangeL(vector * t, vector * s1, vector * s2, double l, double r, double y1, double y2, double tol)
{
    double h, err;
    pair tmp, tmp1, tmp2, tmp3;
    h = -0.001;
    tmp.y1 = y1;
    tmp.y2 = y2;
    push_back(t, r);
    push_back(s1, y1);
    push_back(s2, y2);
    while(l < r)
    {
        tmp1 = rungeKutta(h, tmp.y1, tmp.y2, r);
        tmp2 = rungeKutta(h, tmp1.y1, tmp1.y2, r + h);
        tmp3 = rungeKutta(2*h, tmp.y1, tmp.y2, r);
        err = norm(tmp2, tmp3)/(pow(2, 7) - 1);
        if(err < tol)
        {
            push_back(t, r + h);
            push_back(t, r + 2*h);
            push_back(s1, tmp1.y1);
            push_back(s1, tmp2.y1);
            push_back(s2, tmp1.y2);
            push_back(s2, tmp2.y2);
            tmp = tmp2;
            r += 2*h;
            h = getH(h, err, tol);
        }
        else
        {
            h = getH(h, err, tol);
        }        
    }
}
void solveChange(vector * t, vector * s1, vector * s2, double l, double r, double y1, double y2, double tol)
{
   double h, err; 
   pair tmp, tmp1, tmp2, tmp3;
   h = 0.001;
   tmp.y1 = y1;
   tmp.y2 = y2;
   push_back(t, l);
   push_back(s1, y1);
   push_back(s2, y2);
   while(l < r)
   {
       tmp1 = rungeKutta(h, tmp.y1, tmp.y2, l);
       tmp2 = rungeKutta(h, tmp1.y1, tmp1.y2, l + h);
       tmp3 = rungeKutta(2*h, tmp.y1, tmp.y2, l);
       printf("%17g %lf \n", h, l);
       err = norm(tmp2, tmp3)/(pow(2, 7) - 1);
       if(err < tol)
       {
           if(2*h > (r - l))
           {
               h = (r - l)/2;
               tmp1 = rungeKutta(h, tmp.y1, tmp.y2, l);
               tmp2 = rungeKutta(h, tmp1.y1, tmp1.y2, l + h);
               tmp3 = rungeKutta(2*h, tmp.y1, tmp.y2, l);
           }
           push_back(t, l + h);
           push_back(s1, tmp1.y1);
           push_back(s2, tmp1.y2);
           push_back(t, l + 2*h);
           push_back(s1, tmp2.y1);
           push_back(s2, tmp2.y2);
           tmp = tmp2;
           l += 2*h;
           h = getH(h, err, tol);
           //h = h*min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));
       }
       else
       {
           h = getH(h, err, tol);
           //h = h*min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));
       }       
   }
}
void write_data(vector * tr, vector * s1r, vector * s2r, vector * tl, vector * s1l, vector * s2l)
{
    int i;
    FILE * output = fopen("data.dat", "w");
    if(!output)
    {
        printf("Cannot open file \n");
        exit(1);
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
int main()
{
    int n,i;
    double l, r, y1_0, y2_0;
    double t_0;
    vector tr, s1r, s2r, tl, s1l, s2l;
    init(&s1r, 10);
    init(&tr, 10);
    init(&s2r, 10);
    init(&s1l, 10);
    init(&s2l, 10);
    init(&tl, 10);
    printf("Input n\n");
    scanf("%d", &n);
    printf("Input l and r\n");
    scanf("%lf %lf", &l, &r);
    printf("Input t_0 and  y_0\n");
    scanf("%lf %lf %lf", &t_0, &y1_0, &y2_0);
    a = (double **)malloc(sizeof(double *)*6);
    if(!a)
    {
        printf("Cannot allocate memory \n");
        return 1;
    }
    for(i = 0; i < 6; i++)
    {
        a[i] = (double *)malloc(sizeof(double)*(i + 1));
        if(!a[i])
        {
            printf("Cannot allocate memory \n");
            return 1;
        }
    }
    fillArrays();
    solveChange(&tr, &s1r, &s2r, t_0, r, y1_0, y2_0, 0.000000001);
    solveChangeL(&tl, &s1l, &s2l, l, t_0, y1_0, y2_0, 0.000000001);
    //solveFixed(n, l, r, y1_0, y2_0, c, a, b);
    //solveChange(&t, &s1, &s2, l, r, y1_0, y2_0, 0.000000001, c, a, b);
    write_data(&tr, &s1r, &s2r, &tl, &s1l, &s2l);
    for(i = 0; i < 6; i++)
        free(a[i]);
    free(a);
    free_vec(&tl);
    free_vec(&s1l);
    free_vec(&s2l);
    free_vec(&tr);
    free_vec(&s1r);
    free_vec(&s2r);
    return 0;
}