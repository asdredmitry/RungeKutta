#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
const double omega = 2;
const double EPS = 0.00001;
typedef struct
{
    double y1;
    double y2;
}pair;
double max(double a, double b)
{
    return (a > b) ? a : b;
}
double min(double a, double b)
{
    return (a > b) ? b : a;
}
void fillArrays(double * c, double ** a, double * b)
{
    c[0] = 0;
    c[1] = 0.5;
    c[2] = 2.0/3.0;
    c[3] = 1.0/3.0;
    c[4] = 5.0/6.0;
    c[5] = 1.0/6.0;
    c[6] = 1.0;
    b[0] = 13.0/200.0;
    b[1] = 0.0;
    b[2] = 11.0/40.0;
    b[3] = 11.0/40.0;
    b[4] = 4.0/25.0;
    b[5] = 4.0/25.0;
    b[6] = 13.0/200.0;
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
double * reserve(double * data, int n)
{
    double * tmp = (double *)malloc(sizeof(double)*2*n);
    memcpy(tmp, data, n*sizeof(double));
    free(data);
    return tmp;
}
double f1(double t, double y1, double y2)
{
    return y2;
}
double f2(double t, double y1, double y2) 
{
    return -sin(t);
}
double checkSol(double t)
{
    return sin(t);
}
pair rungeKutta(double h, double y1, double y2, double t, double * c, double ** a, double * b)
{
    double k11, k21, k31, k41, k51, k61, k71;
    double k12, k22, k32, k42, k52, k62, k72; 
    pair output;
    double yn1, yn2;
    k11 = f1(t, y1, y2);
    k12 = f2(t, y1, y2);
    k21 = f1(t + c[1]*h, y1 + a[0][0]*h*k11, y2 + a[0][0]*h*k12);
    k22 = f2(t + c[1]*h, y1 + a[0][0]*h*k11, y2 + a[0][0]*h*k12);
    k31 = f1(t + c[2]*h, y1 + a[1][0]*h*k11 + a[1][1]*h*k21, y2 + a[1][0]*h*k12 + a[1][1]*h*k22);
    k32 = f2(t + c[2]*h, y1 + a[1][0]*h*k11 + a[1][1]*h*k21, y2 + a[1][0]*h*k12 + a[1][1]*h*k22);
    k41 = f1(t + c[3]*h, y1 + a[2][0]*h*k11 + a[2][1]*h*k21 + a[2][2]*h*k31, y2 + a[2][0]*h*k12 + a[2][1]*h*k22 + a[2][2]*h*k32);
    k42 = f2(t + c[3]*h, y1 + a[2][0]*h*k11 + a[2][1]*h*k21 + a[2][2]*h*k31, y2 + a[2][0]*h*k12 + a[2][1]*h*k22 + a[2][2]*h*k32);
    k51 = f1(t + c[4]*h, y1 + a[3][0]*h*k11 + a[3][1]*h*k21 + a[3][2]*h*k31 + a[3][3]*h*k41, y2 + a[3][0]*h*k12 + a[3][1]*h*k22 + a[3][2]*h*k32 + a[3][3]*h*k42);
    k52 = f2(t + c[4]*h, y1 + a[3][0]*h*k11 + a[3][1]*h*k21 + a[3][2]*h*k31 + a[3][3]*h*k41, y2 + a[3][0]*h*k12 + a[3][1]*h*k22 + a[3][2]*h*k32 + a[3][3]*h*k42);
    k61 = f1(t + c[5]*h, y1 + a[4][0]*h*k11 + a[4][1]*h*k21 + a[4][2]*h*k31 + a[4][3]*h*k41 + a[4][4]*h*k51,  y2 + a[4][0]*h*k12 + a[4][1]*h*k22 + a[4][2]*h*k32 + a[4][3]*h*k42 + a[4][4]*h*k52);
    k62 = f2(t + c[5]*h, y1 + a[4][0]*h*k11 + a[4][1]*h*k21 + a[4][2]*h*k31 + a[4][3]*h*k41 + a[4][4]*h*k51,  y2 + a[4][0]*h*k12 + a[4][1]*h*k22 + a[4][2]*h*k32 + a[4][3]*h*k42 + a[4][4]*h*k52);
    k71 = f1(t + c[6]*h, y1 + a[5][0]*h*k11 + a[5][1]*h*k21 + a[5][2]*h*k31 + a[5][3]*h*k41 + a[5][4]*h*k51 + a[5][5]*h*k61, y2 + a[5][0]*h*k12+ a[5][1]*h*k22 + a[5][2]*h*k32 + a[5][3]*h*k42 + a[5][4]*h*k52 + a[5][5]*h*k62);
    k72 = f2(t + c[6]*h, y1 + a[5][0]*h*k11 + a[5][1]*h*k21 + a[5][2]*h*k31 + a[5][3]*h*k41 + a[5][4]*h*k51 + a[5][5]*h*k61, y2 + a[5][0]*h*k12+ a[5][1]*h*k22 + a[5][2]*h*k32 + a[5][3]*h*k42 + a[5][4]*h*k52 + a[5][5]*h*k62);
    yn1 = y1 + h*(b[0]*k11 + b[1]*k21 + b[2]*k31 + b[3]*k41 + b[4]*k51 + b[5]*k61 + b[6]*k71);
    yn2 = y2 + h*(b[0]*k12 + b[1]*k22 + b[2]*k32 + b[3]*k42 + b[4]*k52 + b[5]*k62 + b[6]*k72);
    output.y1 = yn1;
    output.y2 = yn2;
    return output;
}   
double norm(pair tmp1, pair tmp2)
{
    return max(fabs(tmp1.y1 - tmp2.y1), fabs(tmp1.y2 - tmp2.y2));
    //return sqrt(pow(tmp1.y1 - tmp2.y1,2) + pow(tmp2.y2 - tmp1.y2, 2));
}
void solveFixed(double n, double l, double r, double y1, double y2, double * c, double ** a, double * b)
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
        tmp = rungeKutta(h, tmp.y1, tmp.y2, l, c, a, b);
        l += h;
        fprintf(out, "%lf %lf %lf \n", l, tmp.y1, tmp.y2);
    }
    fclose(out);
}
void solveChange(double l, double r, double y1, double y2, double eps, double eps1, double * c, double ** a, double * b)
{
    double h;
    double n;
    pair tmp, tmp1, tmp2, tmp3;
    FILE * out = fopen("data.dat", "w");
    tmp.y1 = y1;
    tmp.y2 = y2;
    n = 0;
    fprintf(out, "%lf %lf %lf \n", l, tmp.y1, tmp.y2);
    while(fabs(l - r) > EPS)
    {
        h = 1;
        hell:
        //printf("%lf \n", h);
        tmp1 = rungeKutta(h, tmp.y1, tmp.y2, l, c, a, b);
        tmp2 = rungeKutta(h/2.0, tmp.y1, tmp.y2, l, c, a, b);
        tmp3 = rungeKutta(h/2.0, tmp2.y1, tmp2.y2, l + h/2.0, c, a, b);
        //printf("%lf %lf %lf %lf \n", tmp1.y1, tmp1.y2, tmp2.y1, tmp2.y2);
        //printf("%lf \n", norm(tmp1, tmp2));
        if(norm(tmp1, tmp3)/(pow(2, 6) - 1) > eps)
        {
            h /= 2.0;
            goto hell;
        }
        else if(norm(tmp1, tmp3)/(pow(2, 6) - 1) < eps1)
        {
            h *= 2;
            goto hell;
        }
        else 
        {  
            tmp = tmp3;
            l += h;
            n = max(fabs(tmp.y1 - checkSol(l)), n);
            fprintf(out, "%lf %lf %lf \n", l, tmp.y1, tmp.y2);
        }
    }
    printf("norm %.17g", n);
    fclose(out);
}
int main()
{
    pair tmp1, tmp2, tmp;
    double * b;
    double * c;
    double ** a;
    int n,i;
    double eps;
    double h;
    double l, r, y1_0, y2_0;
    double t, y, yn;
    double s;
    printf("Input n\n");
    scanf("%d", &n);
    printf("Input l and r\n");
    scanf("%lf %lf", &l, &r);
    printf("Input y_0\n");
    scanf("%lf %lf", &y1_0, &y2_0);
    h = (r - l)/n;
    b = (double *)malloc(sizeof(double)*7);
    c = (double *)malloc(sizeof(double)*7);
    a = (double **)malloc(sizeof(double *)*6);
    if(!a || !b || !c)
    {
        printf("Cannot allocate memory \n");
        return 1;
    }
    for(i = 0; i < 6; i++)
        a[i] = (double *)malloc(sizeof(double)*(i + 1));
    fillArrays(c, a, b);
    //solveFixed(n, l, r, y1_0, y2_0, c, a, b);
    solveChange(l, r, y1_0, y2_0, 0.000001, 0.0000001, c, a, b);
    free(c);
    free(b);
    for(i = 0; i < 6; i++)
        free(a[i]);
    free(a);
    return 0;
}