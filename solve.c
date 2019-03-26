#include "solve.h"
const double alpha = 0.2;
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
void fillArrays()
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
double getH(double h, double err, double tol)
{
    double fac, facmax, facmin;
    fac = 0.8;
    facmax = 1.5;
    facmin = 0;
    h *= min(facmax, max(facmin, fac*pow(tol/err, 1./7.)));
    return h;
}
pair runge_kutta(double h, double y1, double y2, double t)
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
        tmpt = t + c[i]*h;
        for(j = 0; j < i; j++)
        {
            tmpy1 += a[i - 1][j]*h*k1[j];
            tmpy2 += a[i - 1][j]*h*k2[j];
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
void solve_fixed(double n, double l, double r, double y1, double y2, vector * t, vector * s1, vector * s2)
{
    double h;
    pair tmp;
    h = (r - l)/n;
    push_back(t, l);
    push_back(s1, y1);
    push_back(s2, y2);
    while(l <= r)
    {
        tmp = runge_kutta(h, s1->data[s1->i - 1], s2->data[s2->i - 1], t->data[t->i - 1]);
        l += h;
        push_back(t, l);
        push_back(s1, tmp.y1);
        push_back(s2, tmp.y2);
    }    
}
void solve_change(double l, double r, double y1, double y2, double tol, vector * t, vector * s1, vector * s2)
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
        tmp1 = runge_kutta(h, tmp.y1, tmp.y2, l);
        tmp2 = runge_kutta(h, tmp1.y1, tmp1.y2, l + h);
        tmp3 = runge_kutta(2*h, tmp.y1, tmp.y2, l);
        err = norm(tmp2, tmp3)/(pow(2., 7) - 1);
        if(err < tol)
        {
            if(l + h > r)
            {
                tmp = runge_kutta(r - l, tmp.y1, tmp.y2, l);
                push_back(t, r);
                push_back(s1, tmp.y1);
                push_back(s2, tmp.y2);
                break;
            }
            push_back(t, l + h);
            push_back(s1, tmp1.y1);
            push_back(s2, tmp1.y2);
            if(l + 2*h > r)
            {
                tmp = runge_kutta(r - (l + h), tmp1.y1, tmp1.y2, l + h);
                push_back(t, r);
                push_back(s1, tmp.y1);
                push_back(s2, tmp.y2);
                break;
            }
            push_back(t, l + 2*h);
            push_back(s1, tmp2.y1);
            push_back(s2, tmp2.y2);
            tmp = tmp2;
            l += 2*h;
            h = getH(h, err, tol);
        }
        else 
        {
            h = getH(h, err, tol);
        }
    }
}