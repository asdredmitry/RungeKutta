#include <stdio.h>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
using namespace std;
double f(double x, double y)
{
    return 1/(1 + pow(x,2));
}
double y_check(double x)
{
    return atan(x);
}
vector<pair<double, double> > rungeKutta4(double l, double r, double h, double n, double y_0)
{
    vector<pair<double, double> > result;
    result.reserve(n);
    result.push_back(make_pair(l, y_0));
    for(double i = l + h; i <= r; i += h)
    {
        double xn = result[result.size() - 1].first;
        double yn = result[result.size() - 1].second;
        double k1 = h*f(xn, yn);
        double k2 = h*f(xn + h/2, yn + k1/2);
        double k3 = h*f(xn + h/2, yn + k2/2);
        double k4 = h*f(xn + h, yn + k3);
        double ynn = yn + k1/6 + k2/3 + k3/3 + k4/6;
        result.push_back(make_pair(i, ynn));
    }
    return result;
}
vector<pair<double, double> > rungeKutta6(double l, double r, double h, int n, double y_0)
{
    vector<pair<double, double> > result;
    result.reserve(n);
    result.push_back(make_pair(l, y_0));
    //for(double i = l + h; i <= r; i += h)
    //{
        /*double xn = result[result.size() - 1].first;
        double yn = result[result.size() - 1].second;
        double k1 = f(xn, yn);
        double k2 = f(xn + h/3, yn + (h*k1)/3);
        double k3 = f(xn + 2*h/3, yn + 2*h*k2/3);
        double k4 = f(xn + h/3, yn + h/12*(k1 + 4.0*k2 - k3));
        double k5 = f(xn + h/2, yn + h/16*(-k1 + 18.0*k2 - 3.0*k3 -6.0*k4));
        double k6 = f(xn + h/2, yn + h/8*(9.0*k2 - 3.0*k3  - 6.0*k4 + 4.0*k5));
        double k7 = f(xn + h, yn + h/44*(9.0 * k1 - 36.0*k2 + 63.0 *k3 + 72.0*k4 - 64.0*k6));
        double ynn = h/120*(11.0*(k1 + k7) + 81.0*(k3 + k4) - 32.0*(k5 + k6));
        */
          double k1, k2, k3, k4, k5, k6, k7;
   double h3 = h/3;
   double h2_3 = 2*h/3;
   double h2 = h/2;
   double h12 = h/12;
   double h16 = h/16;
   double h8 = h/8;
   double h44 = h/44;
   double h120 = h/120;

   for(double i = l + h; i <= r; i += h){
       double x0 = result[result.size() - 1].first;
       double y0 = result[result.size() - 1].second;
        k1 = f(x0,y0);
        k2 = f(x0 + h3, y0 + h3 * k1);
        k3 = f(x0 + h2_3, y0 + h2_3 * k2 );
        k4 = f(x0 + h3, y0 + h12 * ( k1 + 4.0 * k2 - k3 ) );
        k5 = f(x0 + h2, y0 
                           + h16 * ( -k1 + 18.0 * k2 - 3.0 * k3 - 6.0 * k4 ) );
        k6 = f(x0 + h2, y0 + h8 * ( 9.0 * k2 - 3.0 * k3 - 6.0 * k4
                                                                + 4.0 * k5 ) );
        k7 = f(x0 + h, y0 + h44 * ( 9.0 * k1 - 36.0 * k2 + 63.0 * k3
                                                   + 72.0 * k4 - 64.0 * k6) );
        y0 += h120 * ( 11.0 * (k1 + k7) + 81.0 * (k3 + k4) - 32.0 * (k5 + k6) );
        x0 += h;
        result.push_back(make_pair(x0, y0));
    }
    return result;
}
vector<pair<double, double> > rungeKutta6_1(double l, double r, double h, int n, double y_0)
{
    vector<pair<double, double> > result;
    result.reserve(n);
    result.push_back(make_pair(l, y_0));
    double sqrt21 = sqrt(21);
    for(double i = l + h; i <= r; i += h)
    {
        double x = result[result.size() - 1].first;
        double y = result[result.size() - 1].second;
        double k1 = h*f(x, y);
        double k2 = h*f(x + h, y + k1);
        double k3 = h*f(x + h/2.0, y + (3.0*k1 + k2)/8.0);
        double k4 = h*f(x + (2.0*h)/3.0, y + (8.0*k1 + 2.0*k2 + 8.0*k3)/27.0);
        double k5 = h*f(x + (7.0 - sqrt21)*h/14.0, y + (3.0*(3.0*sqrt21 - 7.0)*k1 - 8.0*(7.0 - sqrt21)*k2 + 48.0*(7.0 - sqrt21)*k3 - 3.0*(21.0 - sqrt21*k4))/392.0);
        double k6 = h*f(x + (7.0 + sqrt21)*h/14.0, y + (-5.0*(231.0 + 51.0*sqrt(21.0))*k1 - 40.0*(7.0 + sqrt21)*k2 - 320.0*sqrt21*k3 + 3.0*(21.0 + 121.0*sqrt21)*k4 + 392.0*(6.0 + sqrt21)*k5)/1960.0);
        double k7 = h*f(x + h, y + (15.0*(22.0 + 7.0*sqrt21)*k1 + 120.0*k2 + 40.0*(7.0*sqrt21 - 5.0)*k3 - 63.0*(3.0*sqrt21 - 2.0)*k4 - 14.0*(49.0 + 9.0*sqrt21)*k5 + 70.0*(7.0 - sqrt21)*k6)/180.0);
        double yn = y + (9.0*k1 + 64.0*k3 + 49.0*k5 + 49.0*k6 + 9.0*k7)/180.0;
        result.push_back(make_pair(i, yn));
    }
    return result;
}
void writeData(vector<pair<double, double> > & result)
{
    FILE * output = fopen("data.dat", "w");
    for(pair<double, double> t : result)
    {
        fprintf(output, "%lf %lf %lf\n", t.first, t.second, y_check(t.first));
    }
    fclose(output);
}
int main()
{
    double l, r;
    int n;
    printf("enter l and r \n");
    scanf("%lf %lf", &l, &r);
    printf("enter n \n");
    scanf("%d", &n);
    double y_0 = 0;
    printf("enter y_0 \n");
    scanf("%lf", &y_0);
    double h = (r - l)/n;
    vector<pair<double, double> > result = rungeKutta6_1(l, r, h, n, y_0);
    writeData(result);
    return 0;
}