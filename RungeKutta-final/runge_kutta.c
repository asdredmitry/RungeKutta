#include <stdio.h>
#include <math.h>
#include <stdlib.h>
const double alpha = -0.2;
double f(double t, double x, double y) {

    return y;
}

double g(double t, double x, double y) {
    return -(1 + alpha*x*x)*x + cos(t);
    //return (1. + 3. * cos(2.*t)) / 2. - y - x * fabs(x);
    //return -1*x;
}

void step(double t, double x, double y, int i, double h, double *res) {
    double k1 = f(t, x, y);
    double m1 = g(t, x, y);

    double k2 = f(t + 0.5 * h, x + h * 0.5 * k1, y + h * 0.5 * m1);
    double m2 = g(t + 0.5 * h, x + h * 0.5 * k1, y + h * 0.5 * m1);

    double k3 = f(t + 2./3. * h, x + h * (2./9. * k1 + 4./9. * k2), y + h * (2./9. * m1 + 4./9. * m2));
    double m3 = g(t + 2./3. * h, x + h * (2./9. * k1 + 4./9. * k2), y + h * (2./9. * m1 + 4./9. * m2));

    double k4 = f(t + 1./3. * h, x + h * (7./36. * k1 + 2./9. * k2 - 1./12. * k3), y + h * (7./36. * m1 + 2./9. * m2 - 1./12. * m3));
    double m4 = g(t + 1./3. * h, x + h * (7./36. * k1 + 2./9. * k2 - 1./12. * k3), y + h * (7./36. * m1 + 2./9. * m2 - 1./12. * m3));

    double k5 = f(t + 5./6. * h,
                    x + h * (-35./144. * k1 - 55./36. * k2 + 35./48. * k3 + 15./8. * k4),
                    y + h * (-35./144. * m1 - 55./36. * m2 + 35./48. * m3 + 15./8. * m4)
                );
    double m5 = g(t + 5./6. * h,
                    x + h * (-35./144. * k1 - 55./36. * k2 + 35./48. * k3 + 15./8. * k4),
                    y + h * (-35./144. * m1 - 55./36. * m2 + 35./48. * m3 + 15./8. * m4)
                );

    double k6 = f(t + 1./6. * h,
                    x + h * (-1./360. * k1 - 11./36. * k2 - 1./8. * k3 + 1./2. * k4 + 1./10. * k5),
                    y + h * (-1./360. * m1 - 11./36. * m2 - 1./8. * m3 + 1./2. * m4 + 1./10. * m5)
                );
    double m6 = g(t + 1./6. * h,
                    x + h * (-1./360. * k1 - 11./36. * k2 - 1./8. * k3 + 1./2. * k4 + 1./10. * k5),
                    y + h * (-1./360. * m1 - 11./36. * m2 - 1./8. * m3 + 1./2. * m4 + 1./10. * m5)
                );

    double k7 = f(t + h,
                    x + h * (-41./260. * k1 + 22./13. * k2 + 43./156. * k3 - 118./39. * k4 + 32./195. * k5 + 80./39. * k6),
                    y + h * (-41./260. * m1 + 22./13. * m2 + 43./156. * m3 - 118./39. * m4 + 32./195. * m5 + 80./39. * m6)
                );
    double m7 = g(t + h,
                    x + h * (-41./260. * k1 + 22./13. * k2 + 43./156. * k3 - 118./39. * k4 + 32./195. * k5 + 80./39. * k6),
                    y + h * (-41./260. * m1 + 22./13. * m2 + 43./156. * m3 - 118./39. * m4 + 32./195. * m5 + 80./39. * m6)
                );

    res[0] = t + h; 
    res[1] = x + h * (13./200. * k1 + 0. * k2 + 11./40. * k3 + 11./40. * k4 + 4./25. * k5 + 4./25. * k6 + 13./200. * k7);
    res[2] = y + h * (13./200. * m1 + 0. * m2 + 11./40. * m3 + 11./40. * m4 + 4./25. * m5 + 4./25. * m6 + 13./200. * m7);
}

double error(double *res_h_second, double *res_2h) {
    return fabs(res_h_second[1] - res_2h[1]) / (pow(2., 6.) - 1);
}

double control(double h, double err, double tol) {
    return h * fmin(1.5, fmax(0, pow((tol/err), 1./(6. + 1.))));
    //return h;
}

int main() {
    double N = 100;
    double b = 3.9;
    double tol = 0.0000000001;
    double h = 0.1;

    double err = 0;
    double eps = 0.000000001;

    FILE *fout = fopen("res.txt", "w");

    double *T = (double*)malloc(sizeof(double) * N);
    double *X = (double*)malloc(sizeof(double) * N);
    double *Y = (double*)malloc(sizeof(double) * N);
    double *res_h_first = (double*)malloc(sizeof(double) * 3);
    double *res_h_second = (double*)malloc(sizeof(double) * 3);
    double *res_2h = (double*)malloc(sizeof(double) * 3);
    //double *tmp_T;
    //double *tmp_X;
    //double *tmp_Y;

    int i = 0;

    T[0] = 0;
    X[0] = 0;
    Y[0] = 1;

    while ((fabs(T[i] + h - b) > eps) && (T[i] + h - b < eps)) {
        printf("%f\n", T[i]);
        if (i + 3 >= N) {
            N *= 2;
            //tmp_T = T;
            //tmp_X = X;
            //tmp_Y = Y;
            T = (double*)realloc(T, sizeof(double) * N);
            X = (double*)realloc(X, sizeof(double) * N);
            Y = (double*)realloc(Y, sizeof(double) * N);
            //free(tmp_T);
            //free(tmp_X);
            //free(tmp_Y);
        }

        step(T[i], X[i], Y[i], i, h, res_h_first);
        step(res_h_first[0], res_h_first[1], res_h_first[2], i, h, res_h_second);
        step(T[i], X[i], Y[i], i, 2 * h, res_2h);

        err = error(res_h_second, res_2h);
        if (err <= tol) {
            i++;

            T[i] = res_h_first[0];
            X[i] = res_h_first[1];
            Y[i] = res_h_first[2];

            i++;

            T[i] = res_h_second[0];
            X[i] = res_h_second[1];
            Y[i] = res_h_second[2];

        } else {
            h = control(h, err, tol);
            step(T[i], X[i], Y[i], i, h, res_h_first);

            i++;

            T[i] = res_h_first[0];
            X[i] = res_h_first[1];
            Y[i] = res_h_first[2];
        }
    }

    for (int j = 0; j <= i; j++) {
        fprintf(fout, "%f %f\n", T[j], X[j]);
    }

    free(T);
    free(X);
    free(Y);
    free(res_h_first);
    free(res_h_second);
    free(res_2h);
    fclose(fout);

    return 0;
}