#include "help.h"
#include "solve.h"

const double EPS = 0.000001;
const double tol = 0.000001;
int main()
{
    pair trash;
    int i;
    double l, r, y01, y02;
    vector t, s1, s2;
    printf("enter l, r, y01, y02 : ");
    scanf("%lf %lf %lf %lf", &l, &r, &y01, &y02);
    init(&t, 10);
    init(&s1, 10);
    init(&s2, 10);
    b = (double *)malloc(sizeof(double)*7);
    c = (double *)malloc(sizeof(double)*7);
    a = (double **)malloc(sizeof(double *)*6);
    if(!a || !c || !b)
    {
        printf("Cannot allocate memory \n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i < 6; i++)
    {
        a[i] = (double *)malloc(sizeof(double)*(i + 1));
        if(!a[i])
        {
            printf("Cannot allocate memory \n");
            exit(EXIT_FAILURE);
        }
    }
    fillArrays();
    solve_fixed(10000, l, r, y01, y02, &t, &s1, &s2);
    for(i = 1 ; i < 100; i++)
        trash = find_val(i, &t, &s1, &s2);
    printf("%lf %lf\n", trash.y1, trash.y2);
    //solve_change(l, r, y01, y02, tol, &t, &s1, &s2);
    printf("norm full %17g", norm_full(&t, &s1));
    write_data_1(&t, &s1, &s2);
    free_vec(&t);
    free_vec(&s1);
    free_vec(&s2);
    for(i = 0; i < 6; i++)
        free(a[i]);
    free(a);
    free(b);
    free(c);
    return 0;
}