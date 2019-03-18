#include <stdio.h>
#include <math.h>
typedef struct
{
    double y1;
    double y2;
}pair;
pair make_pair(double y1, double y2)
{
    pair tmp;
    tmp.y1 = y1;
    tmp.y2 = y2;
    return tmp;
}
pair add(int n, void * pairs)
{
    int i;
    pair out;
    pair * tmp = (pair *)pairs;
    for(i = 0; i < n; i++)
    {
        pairs[i].
    }
}
int main()
{
    pair t;
    t.y1 = 10;
    t.y2 = 100;
}