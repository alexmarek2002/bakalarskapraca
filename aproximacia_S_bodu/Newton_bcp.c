#include <stdio.h>
#include <math.h>

double a, c, ell;

double f(double x){
    return 2.0*x/sqrt(a*a + x*x) - (ell - x)/sqrt(c*c + (ell - x)*(ell - x));
}

double df(double x){
    return 2.0*a*a/pow(a*a + x*x, 1.5) + c*c/pow(c*c + (ell - x)*(ell - x), 1.5);
}

double newton(double x0){
    double x = x0;
    double eps = 1e-10;
    int max_i = 100, i = 0;

    while (fabs(f(x)) > eps && i < max_i){
        x = x - f(x)/df(x);
        i++;
    }
    return x;
}

int main(void){
    a=0.604152;
    c=1.7297;
    ell=2.23454;

    double x0 = 0.0;              // tvoj start
    double x = newton(x0);

    printf("x = %.12f\n", x);
    printf("f(x) = %.3e\n", f(x));
    return 0;
}
