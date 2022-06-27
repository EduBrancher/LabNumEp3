#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_evaluation(double x, double fx, char* func){
    printf("O valor do polinomio para %s avaliado em %.4lf foi %.4lf.\n", func, x, fx);
}

double evaluate_newton_polynomial(double x, double* xarr, int length, double* coefs){
    double eval = 0;
    double coefficient = 0;
    double xval = 0;
    for (int i = 0; i < length; i++){
        coefficient = coefs[i];
        xval = 1;
        for (int j = 0; j < i; j++){
            xval = xval*(x - xarr[j]);
        }
        eval = eval + coefficient*xval;
    }
    return eval;
}

void newton_interpolant(double* xarr, double* fxarr, int length, double div_differences[length][length], double coefs[length]){
    for (int i = 0; i < length; i++){
        div_differences[i][0] = fxarr[i];
    }

    for (int j = 1; j < length; j++){
        for (int i = j; i < length; i++){
            double numerator = div_differences[i][j-1] - div_differences[i-1][j-1];
            double denominator = xarr[i] - xarr[i-j];
            div_differences[i][j] = numerator/denominator;
        }
    }

    for (int i = 0; i < length; i++){
        coefs[i] = div_differences[i][i];
    }
}

//calculate area for a given segment with simpson interpolation
double simpson_area(double start, double end, int length, double coefs[length], double xarr[length]){
    double midpoint = (start + end) / 2;
    double farr[] = {evaluate_newton_polynomial(start, xarr, length, coefs),
                     evaluate_newton_polynomial(midpoint, xarr, length, coefs),
                     evaluate_newton_polynomial(end, xarr, length, coefs)};
    return (((end - start)/6)*(farr[0] + 4*farr[1] + farr[2]));

}

//calculate area for a given segment with trapezoidal interpolation
double trapezoidal_area(double start, double end, int length, double coefs[length], double xarr[length]){
    double farr[] = {evaluate_newton_polynomial(start, xarr, length, coefs),
                     evaluate_newton_polynomial(end, xarr, length, coefs)};
    return ((end - start)/2)*(farr[0] + farr[1]);
}

double integrate_simpson(double start, double end, double h, int length, double coefs[length], double xarr[length]){
    int n = round((end-start)/h);
    printf("Integrando Simpson com %d intervalos.\n", n);
    double area = 0;
    double interval_start = start;
    double interval_end = start + h;
    for (int i = 0; i < n; i++){
        area = area + simpson_area(interval_start, interval_end, length, coefs, xarr);
        interval_start = interval_end;
        interval_end = interval_end + h;
    }
    return area;
}

double integrate_trapezoidal(double start, double end, double h, int length, double coefs[length], double xarr[length]){
    int n = round((end-start)/h);
    printf("Integrando trapezio com %d intervalos.\n", n);
    double area = 0;
    double interval_start = start;
    double interval_end = start + h;
    for (int i = 0; i < n; i++){
        area = area + trapezoidal_area(interval_start, interval_end, length, coefs, xarr);
        interval_start = interval_end;
        interval_end = interval_end + h;
    }
    return area;
}

int main(){
    //initialize stuff
    int length = 7;
    double xarr[] = {0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0};
    double fxarr[] = {0.0, 1.5297, 9.5120, 8.7025, 2.8087, 1.0881, 0.3537};
    double coefs[length];
    for (int i = 0; i < length; i++){
        coefs[i] = 0;
    }
    double divided_differences[length][length];
    for (int i = 0; i < length; i++){
        for (int j = 0; j < length; j++){
            divided_differences[i][j] = 0;
        }
    }
    newton_interpolant(xarr, fxarr, length, divided_differences, coefs);
    for (int i = 0; i < length; i++){
        for (int j = 0; j < length; j++){
            divided_differences[i][j] = 0;
        }
    }
    double val;
    //verificar corretude do polinomio
    for (int x = 0; x <= 30; x++){
        val = evaluate_newton_polynomial(x, xarr, length, coefs);
        if (x % 5 == 0){
            printf("------------------------------\n");
        }
        print_evaluation(x, val, "F");
    }

    //experimentos com n variando
    double area;
    double h = 5;
    for (int i = 0; i < 5; i++){
        area = integrate_simpson(0.0, 30.0, 2*h, length, coefs, xarr);
        printf("A area obtida foi %.4lf.\n", area);
        h = h/2;
    }
    h = 5;
    for (int i = 0; i < 10; i++){
        area = integrate_trapezoidal(0.0, 30.0, h, length, coefs, xarr);
        printf("A area obtida foi %.4lf.\n", area);
        h = h/2;
    }
}