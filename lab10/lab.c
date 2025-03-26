#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 150
#define NT 1000
#define DELTA 0.1
#define DT 0.05
#define XA 7.5
#define SIGMA 0.5

double u0[NX + 1];
double u[NX + 1];
double v[NX + 1];
double vp[NX + 1];
double a[NX + 1];
double E;
FILE *ufile, *Efile;

void warunki_brzegowe() {
    u[0] = u[NX] = 0;
    v[0] = v[NX] = 0;
}

void warunki_poczatkowe(bool warunekPoczatkowy) {
    for(int i = 0; i <= NX; i++) {
        double x = i * DELTA;
        if (warunekPoczatkowy) {
            u[i] = v[i] = 0;
        } else {
            u[i] = exp(-pow(x - XA, 2) / (2 * SIGMA * SIGMA));
            v[i] = 0;
        }
    }
}

double delta_kroneckera(double x, double xf) {
    return (x == xf) ? 1 : 0;
}

double a_F(double x, double t, double xf) {
    if (xf == -1.0)
        return cos( (50 * t) / (DT * NT) );
    else
        return cos( (50 * t) / (DT * NT) ) * delta_kroneckera(x, xf);
}

void wartosc_a(double alpha, double beta, double t, double xf) {
    double x;
    for (int i = 0; i < NX; i++) {
        x = DELTA * i;
       a[i] = (u[i + 1] - 2 * u[i] + u[i - 1]) / (DELTA * DELTA) -
              beta * (u[i] - u0[i]) / DT +
              alpha * a_F(x, t, xf);
    }
}

void energia() {
    double sum = 0.0;
    for (int i = 1; i < NX; i++) {
        sum += v[i] * v[i] + pow( (u[i + 1] - u[i - 1]) / (2 * DELTA), 2 );
    }
    E = (DELTA / 4) * ( pow( (u[1] - u[0]) / DELTA, 2 ) + pow( (u[NX] - u[NX - 1] ) / DELTA, 2 ) ) + DELTA / 2 * sum;
}

char* stworz_plik(char symbol, double alpha, double beta, double xf) {
    char *filename = malloc(30 * sizeof(char));
    if (xf != -1.0)
        sprintf(filename, "%c_a_%.2f_b_%.2f_xf_%.2f.dat", symbol, alpha, beta, xf);
    else 
        sprintf(filename, "%c_a_%.2f_b_%.2f.dat", symbol, alpha, beta);

    if (symbol == 'u') {
        ufile = fopen(filename, "w");
        fclose(ufile);
    } else {
        Efile = fopen(filename, "w");
        fclose(Efile);
    }

    return filename;
}

void zapis_E(char* filename, int n) {
    static int counter = 0;
    if (n == 1) counter = 0;
    Efile = fopen(filename, "a");
    fprintf(Efile, "%f\n", E);
    fclose(Efile);
    printf("Saved E, counter: %d\n", counter);
    counter++;
}

void zapis_u(char* filename, int n) {
    static int counter = 0;
    if (n == 1) counter = 0;
    ufile = fopen(filename, "a");
    for (int i = 0; i <= NX; i++) {
        fprintf(ufile, "%f\n", u[i]);
    }
    fclose(ufile);
    printf("Saved u, counter: %d\n", counter);
    counter++;
}

void algorytm(double alpha, double beta, double xf, double warunekPoczatkowy) {
    char* E_filename = stworz_plik('e', alpha, beta, xf);
    char* u_filename = stworz_plik('u', alpha, beta, xf);
    warunki_brzegowe();
    warunki_poczatkowe(warunekPoczatkowy);

    // zachowanie poprzedniego wyniku
    for (int i = 0; i <= NX; i++) u0[i] = u[i];

    double t = 0.0;
    // inicjalizacja a[]
    wartosc_a(alpha, beta, t, xf);

    for (int n = 1; n <= NT; n++) {
        t = n * DT;
        for (int i = 0; i <= NX; i++) {
            vp[i] = v[i] + DT / 2 * a[i];
            u0[i] = u[i];
            u[i] = u[i] + DT * vp[i];
        }
        for (int i = 0; i <= NX; i++) {
            wartosc_a(alpha, beta, t, xf);
        }
        for (int i = 0; i <= NX; i++) {
            v[i] = vp[i] + DT / 2 * a[i];
        }
        energia();
        
        zapis_E(E_filename, n);
        zapis_u(u_filename, n);
    }
}

int main() {
    algorytm(0, 0, -1.0, false);
    algorytm(0, 0.1, -1.0, false);
    algorytm(0, 1, -1.0, false);
    algorytm(1, 1, 2.5, true);
}