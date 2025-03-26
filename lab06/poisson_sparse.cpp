#include <iomanip>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <cmath>
#include "mgmres.h" 

const double delta = 0.1;
int nx, ny; 
double xmax, ymax, sigma;

double epsilon(int l, double eps1, double eps2) {
    int j = l / (nx + 1);
    int i = l - j * (nx + 1);
    return (i > nx / 2) ? eps2 : eps1;
}

double rho1(long double x, long double y) {
    double exp_arg = -pow(x - 0.25 * xmax, 2.0) / pow(sigma, 2.0) - pow(y - 0.5 * ymax, 2.0) / pow(sigma, 2.0);

    if (exp_arg < -745) {  
        printf("Wykladnik zbyt mały\n");
        return 0.0L;        
    }

    return exp(exp_arg);
}

double rho2(double x, double y) {
    return -exp(-pow(x - 0.75 * xmax, 2) / pow(sigma, 2) - pow(y - 0.5 * ymax, 2) / pow(sigma, 2));
}

void algorytm_wypelniania(double V1, double V2, double V3, double V4, int nx, int ny, double a[], int ja[], int ia[], 
                          double b[], double V[], int *nz_num, bool rho_zero, double eps1, double eps2) {
    int N = (nx + 1) * (ny + 1);
    int k = -1;
    *nz_num = 0;

    std::ofstream test_file;
    if (nx == 4) {
        test_file.open("test.txt");
        test_file << "l\ti_l\tj_l\tb[l]\n";
    }

    for (int l = 0; l < N; l++) {
        int j = l / (nx + 1);
        int i = l % (nx + 1);

        int brzeg = 0; // Boundary indicator: 0 - inside; 1 - boundary
        double vb = 0.0;

        // Boundary conditions
        if (i == 0) { brzeg = 1; vb = V1; }
        if (j == ny) { brzeg = 1; vb = V2; } 
        if (i == nx) { brzeg = 1; vb = V3; } 
        if (j == 0) { brzeg = 1; vb = V4; }

        // Fill RHS vector
        double x = delta * i;
        double y = delta * j;
        b[l] = (rho_zero) ? 0.0 : -(rho1(x, y) + rho2(x, y));
        if (brzeg == 1) b[l] = vb;

        ia[l] = -1;

        // Fill sparse matrix
        if (l - (nx + 1) >= 0 && brzeg == 0) {
            k++;
            if (ia[l] < 0) ia[l] = k;
            a[k] = epsilon(l - (nx + 1), eps1, eps2) / (delta * delta);
            ja[k] = l - (nx + 1);
        }
        if (l - 1 >= 0 && brzeg == 0) {
            k++;
            if (ia[l] < 0) ia[l] = k;
            a[k] = epsilon(l, eps1, eps2) / (delta * delta);
            ja[k] = l - 1;
        }
        k++;
        if (ia[l] < 0) ia[l] = k; 
        a[k] = (brzeg == 0) 
             ? -(2 * epsilon(l, eps1, eps2) + epsilon(l + 1, eps1, eps2) + epsilon(l + (nx + 1), eps1, eps2)) / (delta * delta)
             : 1.0;
        ja[k] = l;

        if (l < N && brzeg == 0) {
            k++;
            a[k] = epsilon(l + 1, eps1, eps2) / (delta * delta);
            ja[k] = l + 1;
        }
        if (l < N - nx - 1 && brzeg == 0) {
            k++;
            a[k] = epsilon(l + (nx + 1), eps1, eps2) / (delta * delta);
            ja[k] = l + (nx + 1);
        }

        if (nx == 4) test_file << l << "\t" << i << "\t" << j << "\t" << b[l] << "\n";
    }

    *nz_num = k + 1;
    ia[N] = *nz_num;

    if (nx == 4) {
        test_file << "\nk\ta[k]\n";
        for (int count = 0; count <= k; count++) {
            test_file << count << "\t" << a[count] << "\n";
        }
        test_file.close();
    }
}

void compute_potential(double V1, double V2, double V3, double V4, int grid_size, const char *output_file, bool rho_zero, double eps1, double eps2) {
    nx = ny = grid_size;
    xmax = delta * nx;
    ymax = delta * ny;
    sigma = xmax / 10;

    int N = (nx + 1) * (ny + 1);

    double *a = (double *)malloc(5 * N * sizeof(double));
    int *ja = (int *)malloc(5 * N * sizeof(int));
    int *ia = (int *)malloc((N + 1) * sizeof(int));
    double *b = (double *)malloc(N * sizeof(double));
    double *V = (double *)malloc(N * sizeof(double));
    int nz_num;

    for (int i = 0; i < 5 * N; i++) {
        a[i] = 0.0;
    }
    for (int i = 0; i < N + 1; i++) {
        ia[i] = -1;
    }
    for (int i = 0; i < N; i++) {
        b[i] = 0.0;
        V[i] = 0.0;
    }

    algorytm_wypelniania(V1, V2, V3, V4, nx, ny, a, ja, ia, b, V, &nz_num, rho_zero, eps1, eps2);

    int itr_max = 500;
    int mr = 500;
    double tol_abs = 1e-8, tol_rel = 1e-8;
    pmgmres_ilu_cr(N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

    std::ofstream output(output_file);
    for (int l = 0; l < N; l++) {
        int j = l / (nx + 1);
        int i = l - j * (nx + 1);
        output << i << " " << j << " " << std::fixed << std::setprecision(6) << V[l] << std::endl;
    }
    output.close();
    
    free(a);
    free(ja);
    free(ia);
    free(b);
    free(V);

    printf("Potential map for nx = ny = %d saved to '%s'.\n", grid_size, output_file);
}

int main() {
    double eps1 = 1.0, eps2 = 1.0;
    double V1 = 10.0, V2 = -10.0, V3 = 10.0, V4 = -10.0;
    compute_potential(V1, V2, V3, V4, 4, "test.dat", true, eps1, eps2);

    compute_potential(V1, V2, V3, V4, 50, "potential_50.dat", true, eps1, eps2);
    compute_potential(V1, V2, V3, V4, 100, "potential_100.dat", true, eps1, eps2);
    compute_potential(V1, V2, V3, V4, 200, "potential_200.dat", true, eps1, eps2);

    V1 = V2 = V3 = V4 = 0.0;
    compute_potential(V1, V2, V3, V4, 100, "potential_eps1.dat", false, eps1, eps2);
    eps2 = 2.0;
    compute_potential(V1, V2, V3, V4, 100, "potential_eps2.dat", false, eps1, eps2);
    eps2 = 10.0;
    compute_potential(V1, V2, V3, V4, 100, "potential_eps10.dat", false, eps1, eps2);

    return 0;
}
