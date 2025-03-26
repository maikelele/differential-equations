#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define NX 400
#define NY 90
#define I1 200
#define I2 210
#define J1 50
#define DELTA 0.01
#define SIGMA 0.1
#define XA 0.45
#define YA 0.45
#define IT_MAX 10000

double psi[NX + 1][NY + 1];
double vx[NX + 1][NY + 1]; 
double vy[NX + 1][NY + 1];
double u0[NX + 1][NY + 1];
double u1[NX + 1][NY + 1];

double vmax;
double dt;
double c;
double xsr;
bool write_u = false;


void funkcja_strumienia() {
    FILE *ifile = fopen("psi.txt", "r");
    int i, j;
    double value;
    while (fscanf(ifile, "%d %d %lf", &i, &j, &value) == 3) 
    {
        psi[i][j] = value;
        
    }
    fclose(ifile);
}


void pole_predkosci() {
    
    for (int i = 1; i < NX; i++) {
        for (int j = 1; j < NY; j++) {
            vx[i][j] = (psi[i][j + 1] - psi[i][j - 1]) / (2 * DELTA);
            vy[i][j] = -(psi[i + 1][j] - psi[i - 1][j]) / (2 * DELTA);
        }
    }
    for (int i = I1; i <= I2; i++) {
        for (int j = 0; j <= J1; j++) {
            vx[i][j] = vy[i][j] = 0;
        }
    }
    
    for (int i = 1; i < NX; i++) {
        vx[i][0] = vy[i][NY] = 0;
    }
    
    for (int j = 0; j <= NY; j++) {
        vx[0][j] = vx[1][j];
        vx[NX][j] = vx[NX - 1][j];
    }
    FILE *vxfile = fopen("vx.dat", "w");
    FILE *vyfile = fopen("vy.dat", "w");
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            fprintf(vxfile, "%f\n", vx[i][j]);
            fprintf(vyfile, "%f\n", vy[i][j]);
        }
    }
    fclose(vxfile);
    fclose(vyfile);
}

double calculate_vmax() 
{
    double tmp = 0.0;
    for(int i = 0; i <= NX; i++) 
    {
        for(int j = 0; j <= NY; j++) 
        {
            if(sqrt(vx[i][j] * vx[i][j] + vy[i][j] * vy[i][j]) > tmp) 
            {
                tmp = sqrt(vx[i][j] * vx[i][j] + vy[i][j] * vy[i][j]);
            }
        }
    }
    
    return tmp;
}

void reset()
{
    for(int i = 0; i <= NX; ++i)
    {
        for(int j = 0; j <= NY; ++j)
        {
            u0[i][j] = 0;
            u1[i][j] = 0;
        }
    }
}

double calculate_xsr() {
    double xsr = 0.0;
    double xi;
    for (int i = 0; i <= NX; i++) {
        xi = i * DELTA;
        for (int j = 0; j <= NY; j++) {
            xsr += xi * u0[i][j] * DELTA * DELTA;
        }
    }

    return xsr;
}

double calculate_c() {
    double c = 0.0;
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            c += u0[i][j] * DELTA * DELTA;
        }
    }

    return c;
}

void algorytm(const char *filename1, const char *filename2, double D) {

    FILE *cfile = fopen(filename1, "w");
    FILE *xfile = fopen(filename2, "w");
    FILE *ufile;
    char filename[20];
    int saves_counter = 0;

    double x, y;
    for(int i = 0; i <= NX; ++i)
     {
        for(int j = 0; j <= NY; ++j) 
        {
            x = i * DELTA;
            y = j * DELTA;
            u0[i][j] = (1 / (2 * M_PI * SIGMA * SIGMA)) * exp(-((x - XA) * (x - XA) + (y - YA) * (y - YA)) / (2 * SIGMA * SIGMA));
        }
    }

    int ip;
    int im;
    int jp;
    int jm;
    static int fname = 0;

    for (int it = 1; it <= IT_MAX; ++it) 
    {
        for (int i = 0; i <= NX; i++) 
        {
            for (int j = 0; j <= NY; j++) 
            {
                u1[i][j] = u0[i][j];
            }
        }

        for (int step = 0; step < 20; ++step) {
            for(int i = 0; i <= NX; ++i) 
            {
                for(int j = 1; j < NY ; j++) 
                {
                    if(i >= I1 && i <= I2 && j <= J1)
                    {
                        continue;
                    }

                    ip = (i + 1) % (NX + 1);
                    im = (i - 1 ) >= 0 ? i - 1 : NX;
                    jp = (j + 1) % (NY + 1);
                    jm = (j - 1 ) >= 0 ? j - 1 : NY;

                    u1[i][j] = (1.0 / (1 + 2 * D * dt / (DELTA * DELTA))) * (u0[i][j] - dt / 2 * vx[i][j] * ((u0[ip][j] - u0[im][j]) / (2 * DELTA) + (u1[ip][j] - u1[im][j]) / (2 * DELTA)) - dt / 2 * vy[i][j] * ((u0[i][jp] - u0[i][jm]) / (2 * DELTA) + (u1[i][jp] - u1[i][jm]) / (2 * DELTA)) + dt / 2 * D * (((u0[ip][j] + u0[im][j] + u0[i][jp] + u0[i][jm] - 4*u0[i][j]) / (DELTA * DELTA)) + ((u1[ip][j] + u1[im][j] + u1[i][jp] + u1[i][jm])/(DELTA * DELTA))));
                
                }
            }
        }

        
        if (!(it % 200)) {
            if (!D) {
                snprintf(filename, sizeof(filename), "u0_%d_%d.dat", saves_counter, it);
            } else {
                snprintf(filename, sizeof(filename), "u1_%d_%d.dat", saves_counter, it);
            }
            ufile = fopen(filename, "w");
            write_u = true;
        }

        for (int i = 0; i <= NX; ++i) 
        {
            for (int j = 0; j <= NY; ++j) 
            {
                u0[i][j] = u1[i][j];
                if (write_u) {
                    fprintf(ufile, "%f\n", u0[i][j]);
                    printf("it: %d, filename: %s\n", it, filename);
                }
            }
        }

        printf("it: %d\n", it);
        if (write_u) {
            fclose(ufile);
            write_u = false;
            saves_counter++;
        }         

        c = calculate_c();
        fprintf(cfile, "%f\n", c);
        xsr = calculate_xsr();
        fprintf(xfile, "%f\n", xsr);
    }
    fclose(cfile);
    fclose(xfile);
}

int main() {
    funkcja_strumienia();
    pole_predkosci();
    vmax = calculate_vmax();
    dt = DELTA / (4 * vmax);

    algorytm("c0.dat", "x0.dat", 0.0);
    reset();
    algorytm("c1.dat", "x1.dat", 0.1);

    return 0;
}