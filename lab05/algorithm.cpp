#include <cmath>
#include <omp.h>
#include <string>
#include <cstdio>

const double delta = 0.2;
const int n_x = 128;
const int n_y = 128;
const double x_max = delta * n_x;
const double y_max = delta * n_y;
const double TOL = 1e-8;
int counter = 0;

inline double pow2(double x){
    return x * x;
}

inline double VB1(double y){
    return sin(M_PI * y / y_max);
}

inline double VB2(double x){
    return -sin(2 * M_PI * x / x_max);
}

inline double VB3(double y){
    return sin(M_PI * y / y_max);
}

inline double VB4(double x){
    return sin(2 * M_PI * x / x_max);
}

void setupBorders(double V[n_x + 1][n_y + 1]){
    #pragma omp parallel for
    for(int i=0; i<=n_x; i++){
        V[i][0] = VB4(i * delta);
        V[i][n_y] = VB2(i * delta);
    }
    #pragma omp parallel for
    for(int j=0; j<=n_y; j++){
        V[0][j] = VB1(j * delta);
        V[n_x][j] = VB3(j * delta);
    }
}

void saveToFile(double V[n_x + 1][n_y + 1], int k){
    FILE* file = fopen(("map" + std::to_string(k) + ".txt").c_str(), "w");
    for(int j=0; j <=n_y; j+=k){
        for(int i=0; i<=n_x; i+=k){
            fprintf(file, "%e ", V[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);
}

double calcS(double V[n_x + 1][n_y + 1], int k) {
    double S = 0;
    #pragma omp parallel for collapse(2) reduction(+:S)
    for(int i=0; i<=n_x - k; i+=k){
        for(int j=0; j<=n_y - k; j+=k){
            S += pow2(k * delta)/2 * (pow2((V[i+k][j]-V[i][j])/(2*k*delta) + (V[i+k][j+k]-V[i][j+k])/(2*k*delta)) 
            + pow2((V[i][j+k]-V[i][j])/(2*k*delta) + (V[i+k][j+k]-V[i+k][j])/(2*k*delta)));
        }
    }
    return S;
}

double calcSTOL(double S_new, double S_old){
    return fabs((S_new - S_old)/S_old);
}

void tightNetwork(double V[n_x + 1][n_y + 1], int k){
    #pragma omp parallel for collapse(2)
    for(int i=0; i<=n_x-k; i+=k){
        for(int j=0; j<=n_y-k; j+=k){
            V[i + k/2][j + k/2] = 0.25 * (V[i][j] + V[i+k][j] + V[i][j+k] + V[i+k][j+k]);
            V[i + k][j + k/2] = 0.5 * (V[i + k][j] + V[i+k][j+k]);
            V[i + k/2][j+k] = 0.5 * (V[i][j + k] + V[i+k][j+k]);
            V[i + k/2][j] = 0.5 * (V[i][j] + V[i+k][j]);
            V[i][j+k/2] = 0.5 * (V[i][j] + V[i][j+k]);
        }
    }
    setupBorders(V);
}

void performRelaxasion(double V[n_x + 1][n_y + 1], int k){
    double S_old = calcS(V, k);
    double S_new = S_old;
    FILE* file = fopen(("Sit" + std::to_string(k) + ".txt").c_str(), "w");
    counter++;

    while(true){
        S_old = S_new;
        #pragma omp parallel for collapse(2)
        for(int i=k; i<=n_x-k; i+=k){
            for(int j=k; j<=n_y-k; j+=k){
                V[i][j] = 0.25 * (V[i+k][j] + V[i-k][j] + V[i][j+k] + V[i][j-k]);
            }
        }
        S_new = calcS(V, k);
        fprintf(file, "%i %e\n", counter, S_new);
        counter++;
        if(calcSTOL(S_new, S_old) < TOL) break;
    }

    fclose(file);
}



int main(){
    double V[n_x + 1][n_y + 1];
    int k = 16;
    setupBorders(V);
    #pragma omp parallel for collapse(2)
    for(int i=1; i<=n_x-1; i++){
        for(int j=1; j<=n_y-1; j++){
            V[i][j] = 0;
        }
    }
    while(k){
        performRelaxasion(V, k);
        saveToFile(V, k);
        if(k>1){
            tightNetwork(V, k);
        }
        k /= 2;
    }

    
    return 0;
}