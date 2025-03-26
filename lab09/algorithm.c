#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define NX 40
#define NY 40
#define N (NX + 1) * (NY + 1)
#define DELTA 1.0
#define DT 1.0
#define TA 40
#define TB 0
#define TC 30
#define TD 0
#define KB 0.1
#define KD 0.6
#define IT_MAX 2000

gsl_matrix* A;
gsl_matrix* B;
gsl_vector* c;
gsl_vector* T;
gsl_permutation *p;
gsl_vector* d;

char filename[11];

int l;
FILE *tfile;
FILE *lfile;


void wypelnienie_macierzy() {
    // wnętrze obszaru
    for (int i = 1; i < NX; i++) {
        for (int j = 1; j < NY; j++) {
            l = i + j * (NX + 1);
            gsl_matrix_set(A,l,l-NX-1,DT/(2*DELTA*DELTA));
            gsl_matrix_set(A,l,l-1,DT/(2*DELTA*DELTA));
            gsl_matrix_set(A,l,l+1,DT/(2*DELTA*DELTA));
            gsl_matrix_set(A,l,l+NX+1,DT/(2*DELTA*DELTA));
            gsl_matrix_set(A,l,l,-2*DT/(DELTA*DELTA)-1);

            gsl_matrix_set(B,l,l-NX-1,-DT/(2*DELTA*DELTA));
            gsl_matrix_set(B,l,l-1,-DT/(2*DELTA*DELTA));
            gsl_matrix_set(B,l,l+1,-DT/(2*DELTA*DELTA));
            gsl_matrix_set(B,l,l+NX+1,-DT/(2*DELTA*DELTA));
            gsl_matrix_set(B,l,l,2*DT/(DELTA*DELTA)-1);
            
            printf("%f\t%f\n", DT/(2*DELTA*DELTA), gsl_matrix_get(A,l,l-1));
            // printf("A[%d][%d] = %f, B[%d][%d] = %f\n", i, j, A[i][j], i, j, B[i][j]);
            // printf(" %.2f ", A[i][j]);
        }
        // printf("\n");
    }


    for(int i=1;i<NX;i++){
        int l=i+NY*(NX+1);
        gsl_matrix_set(A,l,l-NX-1,-1.0/(KB*DELTA));
        gsl_matrix_set(A,l,l,1.0+1.0/(KB*DELTA));
        gsl_vector_set(c,l,TB);
    }

    for(int i=1;i<NX;i++){
        int l=i+0*(NX+1);
        gsl_matrix_set(A,l,l,1.0+1.0/(KD*DELTA));
        gsl_matrix_set(A,l,l+NX+1,-1.0/(KD*DELTA));
        gsl_vector_set(c,l,TD);
    }

    for(int i=0;i<=NX;i+=NX){
        for(int j=0;j<=NY;j++){
            int l=i+j*(NX+1);
            gsl_matrix_set(A,l,l,1);
            gsl_matrix_set(B,l,l,1);
            gsl_vector_set(c,l,0);
        }
    }

}

void warunki_poczatkowe() {
    for(int j=0;j<=NY;j++){
        gsl_vector_set(T,j*(NX+1),TA);
        gsl_vector_set(T,NX+j*(NX+1),TC);
    }
}

void rozklad_lu() {
    p = gsl_permutation_alloc(N);
    int signum = 0;
    gsl_linalg_LU_decomp(A, p, &signum);
}

void zapis_t(int it) {
    snprintf(filename, 11, "t_%d.dat", it);
    tfile = fopen(filename, "w");
    for (int i = 0; i <= NX; i++) {
        for (int j = 0; j <= NY; j++) {
            l = i + j * (NX + 1);
            fprintf(tfile, "%f\n", gsl_vector_get(T,l));
        }
    }
    printf("Saved to %s\n", filename);
    fclose(tfile);
}

void zapis_l(int it) {
    snprintf(filename, 11, "l_%d.dat", it);
    lfile = fopen(filename, "w");
    for(int i=1;i<NX;i++){
        for(int j=1;j<NY;j++){
            fprintf(lfile, "%.10f\n", gsl_vector_get(T,i-1+j*(NX+1))+gsl_vector_get(T,i+1+j*(NX+1))+gsl_vector_get(T,i+(j-1)*(NX+1))+gsl_vector_get(T,i+(j+1)*(NX+1))-4*gsl_vector_get(T,i+j*(NX+1)));
        }
    }
    printf("Saved to %s\n", filename);
    fclose(tfile);
}

void algorytm_CN() {
    wypelnienie_macierzy();
    warunki_poczatkowe();
    for (int i = 0; i < N; i++) {
        for (int j= 0; j < N; j++) {
            printf("%.1f ", gsl_matrix_get(A, i, j));
        }
        printf("\n");
    }
    rozklad_lu();

    int it_plot[] = {100, 200, 500, 1000, 2000};
    int size_it_plot = sizeof(it_plot) / sizeof(*it_plot);
    for (int it = 0; it <= IT_MAX; it++) {

        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0.0, d);
        gsl_blas_daxpy(1.0,c,d);

        // rozwiaż ukł. równ (LU): A * T = vec_d // T = T ^ (n + 1)
        gsl_linalg_LU_solve(A,p,d,T);


        // printf("%d/%d\n", it, IT_MAX);

        for (int idx = 0; idx < size_it_plot; idx++) {
            if (it == it_plot[idx]) {
                zapis_t(it);
                zapis_l(it);
            }
        }
    }

    gsl_vector_free(d);
    gsl_permutation_free(p);
    gsl_vector_free(T);
    gsl_vector_free(c);
    gsl_matrix_free(B);
    gsl_matrix_free(A);
}

int main() {
    A = gsl_matrix_alloc(N,N);
    B = gsl_matrix_alloc(N,N);
    c = gsl_vector_alloc(N);
    T=gsl_vector_alloc(N);
    d=gsl_vector_alloc(N);
    algorytm_CN();
}