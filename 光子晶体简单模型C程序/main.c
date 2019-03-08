#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "svd.c"

#define M 1
#define N 801        
#if (M>N)      
#define MN M
#else
#define MN N
#endif

const int num = 801;
const double pi = 3.1416;
const int NN = 10;


double **dmatrix(int, int, int, int);
double *dvector(int, int);
void svdcmp(double **, int, int, double *, double **);


double _Complex **init(int m, int n) {
    int i;
    double _Complex **Matrix;
    Matrix = (double _Complex**) malloc(sizeof(double _Complex*) * (m+1));
    for(i = 1; i <= m; i++) {
        Matrix[i] = (double _Complex*) malloc(sizeof(double _Complex) * (n+1));
    }
    return Matrix;
}


void print_r(double **a, int m, int n) {
    int i, j;
    for (i = 1; i <= m; i++) {
        for (j = 1; j <= n; j++) {
            printf("%le ", a[i][j]);
        }
        printf("\n");
    }
}

double _Complex **complex_matmul(double _Complex **k1, double _Complex **k2, int r1, int c1, int c2) {
    int i, j, k;
    complex double sum;
    double _Complex **res;

    res = init(r1, c2);
    // res[1][1] = (k1[1][1] * k2[1][1] + k1[1][2] * k2[2][1]);
    // res[1][2] = (k1[1][1] * k2[2][1] + k1[1][2] * k2[2][2]);
    // res[2][1] = (k1[2][1] * k2[1][1] + k1[2][2] * k2[2][1]);
    // res[2][2] = (k1[2][1] * k2[2][1] + k1[2][2] * k2[2][2]);

    for (i = 1; i <= r1; i++) {
        for (j = 1; j <= c2; j++) {
            sum = 0.0 + 0 * _Complex_I;
            for (k = 1; k <= c1; k++) {
                sum += k1[i][k] * k2[k][j];
            }

            // printf("M1[1][1]:%.3f+%.3ei\n", creal(sum), cimag(sum));
            res[i][j] = sum;
        }
    }

    return res;
}

double **matmul(double **k1, double **k2, int r1, int c1, int c2) {
    int i, j, k;
    double sum;
    double **res;

    res = dmatrix(1, MN, 1, MN);
    for (i = 1; i <= r1; i++) {
        for (j = 1; j <= c2; j++) {
            sum = 0.0;
            for (k = 1; k <= c1; k++) {
                sum += k1[i][k] * k2[k][j];
            }
            res[i][j] = sum;
        }
    }
    return res;
}

double **transpose(double **k, int r, int c) {
    int i, j;
    double **res;

    res = dmatrix(1, r, 1, c);
    for (i = 1; i <= r; i++) {
        for (j = 1; j <= c; j++) {
            res[i][j] = k[j][i];
        }
    }
    return res;
}

double inv(double **a, double **u, int m, int n, double w[], double **v, double **a2) {
    int i,j,k;
    double t;
    double t1[MN], t2[MN];

    for (i = 1; i <= M; i++) {
        for (j = 1; j <= N; j++)
            u[i][j] = a[i][j];
    }

    for (i = 1; i <= M; i++) {
        for (j = 1; j <= N; j++)
            u[i][j] = a[i][j];
    }

    // svdcmp(u, M, N, w, v);
    svdcmp(u, MN, MN, w, v);
    for (i = 1; i <= N; i++) {
        for (j = i+1; j <= N; j++) {
            if (w[i] < w[j]) { /* ������ */
                t = w[i];
                w[i] = w[j];
                w[j] = t;
                /* �������U,V������ */
                /* ��U */
                for (k = 1; k <= M; k++) {
                    t1[k] = u[k][i];
                }
                for (k = 1; k <= M; k++) {
                    u[k][i] = u[k][j];
                }
                for (k = 1; k <= M; k++) {
                    u[k][j] = t1[k];
                }

                /* ��V */
                for (k = 1; k <= N; k++) {
                    t2[k] = v[k][i];
                }
                for (k = 1; k <= N; k++) {
                    v[k][i] = v[k][j];
                }
                for (k = 1; k <= N; k++) {
                    v[k][j] = t2[k];
                }
            }
        }
    }


    /* U�MxM�� */
    // printf("=== U ===\n");
    // print_r(u, M, M);

    /* ����M����W�������� */
    double **W;
    W = dmatrix(1, MN, 1, MN);
    // printf("=== W ===\n");
    for (i = 1; i<= M; i++) {
        for (j =1; j <= N; j++) {
            if (i==j) {
                W[i][j] = w[i];
            } else {
                W[i][j] = 0.0;
            }
        }
    }
    // print_r(W, M, N);


    /* V�NxN�� */
    // printf("=== V ===\n");
    // print_r(v, N, N);

    /* ��������U*W*V������a */
    // printf("=== U*W*V' ===\n");
    double **temp;
    double **V, **B, **U;
    V = dmatrix(1, MN, 1, MN);
    U = dmatrix(1, MN, 1, MN);
    temp = dmatrix(1, MN, 1, MN);
    B = dmatrix(1, MN, 1, MN);

    double sum;

    // V = transpose(v, N, N);
    V = v;
    U = transpose(u, M, M);

    W[1][1] = 1 / W[1][1];

    // printf("=== V*(1/W)*U'====\n");
    temp = matmul(V, W, N, N, M);
    B = matmul(temp, U, N, M, M);

    // printf("B\n");

    double res = 0.0;
    for (i = 1; i <= N; i++) {
        res += B[i][1] * a2[1][i];
    }

    return res;
}

double k_div(double **k1, double **k2) {

    double *w;
    double **u, **v;

    u = dmatrix(1, MN, 1, MN);
    w = dvector(1, MN);
    v = dmatrix(1, MN, 1, MN);

    return inv(k1, u, MN, MN, w, v, k2);
}
int main() {

    int i, j, k;


    double n1=2.33;
    double n2=1.50;	// ������ʵ�������
    double h1=200e-9;
    double h2=100e-9;	//����2���Ϊ0
    int A=1;
    int B=0;


    double _Complex **M1, **M2;
    double T[num];
    double **k1, **k2;
//    double t1[MN], t2[MN];

    int wl[num];
    double e1[num], e2[num], k1_val[num], k2_val[num];
    // double M1, M;

    complex double c1, c2;

    for (i = 0; i < num; i++) {
        wl[i] = i + 100;
    }


    // ����
    k1 = dmatrix(1, MN, 1, MN);
    k2 = dmatrix(1, MN, 1, MN);

    for (i = 0; i < num; i++) {
        k1[1][i+1] = 2*pi*n1 / (wl[i]*1e-9);
        e1[i] = 2*pi*n1*h1 / (wl[i]*1e-9);
        k2[1][i+1] = 2*pi*n2 / (wl[i]*1e-9);
        e2[i] = 2*pi*n2*h2 / (wl[i]*1e-9);
    }

    double temp1, temp2;
    temp2 = k_div(k1, k2);
    temp1 = k_div(k2, k1);

    printf("k1/k2:%e\n", temp1);
    printf("k2/k1:%e\n", temp2);
    c1 = 0.5 * _Complex_I * (temp1 + temp2);
    c2 = 0.5 * _Complex_I * (temp1 - temp2);

    printf("c1:%.3f+%.3ei\n", creal(c1), cimag(c1));
    printf("c2:%.3f+%.3ei\n", creal(c2), cimag(c2));

    //M1 = dmatrix(1, 2, 1, 2);
    //M2 = dmatrix(1, 2, 1, 2);
    M1 = init(2,2);
    M2 = init(2,2);

    for (i = 0; i < num; i++) {
        printf("sin(e2[i]): %.3e\t", e2[i]);
        M1[1][1] = cexp( _Complex_I * e1[i]) * (cos(e2[i]) + c1 * sin(e2[i]));
        M1[1][2] = cexp(-_Complex_I * e1[i]) * (-c2 * sin(e2[i]));
        M1[2][1] = cexp( _Complex_I * e1[i]) * c2 * sin(e2[i]);
        M1[2][2] = cexp(-_Complex_I * e1[i]) * (cos(e2[i]) - c1 * sin(e2[i]));
        M2 = M1;

//        printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[1][1]), cimag(M1[1][1]));
        //      printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[1][2]), cimag(M1[1][2]));
        //    printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[2][1]), cimag(M1[2][1]));
        //  printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[2][2]), cimag(M1[2][2]));
        for (k = 0; k < NN-1; k++) {
            M2 = complex_matmul(M2, M1, 2, 2, 2);
        }
        T[i] = pow(cabs((M2[1][1] * M2[2][2] - M2[2][1] * M2[1][2]) / ((A * M2[2][2] - B * M2[1][2]))), 2);
//
////		c2[i] = 0.5*1i*(k1[i][0]/k2[i][0]-k2[i][0]/k1[i][0]);
    }

//    print_r(M2, 2, 2);


    printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[1][1]), cimag(M1[1][1]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[1][2]), cimag(M1[1][2]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[2][1]), cimag(M1[2][1]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M1[2][2]), cimag(M1[2][2]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M2[1][1]), cimag(M2[1][1]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M2[1][2]), cimag(M2[1][2]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M2[2][1]), cimag(M2[2][1]));
    printf("M1[1][1]:%.3f+%.3ei\n", creal(M2[2][2]), cimag(M2[2][2]));
    for (i = 0; i < 100; i++) {
        printf("%.3e\t", T[i]);
    }


    printf("%.3e", cimag(cexp(10 + _Complex_I)));
//	printf("k1(1): %.10f", k1[0]);
    // printf("%d", (3 ** 10));
    return 0;
}


