
#include<stdlib.h>
#include<stdio.h>

double energy(double* x, double *dE, int N, int d)
{

    double E, s, J, dJ;
    double r[d];

    E = 0.0;

    for (int n = 0; n < N-1; n++) {
        for (int k = n+1; k < N; k++) {
            s = 0.0;
            for (int i = 0; i < d; i++) {
                r[i] = x[k*d+i] - x[n*d+i];
                s += r[i]*r[i];
            }
            E += 1./(s*s*s*s*s*s) - 2. / (s*s*s);
            dJ = -12.*(1./(s*s*s*s*s*s*s) - 1./(s*s*s*s));
            for (int i = 0; i < d; i++) {
                dE[k*d+i] += dJ * r[i];
                dE[n*d+i] -= dJ * r[i];
            }
        }
    }

    return E;
}
