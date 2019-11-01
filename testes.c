#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

void printIt(double *m, int t)
{
    for (int i = 0; i < t; i++)
        printf("%.1lf  ", m[i]);
    printf("\n");
}

int main(int argc, char const *argv[])
{
    int on = 1;
    int dimensao = 10;
    double valor = 0, Z;

    double *solido = (double*) malloc(dimensao * sizeof(double));
    double *solido2 = (double*) malloc(dimensao * sizeof(double));

    for (int i = 0; i < dimensao; ++i)
    {
        solido[i] = 10;
        solido2[i] = 10;
    }
    solido[0] = 80;
    solido2[0] = 80;

    double fourier = pow(1.0, 2) * (0.1/pow(1, 2));

    // #pragma omp parallel
    {
        int tID = omp_get_thread_num();
        printf("thread = %d\n", tID);

        while(on)
        {
            // #pragma omp for
            for (int i = 1; i < dimensao; i++)
            {
                if (0 < i && i < dimensao-1)
                {
                    // #pragma omp critical
                    {
                        solido2[i] = solido[i] + (fourier * (solido[i+1] - (2 * solido[i]) + solido[i-1])); // Equação
                        printf("thread(%d)[%d] - for(calculo)\n", tID, i);
                        // printf("\t%.1lf - (%.1lf/4)\n", (2 * solido[i]), (solido[i-1] + solido[i+1]));
                    }
                }
            }
            // #pragma omp barrier

            // #pragma omp single
            {
                printIt(solido2, dimensao);
                printIt(solido, dimensao);

                valor = 0;
                
                // #pragma omp for reduction(+: valor)
                for (int i = 1; i < dimensao; i++)
                {
                    valor += abs(solido2[i] - solido[i]);
                    printf("thread(%d) - for(valor)\n", tID);
                }

                Z = valor / (dimensao);

                printf("Z = %lf\n\n", Z);

                if (Z <= 0)
                {
                    printf("%lf < 0\n", Z);
                    on = 0;
                }

                for (int i = 0; i < dimensao; ++i)
                {
                    solido[i] = solido2[i];
                }
            }
        }   
    }

    free(solido);
    free(solido2);

}