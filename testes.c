#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

void printIt(double *m, int t)
{
    for (int i = 1; i <= t; i++)
        printf("%.1lf  ", m[i]);
    printf("\n");
}

int main(int argc, char const *argv[])
{
    int dimensao = 10;
    double valor = 0;

    double *solido = (double*) malloc(dimensao * sizeof(double));
    double *solido2 = (double*) malloc(dimensao * sizeof(double));

    for (int i = 0; i < dimensao; ++i)
    {
        solido[i] = 10;
        solido2[i] = 10;
    }

    #pragma omp parallel
    {
        int tID = omp_get_thread_num();
        printf("thread = %d\n", tID);

        
        #pragma omp for schedule(dynamic)
        for (int i = 1; i <= dimensao; i++)
        {
            if (0 < i && i < dimensao-1)
            {
                #pragma omp critical
                {
                    solido2[i] = (2 * solido[i]) - ((solido[i-1] + solido[i+1])/4); // Equação
                    printf("thread(%d)[%d] - for(calculo)\n", tID, i);
                    printf("\t%.1lf - (%.1lf/4)\n", (2 * solido[i]), (solido[i-1] + solido[i+1]));
                }
            }
        }
        #pragma omp barrier


        #pragma omp single
        {
            printIt(solido2, dimensao);
            
            // #pragma omp for reduction(+: valor)
            for (int i = 1; i <= dimensao; i++)
            {
                valor += abs(solido[i] - solido2[i]);
                printf("thread(%d) - for(valor)\n", tID);
            }
        }
    }
}