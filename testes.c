#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

void printIt(double *m, int tam)
{
    for (int i = 0; i < tam; ++i)
    {
        printf("%.1lf ", m[i]);
    }
    printf("\n");
}

void printItN(double **m, int tam, int prts, int extra, int cores)
{
    int med;
    for (int i = 1; i <= tam; i += prts)
    {
        med = (i <= (cores*prts)) ? prts : extra;
        for (int j = 1; j < med+1; ++j)
        {
            printf("%.1lf\n", m[(int)((i-1)/prts)][j]);
        }
    }
    printf("\n");   
}

int main(int argc, char const *argv[])
{
    clock_t ticks[2];
    ticks[0] = clock();

    int on, tID, nT;
    int dimensao = 1000000,
    d2 = dimensao+2;
    double valor = 0, Z;
    int i;

    double *solido = (double*) malloc(d2 * sizeof(double));
    double *solido2 = (double*) malloc(d2 * sizeof(double));
    double **partes_solido;

    for (i = 0; i < d2; ++i)
    {
        solido[i] = 10;
        solido2[i] = 10;
    }
    solido[0] = 1000;
    solido2[0] = 1000;

    int cores = 8;

    int tamanho_partes = dimensao/cores;
    int tamanho_extra = dimensao%cores;
    int partes = (0 < tamanho_extra) ? cores+1 : cores;
    printf("Cores presentes:%d\ntamanho das partes:%d\ntamanho da parte extra:%d\n", cores, tamanho_partes, tamanho_extra);

    partes_solido = (double**) malloc(partes * sizeof(double*));
    for (i = 0; i < partes; ++i)
    {
        if (i == cores)
        {
            partes_solido[i] = (double*) malloc((tamanho_extra+2) * sizeof(double));
            printf("Tamanho %d para a parte %d.\n", tamanho_extra+2, i);
        }
        else
        {
            partes_solido[i] = (double*) malloc((tamanho_partes+2) * sizeof(double));
            printf("Tamanho %d para a parte %d.\n", tamanho_partes+2, i);
        }
    }

    int med;
    for (i = 1; i <= dimensao; i += tamanho_partes)
    {
        med = (i <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
        for (int j = 0; j < med+2; ++j)
        {
            partes_solido[(int)((i-1)/tamanho_partes)][j] = solido[i+j-1];
            // printf("solido[%d][%d] = %lf\n", (int)((i-1)/tamanho_partes), j, partes_solido[(int)((i-1)/tamanho_partes)][j]);
        }
    }

    double fourier = pow(1.0, 2) * (0.1/pow(1, 2));


    #pragma omp parallel default(none) private(tID, nT, on, dimensao, fourier, Z) shared(cores, partes, partes_solido, solido2, valor, tamanho_partes, tamanho_extra) num_threads(cores)
    {
        dimensao = 1000000;
        fourier = pow(1.0, 2) * (0.1/pow(1, 2));

        on = 1;
        tID = omp_get_thread_num();
        nT = omp_get_num_threads();
        printf("thread %d de %d.\n", tID+1, nT);
        int med;
        int j;
        #pragma omp barrier

        while (on)
        {
            #pragma omp for
            for (int p = 1; p <= dimensao; p += tamanho_partes)
            {
                j = (int)((p-1)/tamanho_partes);
                med = (p <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
                for (int i = 1; i <= med; ++i)
                {
                    solido2[p+i-1] = partes_solido[j][i] + (fourier * (partes_solido[j][i+1] - (2 * partes_solido[j][i]) + partes_solido[j][i-1]));
                    // printf("thread %d ---- %.1lf = %.1lf + (%.1lf * (%.1lf - (2* %.1lf) + %.1lf))\n", tID, solido2[i], partes_solido[j][i], fourier, partes_solido[j][i+1], partes_solido[j][i], partes_solido[j][i-1]);

                    // printf("p_solido[%d][%d] = %.1lf\n", j, i, partes_solido[j][i]);
                    // printf("solido2[%d] = %.1lf\n\n", p+i-1, solido2[p+i-1]);
                }
            }

            // #pragma omp single
            // {
            //     // printIt(solido2, dimensao);
            //     // printItN(partes_solido, dimensao, tamanho_partes, tamanho_extra, cores);
            // }

            valor = 0;
            
            // #pragma omp for reduction(+: valor)
            for (int p = 1; p <= dimensao; p += tamanho_partes)
            {
                j = (int)((p-1)/tamanho_partes);
                med = (p <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
                for (int i = 1; i <= med; ++i)
                {
                    valor += abs(solido2[p+i-1] - partes_solido[j][i]);
                }
            }
            
            Z = valor / (dimensao);

            // printf("Z = %lf, valor = %lf\n\n", Z, valor);

            #pragma omp master
            printf("Z = %lf > 0\n\n", Z);

            if (Z <= 0)
            {
                #pragma omp master
                printf("Z < 0\n");
                on = 0;
            }
            else
            {
                #pragma omp for
                for (int i = 1; i <= dimensao; i += tamanho_partes)
                {
                    med = (i <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
                    for (int j = 0; j < med+2; ++j)
                    {
                        partes_solido[(int)((i-1)/tamanho_partes)][j] = solido2[i+j-1];
                    }
                }
            }
        }
    }   
    free(solido);
    free(solido2);

    ticks[1] = clock();

    double tempoTomado = (double)(ticks[1] - ticks[0]) / (double)(CLOCKS_PER_SEC); 
    printf("Tempo = %.10lf segundos\n", tempoTomado);
}