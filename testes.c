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

int main(int argc, char const *argv[])
{
    int on,
    tID,
    nT;

    int dimensao = 5;

    double valor = 0,
    Z,
    fourier;

    double *solido = (double*) malloc(dimensao * sizeof(double));
    double *solido2 = (double*) malloc(dimensao * sizeof(double));

    for (int i = 0; i < dimensao; ++i)
    {
        solido[i] = 10;
        solido2[i] = 10;
    }
    solido[0] = 80;
    solido2[0] = 80;

    #pragma omp parallel default(none) private(tID, nT, on, dimensao, fourier) shared(solido, solido2) num_threads(3)
    {
        fourier = pow(1.0, 2) * (0.1/pow(1, 2));
        dimensao = 5;

        on = 1;
        tID = omp_get_thread_num();
        nT = omp_get_num_threads();

        printf("thread %d de %d.\n", tID+1, nT);

        #pragma omp barrier

        while(on)
        {
            #pragma omp for
            for (int i = 1; i < dimensao-1; ++i)
            {
                printf("thread %d no for %d da iteracao %d.\n", tID+1, i, on);
                // #pragma omp critical
                // solido2[i] = solido[i] + (fourier * (solido[i+1] - (2 * solido[i]) + solido[i-1]));
            }

            #pragma omp single
            {
                printf("somente a thread %d no bloco single da iteracao %d.\n", tID+1, on);
            }
            if (on < 3)
            {
                on++;
            }
            else
                on = 0;
        }
    }

    free(solido);
    free(solido2);
    return 0;
}


// int main(int argc, char const *argv[])
// {
//     int on, tID, nT;
//     int dimensao = 10;
//     double valor = 0, Z;

//     double *solido = (double*) malloc(dimensao * sizeof(double));
//     double *solido2 = (double*) malloc(dimensao * sizeof(double));

//     for (int i = 0; i < dimensao; ++i)
//     {
//         solido[i] = 10;
//         solido2[i] = 10;
//     }
//     solido[0] = 80;
//     solido2[0] = 80;

//     double fourier = pow(1.0, 2) * (0.1/pow(1, 2));


//     #pragma omp parallel default(none) private(tID, nT, on, dimensao, fourier) shared(solido, solido2) num_threads(3)
//     {
//         dimensao = 10;
//         fourier = pow(1.0, 2) * (0.1/pow(1, 2));

//         on = 1;
//         tID = omp_get_thread_num();
//         nT = omp_get_num_threads();
//         printf("thread %d de %d.\n", tID+1, nT);

//         #pragma omp barrier

//         while (on)
//         {
//             #pragma omp for
//             for (int i = 1; i < dimensao-1; i++)
//             {
//                 printf("thread %d no for %d da iteracao %d.\n", tID+1, i, on);

//                 #pragma omp critical
//                 solido2[i] = solido[i] + (fourier * (solido[i+1] - (2 * solido[i]) + solido[i-1]));
//             }

//             #pragma omp single
//             {
//                 printIt(solido2, dimensao);
//                 printIt(solido, dimensao);
//             }

//             // valor = 0;
            
//             // // #pragma omp for reduction(+: valor)
//             // for (int i = 1; i < dimensao; i++)
//             // {
//             //     valor += abs(solido2[i] - solido[i]);
//             // }

//             // Z = valor / (dimensao);

//             // printf("Z = %lf\n\n", Z);

//             // if (Z <= 0)
//             // {
//             //     printf("%lf < 0\n", Z);
//             //     on = 0;
//             // }

//             if (on < 3)
//             {
//                 on++;
//             }
//             else
//             {
//                 on = 0;
//             }

//             #pragma omp single
//             {
//                 for (int i = 0; i < dimensao; ++i)
//                 {
//                     solido[i] = solido2[i];
//                 }
//             }
//         }
//     }   

//     free(solido);
//     free(solido2);

// }