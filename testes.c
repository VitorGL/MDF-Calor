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


// int main(int argc, char const *argv[])
// {
//     int vetor[10];

//     for (int i = 0; i < 10; ++i)
//     {
//         vetor[i] = 5+i;
//     }

//     int i, j;
//     int me, n;

//     #pragma omp parallel private(i,j,me,n)
//     {
//         me = omp_get_thread_num();
//         n = omp_get_num_threads();
//         printf( "Hello from %d/%d\n", me, n );
//         #pragma omp barrier

//         for ( i = 0; i < 10; i++ )
//         {
//             #pragma omp for
//             for ( j = 0; j < 5; j++ )
//             { 
//                 printf( "thread(%d) processing %d,%d\n", me, i, j );
//                 // printf("thread() - vetor[%d] = %d\n", i, vetor[i]);

//                 #pragma omp atomic
//                 vetor[i] += 5;

//                 // usleep( 500000 );
//             }
//             #pragma omp barrier
//         }
//     }
//     return 0;
// }


int main(int argc, char const *argv[])
{
    int tID, nT;

    int on = 1;

    while(on)
    {
        #pragma omp parallel default(none) private(tID, nT) shared(on) num_threads(3)
        {
            tID = omp_get_thread_num();
            nT = omp_get_num_threads();

            printf("thread(%d) - while[%d]\n", tID, on);
            #pragma omp barrier

            // exit(0);
            #pragma omp single
            {
                on++;
                if (on == 4)
                {
                    on = 0;
                    
                }
            }
            #pragma omp barrier
        }
    }
    return 0;
}

// int main(int argc, char const *argv[])
// {
//     int on = 0, tID;
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


//     #pragma omp parallel private(on, tID)
//     {
//         tID = omp_get_thread_num();
//         int n = omp_get_num_threads();
//         printf( "Hello from %d/%d\n", tID, n );

//         #pragma omp barrier
//         for (on = 0; on < 100; on++)
//         {
//             #pragma omp for
//             for (int i = 1; i < dimensao; i++)
//             {
//                 if (0 < i && i < dimensao-1)
//                 {
//                     printf( "thread(%d) processing %d,%d\n", tID, i, on);

//                     #pragma omp critical
//                     solido2[i] = solido[i] + (fourier * (solido[i+1] - (2 * solido[i]) + solido[i-1]));
//                 }
//             }
//             #pragma omp barrier

//             #pragma omp single
//             {
//                 printIt(solido2, dimensao);
//                 printIt(solido, dimensao);

//                 valor = 0;
                
//                 // #pragma omp for reduction(+: valor)
//                 for (int i = 1; i < dimensao; i++)
//                 {
//                     valor += abs(solido2[i] - solido[i]);
//                 }

//                 Z = valor / (dimensao);

//                 printf("Z = %lf\n\n", Z);

//                 if (Z <= 0)
//                 {
//                     printf("%lf < 0\n", Z);
//                     on = 100;
//                 }
//                 else
//                     on--;

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