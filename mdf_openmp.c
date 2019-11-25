#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>
#include <time.h>

// #include "fila.h"

// #ifdef _WIN32
// #include <Windows.h>
// #else
// #include <unistd.h>
// #endif

double*** criaMatriz(int x, int y, int z, double value);

double*** alocarMatriz(int x, int y, int z);

void desalocarMatriz(double ****matriz, int x, int y, int z);

void copiarMatriz(double ****cola, double ***copia, int x, int y, int z);

double*** modTempPlane(double ***m, int x, int y, int z, int pos, double temp);

void calorMDF(double h, double tempo, int dimensao, double alfa);

void print_matriz(double ***m, int x, int y, int z);

void calorMDFThread(double h, double tempo, int dimensao, double alfa);

void calculoDoPonto(int **coords);


int main(int argc, char const *argv[])
{
    double h = 1;
    double tempo, ap = 0;
    int as;
    double coef_cond = 1.0;
    int dimensao = 10;

    printf("dimensao:\n");

    scanf("%d", &dimensao);

    while (1)
    {
        printf("Precisao:\n"
            "(1) - baixa\n"
            "(2) - media\n"
            "(3) - alta\n");

        // scanf("%d", &as);

        as = 3;

        if (as == 1)
        {
            ap = 2.0;
            break;
        }
        else if (as == 2)
        {
            ap = 4.0;
            break;
        }
        else if (as == 3)
        {
            ap = 6.0;
            break;
        }
        else
            printf("resposta inesperada\n");
    }

    // tempo = pow(h, 2)/ap*coef_cond+1;

    // while (tempo > pow(h, 2)/ap*coef_cond)
    // {
    //     printf("Digite um valor para tempo, tal que 'tempo < %lf'\n", pow(h, 2)/ap*coef_cond);
    //     scanf("%lf", &tempo);
    // }
    tempo = 0.1;

    printf("Considere o valor da condutividade termica em cm^2/s\n");

    calorMDF(h, tempo, dimensao, coef_cond);
    return 0;
}

double*** criaMatriz(int x, int y, int z, double value)
{
    double ***matriz = alocarMatriz(x, y, z); //Aloca um Vetor de Ponteiros

    for (int i = 0; i < x; i++)
    {
        // matriz[i] = (double**) malloc(y * sizeof(double*));
        for (int j = 0; j < y; j++)
        {
            // matriz[i][j] = (double*) malloc(z * sizeof(double));
            for (int k = 0; k < z; k++)
            {
                matriz[i][j][k] = value;
            }
        }
    }

    return matriz;
}

double*** alocarMatriz(int x, int y, int z)
{
    double ***matriz = (double***) malloc(x * sizeof(double**)); //Aloca um Vetor de Ponteiros

    for (int i = 0; i < x; i++)
    {
        matriz[i] = (double**) malloc(y * sizeof(double*));
        for (int j = 0; j < y; j++)
        {
            matriz[i][j] = (double*) malloc(z * sizeof(double));
        }
    }

    return matriz;
}

void desalocarMatriz(double ****matriz, int x, int y, int z)
{
    for (int i = 0; i < x; i++)
    {
        for (int j = 0; j < y; j++)
        {
            free((*matriz)[i][j]);
        }
        free((*matriz)[i]);
    }
    free((*matriz));
}

void copiarMatriz(double ****cola, double ***copia, int x, int y, int z)
{
    for (int i = 0; i < x; ++i)
    {
        for (int j = 0; j < y; ++j)
        {
            for (int k = 0; k < z; ++k)
            {
                (*cola)[i][j][k] = copia[i][j][k];
            }
        }
    }
}

double*** modTempPlane(double ***m, int x, int y, int z, int pos, double temp)
{
    int t = y/2;
    int nj, nk;

    if (t % 2 == 0)
        nj = t/2;
    else
        nj = (t-1)/2;

    t = z / 2;

    if (t % 2 == 0)
        nk = t / 2;
    else
        nk = (t - 1) / 2;

    for (int j = 0; j < y; j++)
    {
        for (int k = 0; k < z; k++)
        {
            if (((nj-1 < j) && (j < z-nj)) && ((nk-1 < k) && (k < z-nk)))
                m[pos][j][k] = temp;
        }
    }

    return m;
}

void calorMDF(double h, double tempo, int dimensao, double alfa)
{
    int d2 = dimensao+2,
    dim;

    double temp = 35,
    mod_temp = 0;

    int on = 1;

    double fourier,
    c_vizinhas_som,
    f;

    double Z = 0.0,
    valor = 0.0;

    int tID,
    nT;

    // printf("Digite a temperatura base do cubo:\n");
    // scanf("%lf", &temp);

    // printf("Digite a temperatura que sera aplicada a uma area do cubo:\n");
    // scanf("%lf", &mod_temp);
    mod_temp = 0;

	clock_t ticks[2];
    ticks[0] = clock();

    double ***solido = modTempPlane(criaMatriz(d2, d2, d2, temp), d2, d2, d2, 0, mod_temp);
    double ***solido2 = alocarMatriz(d2, d2, d2);
    copiarMatriz(&solido2, solido, d2, d2, d2);


    fourier = pow(alfa, 2) * (tempo/pow(h, 2));

    #pragma omp parallel default(none) private(tID, nT, on, c_vizinhas_som, dim, d2, f, Z) shared(fourier, dimensao, solido, solido2, valor) num_threads(4)
    {
        #pragma omp critical
        {
            dim = dimensao;
            f = fourier;
        }

        d2 = dim+2;
        on = 1;
        tID = omp_get_thread_num();
        nT = omp_get_num_threads();

        printf("thread %d de %d.\n", tID+1, nT);
        #pragma omp barrier

        while (on)
        {
            valor = 0;

            #pragma omp for schedule(auto)
            for (int i = 1; i <= dim; i++)
            {
                for (int j = 1; j <= dim; j++)
                {
                    for (int k = 1; k <= dim; k++)
                    {
                        c_vizinhas_som = (solido[i][j+1][k]
                                       + solido[i][j-1][k]
                                       + solido[i-1][j][k]
                                       + solido[i+1][j][k]
                                       + solido[i][j][k+1]
                                       + solido[i][j][k-1]
                        );
                        #pragma omp critical
                        {
                            solido2[i][j][k] = solido[i][j][k] + f * (c_vizinhas_som - (6 * solido[i][j][k])); // Equação
                            valor += abs(solido[i][j][k] - solido2[i][j][k]);
                        }
                    }
                }
            }

            // #pragma omp single
            // print_matriz(solido2, dim, dim, dim);

            // #pragma omp for reduction(+: valor)
            // #pragma omp single
            // for (int i = 1; i <= dim; i++)
            //     for (int j = 1; j <= dim; j++)
            //         for (int k = 1; k <= dim; k++)
            //         {
            //             #pragma omp atomic
            //         }
            // #pragma omp barrier

            Z = valor / ((d2) * (d2) * (d2));

            #pragma omp single
            printf("Z = %lf\n\n", Z);

            if (Z <= 0)
            {
                #pragma omp master
                printf("(Z < 0)\n");
                on = 0;
            }
            else
            {
                #pragma omp single
                {
                    for (int i = 0; i < dimensao; ++i)
                    {
                        copiarMatriz(&solido, solido2, d2, d2, d2);
                    }
                }
            }
        }
    }

    ticks[1] = clock();

    double tempoTomado = (double)(ticks[1] - ticks[0]) / (double)(CLOCKS_PER_SEC); 
    printf("Tempo = %.10lf segundos\n", tempoTomado);

    desalocarMatriz(&solido, d2, d2, d2);
    desalocarMatriz(&solido2, d2, d2, d2);
}

void print_matriz(double ***m, int x, int y, int z)
{
    for (int i = 1; i <= x; i++)
    {
        for (int j = 1; j <= y; j++)
        {
            for (int k = 1; k <= z; k++)
                // printf("M(%d)(%d)(%d) = %lf  ", i, j, k, m[i][j][k]);
                printf("%lf  ", m[i][j][k]);
            printf("\n");
        }
        printf("\n\n");
    }
}


// #include "mpi.h"
// #include <stdio.h>
// #include <stdlib.h>
// extern void draw_heat(int nx, int ny);       /* X routine to create graph */

// #define NXPROB      20                 /* x dimension of problem grid */
// #define NYPROB      20                 /* y dimension of problem grid */
// #define STEPS       100                /* number of time steps */
// #define MAXWORKER   8                  /* maximum number of worker tasks */
// #define MINWORKER   3                  /* minimum number of worker tasks */
// #define BEGIN       1                  /* message tag */
// #define LTAG        2                  /* message tag */
// #define RTAG        3                  /* message tag */
// #define NONE        0                  /* indicates no neighbor */
// #define DONE        4                  /* message tag */
// #define MASTER      0                  /* taskid of first process */

// struct Parms { 
//   float cx;
//   float cy;
// } parms = {0.1, 0.1};

// int main (int argc, char *argv[])
// {
// void inidat(), prtdat(), update();
// float  u[2][NXPROB][NYPROB];        /* array for grid */
// int taskid,                     /* this task's unique id */
//     numworkers,                 /* number of worker processes */
//     numtasks,                   /* number of tasks */
//     averow,rows,offset,extra,   /* for sending rows of data */
//     dest, source,               /* to - from for message send-receive */
//     left,right,        /* neighbor tasks */
//     msgtype,                    /* for message types */
//     rc,start,end,               /* misc */
//     i,ix,iy,iz,it;              /* loop variables */
// MPI_Status status;


// /* First, find out my taskid and how many tasks are running */
//    MPI_Init(&argc,&argv);
//    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
//    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
//    numworkers = numtasks-1;

//    if (taskid == MASTER) {
//       /************************* master code *******************************/
//       /* Check if numworkers is within range - quit if not */
//       if ((numworkers > MAXWORKER) || (numworkers < MINWORKER)) {
//          printf("ERROR: the number of tasks must be between %d and %d.\n",
//                  MINWORKER+1,MAXWORKER+1);
//          printf("Quitting...\n");
//          MPI_Abort(MPI_COMM_WORLD, rc);
//          exit(1);
//          }
//       printf ("Starting mpi_heat2D with %d worker tasks.\n", numworkers);

//       /* Initialize grid */
//       printf("Grid size: X= %d  Y= %d  Time steps= %d\n",NXPROB,NYPROB,STEPS);
//       printf("Initializing grid and writing initial.dat file...\n");
//       inidat(NXPROB, NYPROB, u);
//       prtdat(NXPROB, NYPROB, u, "initial.dat");

//       /* Distribute work to workers.  Must first figure out how many rows to */
//       /* send and what to do with extra rows.  */
//       averow = NXPROB/numworkers;
//       extra = NXPROB%numworkers;
//       offset = 0;
//       for (i=1; i<=numworkers; i++)
//       {
//          rows = (i <= extra) ? averow+1 : averow; 
//          /* Tell each worker who its neighbors are, since they must exchange */
//          /* data with each other. */  
//          if (i == 1) 
//             left = NONE;
//          else
//             left = i - 1;
//          if (i == numworkers)
//             right = NONE;
//          else
//             right = i + 1;
//          /*  Now send startup information to each worker  */
//          dest = i;
//          MPI_Send(&offset, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
//          MPI_Send(&rows, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
//          MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
//          MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
//          MPI_Send(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, dest, BEGIN, 
//                   MPI_COMM_WORLD);
//          printf("Sent to task %d: rows= %d offset= %d ",dest,rows,offset);
//          printf("left= %d right= %d\n",left,right);
//          offset = offset + rows;
//       }
//       /* Now wait for results from all worker tasks */
//       for (i=1; i<=numworkers; i++)
//       {
//          source = i;
//          msgtype = DONE;
//          MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, 
//                   &status);
//          MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//          MPI_Recv(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, source,
//                   msgtype, MPI_COMM_WORLD, &status);
//       }

//       /* Write final output, call X graph and finalize MPI */
//       printf("Writing final.dat file and generating graph...\n");
//       prtdat(NXPROB, NYPROB, &u[0][0][0], "final.dat");
//       printf("Click on MORE button to view initial/final states.\n");
//       printf("Click on EXIT button to quit program.\n");
//       draw_heat(NXPROB,NYPROB);
//       MPI_Finalize();
//    }   /* End of master code */



//    /************************* workers code **********************************/
//    if (taskid != MASTER) 
//    {
//       /* Initialize everything - including the borders - to zero */
//       for (iz=0; iz<2; iz++)
//          for (ix=0; ix<NXPROB; ix++) 
//             for (iy=0; iy<NYPROB; iy++) 
//                u[iz][ix][iy] = 0.0;

//       /* Receive my offset, rows, neighbors and grid partition from master */
//       source = MASTER;
//       msgtype = BEGIN;
//       MPI_Recv(&offset, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//       MPI_Recv(&rows, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//       MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//       MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
//       MPI_Recv(&u[0][offset][0], rows*NYPROB, MPI_FLOAT, source, msgtype, 
//                MPI_COMM_WORLD, &status);

//       /* Determine border elements.  Need to consider first and last columns. */
//       /* Obviously, row 0 can't exchange with row 0-1.  Likewise, the last */
//       /* row can't exchange with last+1.  */
//       start=offset;
//       end=offset+rows-1;
//       if (offset==0) 
//          start=1;
//       if ((offset+rows)==NXPROB) 
//          end--;
//       printf("task=%d  start=%d  end=%d\n",taskid,start,end);

//       /* Begin doing STEPS iterations.  Must communicate border rows with */
//       /* neighbors.  If I have the first or last grid row, then I only need */
//       /*  to  communicate with one neighbor  */
//       printf("Task %d received work. Beginning time steps...\n",taskid);
//       iz = 0;
//       for (it = 1; it <= STEPS; it++)
//       {
//          if (left != NONE)
//          {
//             MPI_Send(&u[iz][offset][0], NYPROB, MPI_FLOAT, left,
//                      RTAG, MPI_COMM_WORLD);
//             source = left;
//             msgtype = LTAG;
//             MPI_Recv(&u[iz][offset-1][0], NYPROB, MPI_FLOAT, source,
//                       msgtype, MPI_COMM_WORLD, &status);
//          }
//          if (right != NONE)
//          {
//             MPI_Send(&u[iz][offset+rows-1][0], NYPROB, MPI_FLOAT, right,
//                       LTAG, MPI_COMM_WORLD);
//             source = right;
//             msgtype = RTAG;
//             MPI_Recv(&u[iz][offset+rows][0], NYPROB, MPI_FLOAT, source, msgtype,
//                       MPI_COMM_WORLD, &status);
//          }
//          /* Now call update to update the value of grid points */
//          update(start,end,NYPROB,&u[iz][0][0],&u[1-iz][0][0]);
//          iz = 1 - iz;
//       }
//       /* Finally, send my portion of final results back to master */
//       MPI_Send(&offset, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
//       MPI_Send(&rows, 1, MPI_INT, MASTER, DONE, MPI_COMM_WORLD);
//       MPI_Send(&u[iz][offset][0], rows*NYPROB, MPI_FLOAT, MASTER, DONE, 
//                MPI_COMM_WORLD);
//       MPI_Finalize();
//    }
// }


// /**************************************************************************
//  *  subroutine update
//  ****************************************************************************/
// void update(int start, int end, int ny, float *u1, float *u2)
// {
//    int ix, iy;
//    for (ix = start; ix <= end; ix++) 
//       for (iy = 1; iy <= ny-2; iy++) 
//          *(u2+ix*ny+iy) = *(u1+ix*ny+iy)  + 
//                           parms.cx * (*(u1+(ix+1)*ny+iy) +
//                           *(u1+(ix-1)*ny+iy) - 
//                           2.0 * *(u1+ix*ny+iy)) +
//                           parms.cy * (*(u1+ix*ny+iy+1) +
//                          *(u1+ix*ny+iy-1) - 
//                           2.0 * *(u1+ix*ny+iy));
// }

// /*****************************************************************************
//  *  subroutine inidat
//  *****************************************************************************/
// void inidat(int nx, int ny, float *u) {
// int ix, iy;

// for (ix = 0; ix <= nx-1; ix++) 
//   for (iy = 0; iy <= ny-1; iy++)
//      *(u+ix*ny+iy) = (float)(ix * (nx - ix - 1) * iy * (ny - iy - 1));
// }

// /**************************************************************************
//  * subroutine prtdat
//  **************************************************************************/
// void prtdat(int nx, int ny, float *u1, char *fnam) {
// int ix, iy;
// FILE *fp;

// fp = fopen(fnam, "w");
// for (iy = ny-1; iy >= 0; iy--) {
//   for (ix = 0; ix <= nx-1; ix++) {
//     fprintf(fp, "%8.1f", *(u1+ix*ny+iy));
//     if (ix != nx-1) 
//       fprintf(fp, " ");
//     else
//       fprintf(fp, "\n");
//     }
//   }
// fclose(fp);
// }