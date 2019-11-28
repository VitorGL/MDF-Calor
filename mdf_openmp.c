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

double*** alocarSolido(int x, int y, int z);

double**** alocarPartesSolido(int y, int z, int partes, int cores, int tamanho_partes, int tamanho_extra);

void desalocarSolido(double ****matriz, int x, int y, int z);

void copiarMatriz(double ****cola, double ***copia, int x, int y, int z);

double*** modTempPlane(double ***m, int x, int y, int z, int pos, double temp);

void calorMDF(double h, double tempo, int dimensao, double alfa);

void print_matriz(double ***m, int x, int y, int z);

void calorMDFThread(double h, double tempo, int dimensao, double alfa);

void calculoDoPonto(int **coords);


int main(int argc, char const *argv[])
{
    double h = 1;
    double tempo;
    // double ap = 0;
    // int as;
    double coef_cond = 1.0;
    int dimensao = atoi(argv[1]);


    // printf("dimensao:\n");
    // scanf("%d", &dimensao);

    // while (1)
    // {
    //     printf("Precisao:\n"
    //         "(1) - baixa\n"
    //         "(2) - media\n"
    //         "(3) - alta\n");

    //     scanf("%d", &as);

    //     if (as == 1)
    //     {
    //         ap = 2.0;
    //         break;
    //     }
    //     else if (as == 2)
    //     {
    //         ap = 4.0;
    //         break;
    //     }
    //     else if (as == 3)
    //     {
    //         ap = 6.0;
    //         break;
    //     }
    //     else
    //         printf("resposta inesperada\n");
    // }

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
    double ***matriz = alocarSolido(x, y, z); //Aloca um Vetor de Ponteiros

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

double*** alocarSolido(int x, int y, int z)
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

double**** alocarPartesSolido(int y, int z, int partes, int cores, int tamanho_partes, int tamanho_extra)
{
    double ****partes_solido = (double****) malloc(partes * sizeof(double***));
    for (int p = 0; p < partes; ++p)
    {
        if (p == cores)
        {
            partes_solido[p] = (double***) malloc((tamanho_extra+2) * sizeof(double**));
            for (int i = 0; i < tamanho_extra+2; i++)
            {
                partes_solido[p][i] = (double**) malloc(y * sizeof(double*));
                for (int j = 0; j < y; j++)
                {
                    partes_solido[p][i][j] = (double*) malloc(z * sizeof(double));
                    // printf("solido[%d][%d][%d] = alloc(%d*)\n", p, i, j, z);
                }
            }
            printf("Tamanho %d para a parte %d.\n", tamanho_extra+2, p);
        }
        else
        {
            partes_solido[p] = (double***) malloc((tamanho_partes+2) * sizeof(double**));
            for (int i = 0; i < tamanho_partes+2; i++)
            {
                partes_solido[p][i] = (double**) malloc(y * sizeof(double*));
                for (int j = 0; j < y; j++)
                {
                    partes_solido[p][i][j] = (double*) malloc(z * sizeof(double));
                    // printf("solido[%d][%d][%d] = alloc(%d*)\n", p, i, j, z);
                }
            }
            printf("Tamanho %d para a parte %d.\n", tamanho_partes+2, p);
        }
    }

    return partes_solido;
}

void desalocarSolido(double ****matriz, int x, int y, int z)
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
    double ***solido2 = alocarSolido(d2, d2, d2);
    double ****partes_solido;
    copiarMatriz(&solido2, solido, d2, d2, d2);

    int cores = 8;

    int tamanho_partes = (cores <= dimensao) ? dimensao/cores : 1;
    int tamanho_extra = (cores <= dimensao) ? dimensao%cores : 0;
    int partes = (cores <= dimensao) ? ((0 < tamanho_extra) ? cores+1 : cores) : dimensao;
    printf("Cores presentes:%d\ntamanho das partes:%d\ntamanho da parte extra:%d\n", cores, tamanho_partes, tamanho_extra);

    partes_solido = alocarPartesSolido(d2, d2, partes, cores, tamanho_partes, tamanho_extra);

    fourier = pow(alfa, 2) * (tempo/pow(h, 2));

    #pragma omp parallel default(none) private(tID, nT, on, c_vizinhas_som, dim, d2, f, Z) shared(cores, partes, partes_solido, tamanho_partes, tamanho_extra, fourier, dimensao, solido, solido2, valor) num_threads(cores)
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
        int m;
        int med;

        printf("thread %d de %d.\n", tID+1, nT);
        
        #pragma omp for
        for (int p = 1; p < dimensao; p += tamanho_partes)
        {
            med = (p <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
            for (int i = 0; i < med+2; ++i)
            {
                for (int j = 0; j < d2; ++j)
                {
                    for (int k = 0; k < d2; ++k)
                    {
                        // printf("solido[%d][%d][%d][%d]\n", (int)((p-1)/tamanho_partes), i, j, k);
                        partes_solido[(int)((p-1)/tamanho_partes)][i][j][k] = solido[p+i-1][j][k];
                        // printf("solido[%d][%d][%d][%d] = %lf\n", (int)((p-1)/tamanho_partes), i, j, k, partes_solido[(int)((p-1)/tamanho_partes)][i][j][k]);
                    }
                }
            }
        }

        while (on)
        {
            valor = 0;

            #pragma omp for
            for (int p = 1; p < dimensao; p += tamanho_partes)
            {
                m = (int)((p-1)/tamanho_partes);
                med = (p <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
                for (int i = 1; i <= med; ++i)
                {
                    for (int j = 1; j <= dim; j++)
                    {
                        for (int k = 1; k <= dim; k++)
                        {
                            // printf("solido[%d][%d][%d][%d]\n", (int)((p-1)/tamanho_partes), i, j, k);
                            // printf("solido[%d][%d][%d][%d] = %lf\n", (int)((p-1)/tamanho_partes), i, j, k, partes_solido[(int)((p-1)/tamanho_partes)][i][j][k]);
                            c_vizinhas_som = (partes_solido[m][i][j+1][k]
                                           + partes_solido[m][i][j-1][k]
                                           + partes_solido[m][i-1][j][k]
                                           + partes_solido[m][i+1][j][k]
                                           + partes_solido[m][i][j][k+1]
                                           + partes_solido[m][i][j][k-1]
                            );
                            solido2[p+i-1][j][k] = partes_solido[m][i][j][k] + f *(c_vizinhas_som - (6 * partes_solido[m][i][j][k])); // Equação
                        }
                    }
                }
            }

            // #pragma omp single
            // print_matriz(solido2, dim, dim, dim);

            #pragma omp for reduction(+: valor)
            for (int p = 1; p < dimensao; p += tamanho_partes)
            {
                m = (int)((p-1)/tamanho_partes);
                med = (p <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
                for (int i = 1; i <= med; ++i)
                {
                    for (int j = 1; j <= dim; j++)
                    {
                        for (int k = 1; k <= dim; k++)
                        {
                            valor += abs(solido2[p+i-1][j][k] - partes_solido[m][i][j][k]);
                        }
                    }
                }
            }

            Z = valor / ((d2) * (d2) * (d2));

            #pragma omp single
            printf("Z = %lf\n\n", Z);

            if (Z <= 0)
            {
                #pragma omp single
                printf("thread %d -> (Z < 0)\n", tID);
                on = 0;
            }
            else
            {
                #pragma omp for
                for (int p = 1; p < dimensao; p += tamanho_partes)
                {
                    med = (p <= (cores*tamanho_partes)) ? tamanho_partes : tamanho_extra;
                    for (int i = 0; i < med+2; ++i)
                    {
                        for (int j = 0; j < d2; ++j)
                        {
                            for (int k = 0; k < d2; ++k)
                            {
                                partes_solido[(int)((p-1)/tamanho_partes)][i][j][k] = solido2[p+i-1][j][k];
                            }
                        }
                    }
                }
            }
        }
    }

    ticks[1] = clock();

    double tempoTomado = (double)(ticks[1] - ticks[0]) / (double)(CLOCKS_PER_SEC); 
    printf("Tempo = %.10lf segundos\n", tempoTomado);

    desalocarSolido(&solido, d2, d2, d2);
    desalocarSolido(&solido2, d2, d2, d2);
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