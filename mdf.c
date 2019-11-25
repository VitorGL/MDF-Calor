#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// #include "fila.h"

#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif

double*** criaMatriz(int x, int y, int z, double value);

double*** alocarMatriz(int x, int y, int z);

void desalocarMatriz(double ****matriz, int x, int y, int z);

void copiarMatriz(double ****cola, double ***copia, int x, int y, int z);

double*** modTempPlane(double ***m, int x, int y, int z, int pos, double temp);

void calorMDF(double h, double tempo, int dimensao, double alfa);

void print_matriz(double ***m, int x, int y, int z);

void calorMDFThread(double h, double tempo, int dimensao, double alfa);

void calculoDoPonto(int **coords);

void threader();

int dimensao;


int main(int argc, char const *argv[])
{
    double h = 1;
    double tempo, ap = 0;
    int as;
    double coef_cond = 1.0;
    // dimensao = 20;

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
    int d2 = dimensao+2;

    double temp = 35,
    mod_temp = 0;

    int on = 1;

    double fourier,
    c_vizinhas_som;

    double Z = 0,
    valor = 0;

    // printf("Digite a temperatura base do cubo:\n");
    // scanf("%lf", &temp);

    // printf("Digite a temperatura que sera aplicada a uma area do cubo:\n");
    // scanf("%lf", &mod_temp);

    clock_t ticks[2];
    ticks[0] = clock();

    double ***solido = modTempPlane(criaMatriz(d2, d2, d2, 35), d2, d2, d2, 0, mod_temp);
    double ***solido2 = alocarMatriz(d2, d2, d2);
    copiarMatriz(&solido2, solido, d2, d2, d2);


    fourier = pow(alfa, 2) * (tempo/pow(h, 2));

    while (on)
    {
        valor = 0;
        for (int i = 1; i <= dimensao; i++)
        {
            for (int j = 1; j <= dimensao; j++)
            {
                for (int k = 1; k <= dimensao; k++)
                {
                    c_vizinhas_som = (solido[i][j+1][k]
                                   + solido[i][j-1][k]
                                   + solido[i-1][j][k]
                                   + solido[i+1][j][k]
                                   + solido[i][j][k+1]
                                   + solido[i][j][k-1]
                    );

                    solido2[i][j][k] = solido[i][j][k] + fourier * (c_vizinhas_som - (6 * solido[i][j][k])); // Equação   
                    valor += abs(solido[i][j][k] - solido2[i][j][k]);
                }
            }
        }
       
        // print_matriz(solido2, dimensao, dimensao, dimensao);


        Z = valor / ((d2) * (d2) * (d2));

        printf("Z = %lf\n\n", Z);

        if (Z <= 0)
        {
            printf("A condicaoo de estabilidade foi atingida com valor (%lf < 0)\n", Z);
            on = 0;
        }
        else
            copiarMatriz(&solido, solido2, d2, d2, d2);
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




// FILE *arq;

//     if ((arq = fopen("log.txt", "w")) == NULL);
//         printf("Erro ao abrir arquivo de log\n");

//     for (int i = 1; i <= dimensao; i++)
//     {
//         for (int j = 1; j <= dimensao; j++)
//         {
//             for (int k = 1; k <= dimensao; k++)
//                 fprintf(arq, "%lf  ", u[i][j][k]);
//             fprintf(arq, "\n");
//         }
//         fprintf(arq, "\n\n");
//     }

//     fclose(arq);

//     desalocarMatriz(&u, d2, d2, d2);
//     desalocarMatriz(&u2, d2, d2, d2);