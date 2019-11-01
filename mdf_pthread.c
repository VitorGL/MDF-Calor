#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>
#include <time.h>
#include "fila.h"

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

double fourier;
double ***solido;
double ***solido2;
List *fila = NULL;
int dimensao;
pthread_mutex_t w_lock;


int main(int argc, char const *argv[])
{
    double h = 1;
    double tempo, ap = 0;
    int as;
    double coef_cond = 1.0;
    dimensao = 10;

    while (1)
    {
        printf("Precisao:\n"
            "(1) - baixa\n"
            "(2) - media\n"
            "(3) - alta\n");

        scanf("%d", &as);

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

    tempo = pow(h, 2)/ap*coef_cond+1;

    while (tempo > pow(h, 2)/ap*coef_cond)
    {
        printf("Digite um valor para tempo, tal que 'tempo < %lf'\n", pow(h, 2)/ap*coef_cond);
        scanf("%lf", &tempo);
    }

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
            if (nj-1 < j < z-nj && nk-1 < k < z-nk)
                m[pos][j][k] = temp;
        }
    }

    return m;
}

void calorMDF(double h, double tempo, int dimensao, double alfa)
{
    int d2 = dimensao+2;
    double temp = 35, mod_temp = 0;

    // printf("Digite a temperatura base do cubo:\n");
    // scanf("%lf", &temp);

    printf("Digite a temperatura que sera aplicada a uma area do cubo:\n");
    scanf("%lf", &mod_temp);

    double ***u = modTempPlane(criaMatriz(d2, d2, d2, 35), d2, d2, d2, 0, mod_temp);
    double ***u2 = alocarMatriz(d2, d2, d2);
    copiarMatriz(&u2, u, d2, d2, d2);

    int on = 1;
    double fourier, c_vizinhas_som;
    double Z = 0, Z2 = 0, valor = 0;

    fourier = pow(alfa, 2) * (tempo/pow(h, 2));

    while (on)
    {
        copiarMatriz(&u2, u, d2, d2, d2);

        for (int i = 1; i <= dimensao; i++)
        {
            for (int j = 1; j <= dimensao; j++)
            {
                for (int k = 1; k <= dimensao; k++)
                {
                    if (0 < i <= dimensao && 0 < j <= dimensao && 0 < k <= dimensao)
                    {
                        c_vizinhas_som = (u[i][j+1][k]
                                       + u[i][j-1][k]
                                       + u[i-1][j][k]
                                       + u[i+1][j][k]
                                       + u[i][j][k+1]
                                       + u[i][j][k-1]
                        );

                        u2[i][j][k] = u[i][j][k] + fourier * (c_vizinhas_som - (6 * u[i][j][k])); // Equação
                        // printf("%lf + (%lf/%lf * %lf) * (%lf - %lf)\n",u[i][j][k], tempo, pow(h, 2), alfa, c_vizinhas_som, 6*u[i][j][k]);
                        // printf("u2[%d][%d][%d] = %lf\n", i, j, k, u2[i][j][k]);
                        // sleep(0.5);
                    }
                }
            }
        }

        print_matriz(u2, dimensao, dimensao, dimensao);

        valor = 0;

        for (int i = 1; i <= dimensao; i++)
            for (int j = 1; j <= dimensao; j++)
                for (int k = 1; k <= dimensao; k++)
                {
                    valor += abs(u[i][j][k] - u2[i][j][k]);
                }


        Z = valor / ((d2) * (d2) * (d2));

        printf("Z = %lf\n\n", Z);
        // sleep(1);

        if (Z <= 0) //  || Z == Z2
        {
            printf("A condicaoo de estabilidade foi atingida com valor (%lf < 0)\n", Z);
            on = 0;
        }
        else
            Z2 = Z;

        copiarMatriz(&u, u2, d2, d2, d2);
    }

    FILE *arq;

    if ((arq = fopen("log.txt", "w")) == NULL);
        printf("Erro ao abrir arquivo de log\n");

    for (int i = 1; i <= dimensao; i++)
    {
        for (int j = 1; j <= dimensao; j++)
        {
            for (int k = 1; k <= dimensao; k++)
                fprintf(arq, "%lf  ", u[i][j][k]);
            fprintf(arq, "\n");
        }
        fprintf(arq, "\n\n");
    }

    fclose(arq);

    desalocarMatriz(&u, d2, d2, d2);
    desalocarMatriz(&u2, d2, d2, d2);
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

void calorMDFThread(double h, double tempo, int dimensao, double alfa)
{
    // printf("Digite a temperatura base do cubo:\n");
    // scanf("%lf", &temp);

    // printf("Digite a temperatura que sera aplicada a uma area do cubo:\n");
    // scanf("%lf", &mod_temp);
    double mod_temp = 20.0;

    clock_t inicio = clock();

    int d2 = dimensao+2;
    double temp = 35;
    int on = 1;
    double Z = 0, Z2 = 0, valor = 0;

    solido = modTempPlane(criaMatriz(d2, d2, d2, 35), d2, d2, d2, 0, mod_temp);
    solido2 = alocarMatriz(d2, d2, d2);
    copiarMatriz(&solido2, solido, d2, d2, d2);

    fourier = pow(alfa, 2) * (tempo/pow(h, 2));

    clock_t ini_thread = clock();

    // global q
    // fila = Queue()

    pthread_t tid[8];
    for (int i = 0; i < 8; i++)
    {
        pthread_create(&(tid[i]), NULL, threader, (void *)&tid); // ultimo parametro é um ponteiro void pro argumento da funcao
    }

    int qntx = 2;
    int qnty = 2;
    int qntz = 2;
    int qntIntervalos = 0;

    int subx = dimensao / qntx;
    int suby = dimensao / qnty;
    int subz = dimensao / qntz;

    int cx;
    int fx;
    int cy;
    int fy;
    int cz;
    int fz;

    int ***intervalos = (int***) malloc((qntx*qnty*qntz) * sizeof(int**));

    for (int i = 0; i < qntx; i++)
    {
        cx = (i * subx)+1;
        fx = ((i * subx) + subx);
        // print(cx, fx)
        for (int j = 0; j < qnty; j++)
        {
            cy = (j * suby)+1;
            fy = ((j * suby) + suby);
            // print(cy, fy)
            for (int k = 0; k < qntz; k++)
            {
                cz = (k * subz)+1;
                fz = ((k * subz) + subz);
                // print(cz, fz)

                intervalos[qntIntervalos] = (int**) malloc(3 * sizeof(int*));
                intervalos[qntIntervalos][0] = (int*) malloc(2 * sizeof(int));
                intervalos[qntIntervalos][1] = (int*) malloc(2 * sizeof(int));
                intervalos[qntIntervalos][2] = (int*) malloc(2 * sizeof(int));
                intervalos[qntIntervalos][0][0] = cx;
                intervalos[qntIntervalos][0][1] = fx;
                intervalos[qntIntervalos][1][0] = cy;
                intervalos[qntIntervalos][1][1] = fy;
                intervalos[qntIntervalos][2][0] = cz;
                intervalos[qntIntervalos][2][1] = fz;
                qntIntervalos++;
            }
        }
    }

    double tempo_thread = clock() - ini_thread;
    
    while (on)
    {
        copiarMatriz(&solido2, solido, d2, d2, d2);

        for (int i = 0; i < (qntx*qnty*qntz); i++)
        {
            insert(&fila, (void *)(intervalos[i]));
        }

        q.join(); // usar o semaforo lock aqui
        // print_matriz(solido2);

        valor = 0;


        for (int i = 1; i <= dimensao; i++)
            for (int j = 1; j <= dimensao; j++)
                for (int k = 1; k <= dimensao; k++)
                {
                    valor += abs(solido[i][j][k] - solido2[i][j][k]);
                    // printf("valor = %lf = |%lf - %lf|\n\n", valor, u[i][j][k], u2[i][j][k]);
                }


        Z = valor / ((d2) * (d2) * (d2));

        printf("Z = %lf\n\n", Z);

        if (Z <= 0) //  || Z == Z2
        {
            printf("A condicaoo de estabilidade foi atingida com valor (%lf < 0)\n", Z);
            on = 0;
        }
        else
            Z2 = Z;

        copiarMatriz(&solido, solido2, d2, d2, d2);


    }

    printf("Tempo= %li", (clock() - inicio));
    printf("Tempo de inicialização a mais das threads= %li", tempo_thread);

    print_matriz(solido2, dimensao, dimensao, dimensao);

    desalocarMatriz(&solido, d2, d2, d2);
    desalocarMatriz(&solido2, d2, d2, d2);
}

void calculoDoPonto(int **coords)
{
    int cx = coords[0][0], fx = coords[0][1];
    int cy = coords[1][0], fy = coords[1][1];
    int cz = coords[2][0], fz = coords[2][1];
    double c_vizinhas_som;

    // Calculo da função
    for (int i = 1; i <= dimensao; i++)
    {
        for (int j = 1; j <= dimensao; j++)
        {
            for (int k = 1; k <= dimensao; k++)
            {
                if (0 < i <= dimensao && 0 < j <= dimensao && 0 < k <= dimensao)
                {
                    c_vizinhas_som = (solido[i][j+1][k]
                                   + solido[i][j-1][k]
                                   + solido[i-1][j][k]
                                   + solido[i+1][j][k]
                                   + solido[i][j][k+1]
                                   + solido[i][j][k-1]
                    );

                    pthread_mutex_lock(&w_lock);
                    solido2[i][j][k] = solido[i][j][k] + fourier * (c_vizinhas_som - (6 * solido[i][j][k])); // Equação
                    pthread_mutex_unlock(&w_lock);
                }
            }
        }
    }
}


void threader()
{
    int **coords;
    while (1)
    {
        coords = (int *)pop(&fila);
        calculoDoPonto(coords);
        q.task_done(); // usar o semaforo unlock aqui
    }
}
