#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include <unistd.h>
#include "fila.h"

List *L = NULL;

void *threader(void *nsei);

int main(int argc, char const *argv[])
{

	int ***intervalos = (int***) malloc(sizeof(int**));
	intervalos[0] = (int**) malloc(sizeof(int*));
    intervalos[0][0] = (int*) malloc(2 * sizeof(int));
    intervalos[0][0][0] = 1;
    intervalos[0][0][1] = 2;
	insert(&L, (void *)(&(intervalos[0][0])));
	int **coisa = (int *)pop(&L);

    printf("%d  %d\n", coisa[0][0], coisa[0][1]);

	// exib(&L);


	// pthread_t tid[8];
	// for (int i = 0; i < 8; i++)
	// {
	// 	pthread_create(&(tid[i]), NULL, threader, (void *)&i);
	// }

	// for (int i = 0; i < 8; ++i)
	// {
	// 	pthread_join(tid[i], NULL);
	// }

	// for (int i = 0; i < 8; ++i)
	// {
	// 	pthread_join(tid[i], NULL);
	// }

	return 0;
}

void calculoDoPonto(int coisas)
{
	printf("%d\n", coisas);
}

void *threader(void *pargs)
{
	int coords = 0;
	while (1)
	{
        calculoDoPonto(coords);
	}
}